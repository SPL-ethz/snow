clear all; close all; clc;

%% Heat conduction model for freezing of packed vials 

tic
%% Next steps: Old
%(1) Make sure that the heat flows are computed based on the temperature readouts of the prior time step. This is currently not the case and leads to some errors regarding voxels that should freeze at the same time
%(2) Add stochastic nucleation temperature distribution: Use an arbitrary
%distribution to compute nucleation propabilities for specific temperatures
%and dt and a random number generator to decide if nucleation occurs. 

%1 is done, but one dimension still behaves anisotropic -> issue fixed
%2 is done via Braatz approach: empirical nucleation relation with
%parameters b and kb

%% Next steps: Update
%Make a plot that shows the total time to total freezing of all vials for
%all cycles (all data already stored in time_nuc, just need to extract

%% Next steps (Dec 14): Include entire process, i.e. crystal formation and sublimation time
%For crystal size assume that it scales with freezing time only
%Use 1D sublimation model to track drying process. (idea: vials with larger
%crystals sublimate faster, thus cool down stronger. consequently they
%receive heat from neighboring vials that sublimate then even slower.)
%After fast vials completely sublimated they heat up thus enhancing
%sublimation in neighboring vials

%Done

%% Next steps: January 20: Improve sublimation model 
%step 1: consider real mass balance with dried fraction to avoid numerical
%issues -> Part is done
%step 2: integrate position-dependent dP and Rp, doing so requires saving
%this info during freezing -> extraction of dP done, but not yet integrated
%in sublimation model

%% Next steps: January 30: Improve sublimation model
%step1: use non-simplified energy balance, i.e. add dm/dt part 
%step2: reevaluate Rp calculation

%% Next steps: March 17: Upscaling and Vectorizing
%Eriks model is capable of running with 50k vials in 10s. He has vectorized
%his model, so I do this too. Currently runtime is 23 min for two cycles.

%% Next steps: Further speedup
%I can try removing the two innermost loops by using vial as variable.
%Problem there is that I have to change all the output information, too.

%% Material and geometrical properties

%Geometry and materials
lambda = 0.8*0.05; % W/[mK] thermal conductivity of glass / times 0.1 to take into account that the actual contact area is smaller
thickness = 0.002; % [m] thickness of glass (2 mm, usually 2 layers of 1 mm)
len = 0.01; % [m] side length of cubic vial volume
A = len*len;
V = A*len;

%Sample properties
rho_l = 1000; %kg/m^3 density of water / assumed constant for all phases
k_f = 1.853; %K kg / mol cryoscopic constant
cp_w = 4187; %J/[Kkg] heat capacity liquid water
cp_i = 2108; %J/[Kkg] heat capacity ice 0°C
cp_s = 1240; %J/Kkg, heat capacity sucrose
solid_fraction = 0.05; %mass fraction of solute in solution
Dcp = cp_i - cp_w;
cp_solution = solid_fraction*cp_s + (1-solid_fraction)*cp_w; %heat capacity of solution
M_s = 342.3/1000; %molar mass sucrose
T_eq = 0; %°C equilibrium freezing temperature
T_nuc = -10;%°C nucleation temperature
Dh = 333550; %J/kg heat of fusion water 

%Nucleation parameters (Braatz)
kb = 0.000000001; %m−3 s−1 K−b value reduced compared to Braatz paper (he had a value of 10)
b = 12;

%Sublimation parameters (sublimation is not part of the model anymore)
C_emp = 0.1*10^-6; % Empirical parameter for Kurz and Fisher model
t_factor = 0.225; %tortuosity factor tau^2/epsilon
Mw = 0.01802; %kg/mol molar mass of water
Rg = 8.3144; %gas constant
%T_sub0 = -45 + 273.15;% Initial temperature of sublimating vials
T_sub0 = 225;% Initial temperature of sublimating vials
T_shelf = -10 + 273.15;
pwc = 13; %Chamber pressure in Pa
solid_fraction = 0.05;
drho = rho_l*(1-solid_fraction); %density change during sublimation
Dhsub = 2.83*10^6; % latent heat of sublimation [J/kg]
Antoine = [4.6543	1435.26 -64.848];

as(1) = -0.212144006*100;
as(2) = 0.273203819*100;
as(3) = -0.610598130*10;

bs(1) = 0.333333333/100;
bs(2) = 0.120666667*10;
bs(3) = 0.170333333*10;

pcrit = 611.657; %[Pa] Critical pressure
Tcrit = 273.16; %[K] Critical temperature

%HT Calculations
kc_xy = 10; %heat transfer in xy plane
kc_z = 10; %heat transfer in z direction
h = 10; %external heat transfer coefficient

mass = rho_l*V;
hl = mass.*cp_solution;
hs = mass.*cp_i;
depression = k_f/M_s *(solid_fraction/(1-solid_fraction));

% sigma_equivalent = Dh/cp_l; % temperature rise required for complete solidification
% sigma_solid = sigma_equivalent.*cp_l./cp_s; %similar temperature rise for a frozen sample, i.e. with solid heat capacity

%% Define model parameters

% Number of vials in each dimension (box: 20/12/3, pallet: 40/36/18)
x = 20;
y = 12;
z = 3; 

%number of cycles
N = 100; 

%Cooling rate
CR_min = 0.5; %K/min
CR = CR_min/60;

Dt = 2; %s timestep
n = 20000; % number of timesteps to carry out 
t_tot = Dt*n/60; %total time in min


t_time = linspace(Dt,Dt*n,n);

parpool(16);

%% Model for heat transfer during packed freezing
% Asumptions: Vials are cubic, experience conductive HT with immediate
% neighbors only, nucleation occurs when vials reach a certain temperature,
% nucleation is followed by recalescence to T_eq, stay at T_eq till
% completely frozen (assume pure compound), no temperature gradient within
% vials

q = zeros(x*y*z,1);
%Define interaction matrix

Hext = A.*h.*Dt;
Hint_xy = A.*kc_xy.*Dt;
Hint_z = A.*kc_z.*Dt;

%self interaction

d0 = -6*ones(x*y*z,1);

%x interaction (left and right)

dx = ones(x*y*z-1,1);

for u = 1:(x*y*z-1)
    if mod(u,x) == 0
        dx(u) = 0;
    end
end

%y interaction (up and down)

dy = ones(x*y*z-x,1);

for u = 1:(x*y*z-x-1)
    if mod(u,(y*x)) == 0
        for k = 1:x
        dy(u+1-k) = 0;
        end
        u = u+x;
    end
end

%z interaction (front and back)

dz = ones(x*y*(z-1),1);

%Final interaction matrix

INT = sparse(diag(d0) + diag(dx,1) + diag(dx,-1) + diag(dy,x) + diag(dy,-x) + diag(dz,x*y) + diag(dz,-x*y));
%Boundaries

EXT = - sum(INT); %number of contact faces with environment of each vial

%Heat transfer matrix

d_red = d0' + EXT; %remove external interactions

%Calculation of self contribution requires separation of external HT, xy HT
%and z HT
for vial = 1:(x*y*z)
    if vial < (x*y+1) || vial > (length(d_red) - x*y)
    d_red(vial) = d_red(vial) +1; 
    Z_vial = 1;
    else
    d_red(vial) = d_red(vial) + 2;
    Z_vial = 2;
    end
    H0(vial) = d_red(vial)*Hint_xy - Z_vial*Hint_z - EXT(vial)*Hext;
    
end

Hx = dx.*Hint_xy;
Hy = dy.*Hint_xy;
Hz = dz.*Hint_z;

HEXT = EXT.*Hext;
HINT = sparse(diag(H0) + diag(Hx,1) + diag(Hx,-1) + diag(Hy,x) + diag(Hy,-x) + diag(Hz,x*y) + diag(Hz,-x*y));

kbDtV = kb*Dt*V;


xyz = zeros(x,y,z);
xyzp = zeros(x,y,z,N);

time_frozen = xyzp;
time_nuc = xyzp;

time_frozen_vial = zeros( (x*y*z),N);
time_nuc_vial = zeros( (x*y*z),N);
temp_nuc_vial = zeros( (x*y*z),N);


Time_nuc = zeros(N,4);
Time_frozen = zeros(N,4);

time_max_nuc = zeros(N,1);
time_min_nuc = zeros(N,1);
time_max_frozen = zeros(N,1);
time_min_frozen = zeros(N,1);

%%
parfor p=1:N
    
rng(p) 

%% Define temperatures within the loop to be faster
T_ext = -40; %Temperature of the cooling medium (°C)
T_init = 20; %Initial temperature of the vials (°C)

T_time = (T_init*ones(x*y*z,1)); %Time_dependent temperature in each vial
sigma_time = zeros(x*y*z,1); %Time_dependent ice fraction in each vial
sigma_init = zeros(x*y*z,1); %initial ice fraction at recalescence



%Define temporary variables to ease data transfer     
time_frozen_temp = zeros(x*y*z,1);
time_nuc_temp = zeros(x*y*z,1);
temp_nuc_temp = zeros(x*y*z,1);

%Output for temperature
temp_vial_corner = zeros(n,1);
sigma_vial_corner = zeros(n,1);
sigma_vial_center = zeros(n,1);
temp_vial_edge = zeros(n,1);
temp_vial_center = zeros(n,1);

solid_counter = zeros(x*y*z,1);

for i=1:n
    
    q_ext = HEXT*T_ext;
    q = HINT*T_time + q_ext';
  
    
    temp_vial_corner(i) = T_time(1);
    sigma_vial_corner(i) = sigma_time(1);
    sigma_vial_center(i) = sigma_time(round(x*y*z/2+x*y/2+x/2));
    temp_vial_edge(i) = T_time(round(x/2));
    temp_vial_center(i) = T_time(round(x*y*z/2+x*y/2+x/2));

    for vial =1:(x*y*z)

               
            %% decision tree: Two states, liquid and ice growth
            %Note that there is currently no implementation for the reverse shift
            %from ice growth state to liquid.
            %Both controlled and stochastic nucleation can be used to shift
            %to ice growth state            
            
            %ice growth state
            if sigma_time(vial) > 0 
                beta = (1-solid_fraction);
                cp_sigma = solid_fraction*cp_s+ beta*(sigma_time(vial)*Dcp + cp_w);
                alpha = -q(vial)/mass;               
                gamma = -cp_sigma*depression;
                delta = 1 - sigma_time(vial);
                
              %omega = beta*delta*(Dh + T_time(vial)*Dcp);
              omega = beta*delta*Dh;                
                A1 =  -alpha*delta -(1+sigma_time(vial))*omega + gamma;
                A0 =  alpha*delta + sigma_time(vial)*(omega-gamma);
                DIS = sqrt(A1^2 - 4*omega*A0);
                sigma_time(vial) = (- A1 - DIS)/(2*omega);

                T_time(vial) = T_eq - depression*(1/delta);
                
                if solid_counter(vial) == 0
                     if sigma_time(vial) > 0.9 %identify completely frzoen vials                    
                        time_frozen_temp(vial) = i+1;
                        solid_counter(vial) = 1;
                     end
                end
           
            
            %liquid vials
            else
            T_time(vial) = T_time(vial) + (q(vial)./hl);  
            if T_time(vial) < T_eq %only test for nucleation if nucleation is thermodynamically feasible
             Xx = kbDtV.*(T_eq - T_time(vial)).^b; %probability of nucleation
             Yy = rand(); %random number for stochastic nucleation
             
            
                 if Yy < Xx %then we have nucleation
                time_nuc_temp(vial) = i+1;
                temp_nuc_temp(vial) = T_time(vial);
                q_0 = (T_eq - T_time(vial))*cp_solution*mass;                                       
                alpha = q_0/mass;
                beta = (1-solid_fraction);
                gamma = -cp_solution*k_f/M_s *(solid_fraction/(1-solid_fraction));
                delta = 1;
                
              %omega = beta*delta*(Dh + T_time(vial)*Dcp);
              omega = beta*delta*Dh;
                A2 = omega;
                A1 =-alpha*delta -(1+sigma_time(vial))*omega + gamma;
                A0 = alpha*delta + sigma_time(vial)*(omega-gamma);
                DIS = sqrt(A1^2 - 4*A2*A0);
                sigma_time(vial) = (- A1 - DIS)/(2*A2);
                T_time(vial) = T_eq - k_f/M_s *(solid_fraction/(1-solid_fraction))*(1/(1-sigma_time(vial)));
                   
                   
                 end
                 
                 
            end
            end
            
            
            
            
    end
end

time_frozen_vial(:,p) = time_frozen_temp;
time_nuc_vial(:,p) = time_nuc_temp;
temp_nuc_vial(:,p) = temp_nuc_temp;

temp_corner_all(:,p) = temp_vial_corner;
sigma_corner_all(:,p) = sigma_vial_corner;
sigma_center_all(:,p) = sigma_vial_center;
temp_edge_all(:,p) = temp_vial_edge;
temp_center_all(:,p) = temp_vial_center;

end
%% Data extraction from loop

%Change format of output variables
for vial = 1:(x*y*z)
    av = floor( (vial-1)/(x*y));
    bv = mod( (vial-1),(x*y));
    cv = floor(bv/x);
    dv = mod(bv,x);
    
    xv = dv + 1;
    yv = cv + 1;
    zv = av + 1;
    
    time_frozen(xv,yv,zv,:) = time_frozen_vial(vial,:);
    time_nuc(xv,yv,zv,:) = time_nuc_vial(vial,:);
    temp_nuc(xv,yv,zv,:) = temp_nuc_vial(vial,:);
    
end



delete(gcp('nocreate')); %close parallel computing
toc;

%%

time_solid = time_frozen - time_nuc;

Time_nuc_all = time_nuc.*Dt./60;
Time_frozen_all = time_frozen.*Dt./60;
Time_solid_all = Time_frozen_all - Time_nuc_all;
Temp_nuc_all = temp_nuc;

%%

for p=1:N

Time_nuc(p,1) = time_nuc(1,1,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,2) = time_nuc(2,1,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,3) = time_nuc(2,2,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,4) = time_nuc(2,2,2,p).*Dt./60; %[min] nucleation times in min

Time_frozen(p,1) = time_frozen(1,1,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,2) = time_frozen(2,1,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,3) = time_frozen(2,2,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,4) = time_frozen(2,2,2,p).*Dt./60; %[min] nucleation times in min


C = max(max(max(time_nuc(:,:,:,p))));
time_max_nuc(p) = C.*Dt./60;

D = min(min(min(time_nuc(:,:,:,p))));
time_min_nuc(p) = D.*Dt./60;

C_frozen = max(max(max(time_frozen(:,:,:,p))));
time_max_frozen(p) = C_frozen.*Dt./60;

D_frozen = min(min(min(time_frozen(:,:,:,p))));
time_min_frozen(p) = D_frozen.*Dt./60;

end

%% Plot temperature vs time 

times = linspace(0,Dt*n,n+1)./60;
iterations = linspace(1,N,N);

TimeN_sort = sort(Time_nuc,1,'ascend');
TimeF_sort = sort(Time_frozen,1,'ascend');
fraction = [0 linspace(0,1,(N+1)) 1];

for g=1:4
TimeN_sort_plot(g,:) = [0 TimeN_sort(1,g) TimeN_sort(:,g)' (10*TimeN_sort(end,g))];
TimeF_sort_plot(g,:) = [0 TimeF_sort(1,g) TimeF_sort(:,g)' (10*TimeF_sort(end,g))];
end
%%
% figure(1)
% box on;
% plot(times,Temp_corner(:,1))
% xlabel('Time [min]');
% ylabel('Temperature [°C]');
% %title('Temperature Evolution of Corner in 5 Cycles');
% axis([0 300 -20 21]);
% ax = gca
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;

%% 
figure(2)
plot(TimeN_sort_plot,fraction)
ylabel('Simulated CDF');
xlabel('Nucleation time [min]');
title('Simulation of Nucleation Times (2000 Cycles)');
legend('(1,1)','(1,4)','(2,2)','(4,4)');
%set(gca,'Xdir','reverse')
axis([0 175 -0.05 1.05]);


%% Find time when last vial nucleates

time_max_sort = sort(time_max_nuc,'ascend');
time_max_sort_plot = [0 time_max_sort(1) time_max_sort(:)' (10*time_max_sort(end))];

time_min_sort = sort(time_min_nuc,'ascend');
time_min_sort_plot = [0 time_min_sort(1) time_min_sort(:)' (10*time_min_sort(end))];

time_max_sort_frozen = sort(time_max_frozen,'ascend');
time_max_sort_plot_frozen = [0 time_max_sort_frozen(1) time_max_sort_frozen(:)' (10*time_max_sort_frozen(end))];

time_min_sort_frozen = sort(time_min_frozen,'ascend');
time_min_sort_plot_frozen = [0 time_min_sort_frozen(1) time_min_sort_frozen(:)' (10*time_min_sort_frozen(end))];

figure(3)
plot(time_max_sort_plot,fraction)
ylabel('Simulated CDF');
xlabel('Nucleation time [min]');
title('Nucleation time of last vial');
%set(gca,'Xdir','reverse')
axis([0 175 -0.05 1.05]);

figure(4)
plot(TimeN_sort_plot,fraction)
hold on;
plot(time_min_sort_plot,fraction)
hold on;
plot(time_max_sort_plot,fraction)
ylabel('Simulated CDF');
xlabel('Nucleation time [min]');
title('Nucleation Times (T$_{ext}$ = -30$^\circ$C)','interpreter','latex');
legend('Corner (1,1,1)','Edge (2,1,1)','Center (2,2,1)','Core (2,2,2)','First vial','Last vial');
%set(gca,'Xdir','reverse')
axis([0 900 -0.05 1.05]);
%%
figure(44)
plot(TimeF_sort_plot,fraction)
hold on;
plot(time_min_sort_plot_frozen,fraction)
hold on;
plot(time_max_sort_plot_frozen,fraction)
ylabel('Simulated CDF');
xlabel('Freezing time [min]');
%title('Simulation of Freezing Times (2000 Cycles)');
legend('Corner','Edge','Center','Core','First vial','Last vial');
%set(gca,'Xdir','reverse')
axis([0 375 -0.05 1.05]);
ax = gca
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
%%
% figure(5)
% hold on;
% plot(times,Temp(:,1:2))
% xlabel('Time [min]');
% ylabel('Temperature [°C]');
% title('Temperature Evolution of Corner Vial (1,1) in 5 Cycles');
% axis([0 150 -40 21]);

%% Nucleation times for first run

time_firstrun = time_nuc(:,:,:,1)./60.*Dt; %min

for a=1:x
    for d=1:y
        for c=1:z
            time_nuc_mean(a,d,c) = sum(time_nuc(a,d,c,:))/(60*N)*Dt;
            time_frozen_mean(a,d,c) = sum(time_frozen(a,d,c,:)-time_nuc(a,d,c,:))/(60*N)*Dt;
            time_nuc_std(a,d,c) = std(time_nuc(a,d,c,:))/60.*Dt;
            time_frozen_std(a,d,c) = std(time_frozen(a,d,c,:)-time_nuc(a,d,c,:))/60.*Dt;
        end
    end
end

%% Individual solidification times for a single batch:
time_frozen_1 = time_frozen(:,:,:,1)./60;
time_nuc_1 = time_nuc(:,:,:,1)./60;
temp_nuc_1 = temp_nuc(:,:,:,1);

time_solid_1 = (time_frozen(:,:,:,1) - time_nuc(:,:,:,1))./60.*Dt;

%% Differentiate temperature profile to identify nucleation events 
% The idea here is that we would see spikes at the times when the
% neighboring vials nucleate (which is indeed true)

% T1 = Temp(:,1); %temperature profile of core vial in first run
% T1_diff = diff(T1); 
% 
% figure(6)
% plot(times(2:end),T1_diff.*Dt)
% xlabel('Time [min]');
% ylabel('Temperature Change [K/s]');
% title('dT/dt Profile of Core Vial (2,2,2) in 1 Cycle');
% axis([0 150 -0.01 0.01]);


%% Freezing distributions for presentations
%creation of histograms to show how the nucleation/solid/freeze times vary

%histograms for all vials for all runs
%corner: 111 113 131 311 331 313 133 333
%edge: 211 121 112 213 123 231 132 321 312 332 323 233
%center: 122 221 212 232 223 322
%core: 222

%nucleation counts
bin_number = 500;
sol_number = 100;
temp_number = 300;

tempedges = linspace(-temp_number/20,0,temp_number+1);

binedges = linspace(0,bin_number,bin_number+1); %bin width of 1 min
CornerCounts = histcounts(time_nuc(1,1,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,1,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,3,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,1,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,3,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,1,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,3,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,3,3,:)./60.*Dt, 'BinEdge',binedges);
EdgeCounts = histcounts(time_nuc(2,1,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,2,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,1,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,1,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,2,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,3,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(1,3,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,2,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,1,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,3,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,2,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,3,3,:)./60.*Dt, 'BinEdge',binedges);
CenterCounts = histcounts(time_nuc(1,2,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,2,1,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,1,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,3,2,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(2,2,3,:)./60.*Dt, 'BinEdge',binedges)+histcounts(time_nuc(3,2,2,:)./60.*Dt, 'BinEdge',binedges);
CoreCounts = histcounts(time_nuc(2,2,2,:)./60.*Dt, 'BinEdge',binedges);


solidedges = linspace(0,sol_number,(2*sol_number+1));
CornerSolids = histcounts(time_solid(1,1,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,1,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,3,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,1,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,3,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,1,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,3,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,3,3,:)./60.*Dt, 'BinEdge',solidedges);
EdgeSolids = histcounts(time_solid(2,1,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,2,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,1,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,1,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,2,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,3,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(1,3,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,2,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,1,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,3,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,2,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,3,3,:)./60.*Dt, 'BinEdge',solidedges);
CenterSolids = histcounts(time_solid(1,2,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,2,1,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,1,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,3,2,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(2,2,3,:)./60.*Dt, 'BinEdge',solidedges)+histcounts(time_solid(3,2,2,:)./60.*Dt, 'BinEdge',solidedges);
CoreSolids = histcounts(time_solid(2,2,2,:)./60.*Dt, 'BinEdge',solidedges);

figure(11)
sgtitle('Position dependent nucleation time histograms (2000 cycles)')
subplot(2,2,1)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0)
xlabel('Nucleation time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
hold on
histogram('BinCounts',CornerCounts,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Corner vials')
grid
subplot(2,2,2)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0)
xlabel('Nucleation time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 250 0 2000]);
hold on
histogram('BinCounts',EdgeCounts,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Edge vials')
grid
subplot(2,2,3)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0)
xlabel('Nucleation time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 250 0 2000]);
hold on
histogram('BinCounts',CenterCounts,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Center vials')
grid
subplot(2,2,4)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0)
xlabel('Nucleation time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 250 0 2000]);
hold on
histogram('BinCounts',CoreCounts,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Core vial')
grid
%%
figure(110)
sgtitle('Position dependent solidification time histograms (2000 cycles)')
subplot(2,2,1)
histogram(time_solid./60.*Dt, 'BinWidth',0.5,'EdgeAlpha',0)
xlabel('Solidification time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
hold on
histogram('BinCounts',CornerSolids,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Corner vials')
axis([20 80 0 3000]);
grid
subplot(2,2,2)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0)
xlabel('Solidification time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 3000]);
hold on
histogram('BinCounts',EdgeSolids,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Edge vials')
grid
subplot(2,2,3)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0)
xlabel('Solidification time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 3000]);
hold on
histogram('BinCounts',CenterSolids,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Center vials')
grid
subplot(2,2,4)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0)
xlabel('Solidification time [min]');
ylabel('Counts');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 3000]);
hold on
histogram('BinCounts',CoreSolids,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Core vial')
grid
%%
%'Normalization','pdf'
figure(12)
histogram(time_solid./60.*Dt, 'BinWidth',.5)
xlabel('Solidification time [min]');
ylabel('Counts');
title('Histogram of all solidification times for 1000 cycles');
axis([0 80 0 3000]);
dim = [0.63 0.5 0.4 0.4];
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah12 = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah12.FontSize = 10;

%%
figure(13)
histogram(time_frozen./60.*Dt, 'BinWidth',1)
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 5000]);
title('Histogram of all freezing times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 0~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 0~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;

%%
% %histograms for all vials for first run: Not really interesting
% 
% figure(14)
% histogram(time_nuc_1, 50)
% xlabel('Nucleation time [min]');
% ylabel('Counts');
% title('Histogram of all nucleation times for a single cycle');
% 
% figure(15)
% histogram(time_solid_1, 50)
% xlabel('Solidification time [min]');
% ylabel('Counts');
% title('Histogram of all solidification times for a single cycle');
% 
% figure(16)
% histogram(time_frozen_1, 50)
% xlabel('Freezing time [min]');
% ylabel('Counts');
% title('Histogram of all freezing times for a single cycle');

%histograms for single corner vials for all runs
%%
figure(17)
histogram(time_nuc(1,1,1,:)./60.*Dt, 'BinWidth',1)
xlabel('Nucleation time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of corner vial (1,1) nucleation times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;
%%
figure(18)
histogram(time_solid(1,1,1,:)./60.*Dt,'BinWidth',1)
xlabel('Solidification time [min]');
ylabel('Counts');
axis([0 60 0 500]);
title('Histogram of corner vial (1,1) solidification times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;

figure(19)
histogram(time_frozen(1,1,1,:)./60.*Dt, 'BinWidth',1)
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of corner vial (1,1) freezing times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;

%histograms for single corner vials for all runs

figure(20)
histogram(time_nuc(2,2,2,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Nucleation time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of center vial (4,4) nucleation times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;

figure(21)
histogram(time_solid(2,2,2,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Solidification time [min]');
ylabel('Counts');
axis([0 60 0 500]);
title('Histogram of center vial (4,4) solidification times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;

figure(22)
histogram(time_frozen(2,2,2,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of center vial (4,4) freezing times for 1000 cycles');
str={'\textbf{HT Parameters}',...
    '$$k_{sh} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{int} = 20~ W m^{-2} K^{-1} $$'...
    '$$k_{ext} = 20~ W m^{-2} K^{-1} $$'};
ah = annotation('textbox',dim,'interpreter','latex','String',str,'FitBoxToText','on');

ah.FontSize = 10;


%% Freezing distributions independent of system size (needed for larger systems)

Counts_nuc_corner = zeros(bin_number,1)';
Counts_nuc_edge = zeros(bin_number,1)';
Counts_nuc_face = zeros(bin_number,1)';
Counts_nuc_core = zeros(bin_number,1)';

Counts_sol_corner = zeros(2*sol_number,1)';
Counts_sol_edge = zeros(2*sol_number,1)';
Counts_sol_face = zeros(2*sol_number,1)';
Counts_sol_core = zeros(2*sol_number,1)';

Counts_temp_corner = zeros(temp_number,1)';
Counts_temp_edge = zeros(temp_number,1)';
Counts_temp_face =  zeros(temp_number,1)';
Counts_temp_core =  zeros(temp_number,1)';

count_corner = 0;
count_edge = 0;
count_face = 0;
count_core = 0;

for k = 1:x
    for l = 1:y
        for m = 1:z
            counter = 0;
            if k == 1 || k == x
                counter = counter +1; 
            end
            if l == 1 || l == y
                counter = counter +1; 
            end
            if m == 1 || m == z
                counter = counter +1; 
            end

           if counter == 3
               Counts_nuc_corner = Counts_nuc_corner + histcounts(time_nuc(k,l,m,:)./60.*Dt, 'BinEdge',binedges);
               Counts_sol_corner = Counts_sol_corner + histcounts(time_solid(k,l,m,:)./60.*Dt, 'BinEdge',solidedges);
               Counts_temp_corner = Counts_temp_corner + histcounts(temp_nuc(k,l,m,:), 'BinEdge',tempedges);
               count_corner = count_corner +1;
           end
           if counter == 2
               Counts_nuc_edge = Counts_nuc_edge + histcounts(time_nuc(k,l,m,:)./60.*Dt, 'BinEdge',binedges);
               Counts_sol_edge = Counts_sol_edge + histcounts(time_solid(k,l,m,:)./60.*Dt, 'BinEdge',solidedges);
               Counts_temp_edge = Counts_temp_edge + histcounts(temp_nuc(k,l,m,:), 'BinEdge',tempedges);
               count_edge = count_edge +1;
           end
           if counter == 1
               Counts_nuc_face = Counts_nuc_face + histcounts(time_nuc(k,l,m,:)./60.*Dt, 'BinEdge',binedges);
               Counts_sol_face = Counts_sol_face + histcounts(time_solid(k,l,m,:)./60.*Dt, 'BinEdge',solidedges);
               Counts_temp_face = Counts_temp_face + histcounts(temp_nuc(k,l,m,:), 'BinEdge',tempedges);
               count_face = count_face +1;
           end
           if counter == 0
               Counts_nuc_core = Counts_nuc_core + histcounts(time_nuc(k,l,m,:)./60.*Dt, 'BinEdge',binedges);
               Counts_sol_core = Counts_sol_core + histcounts(time_solid(k,l,m,:)./60.*Dt, 'BinEdge',solidedges);
               Counts_temp_core = Counts_temp_core + histcounts(temp_nuc(k,l,m,:), 'BinEdge',tempedges);
               count_core = count_core +1;
           end  
           
        end
    end
end

%calculate fractions of population
frac_corner = count_corner / (x*y*z);
frac_edge = count_edge / (x*y*z);
frac_face = count_face / (x*y*z);
frac_core = count_core / (x*y*z);

divisor = x*y*z*N;
sol_length = solidedges(2) - solidedges(1);
temp_length = abs(tempedges(2) - tempedges(1));

Counts_nuc_corner = Counts_nuc_corner ./ divisor;
Counts_nuc_edge = Counts_nuc_edge ./ divisor;
Counts_nuc_face = Counts_nuc_face ./ divisor;
Counts_nuc_core = Counts_nuc_core ./ divisor;

Counts_sol_corner = Counts_sol_corner ./ (divisor*sol_length);
Counts_sol_edge = Counts_sol_edge ./ (divisor*sol_length);
Counts_sol_face = Counts_sol_face ./ (divisor*sol_length);
Counts_sol_core = Counts_sol_core ./ (divisor*sol_length);

Counts_temp_corner = Counts_temp_corner ./ (divisor*temp_length);
Counts_temp_edge = Counts_temp_edge ./ (divisor*temp_length);
Counts_temp_face = Counts_temp_face ./ (divisor*temp_length);
Counts_temp_core = Counts_temp_core ./ (divisor*temp_length);


figure(200)
sgtitle('Position dependent nucleation time histograms (2000 cycles)')
subplot(2,2,1)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0,'Normalization','probability')
xlabel('Nucleation time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 bin_number 0 .012]);
hold on
histogram('BinCounts',Counts_nuc_corner,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Corner vials')
grid
subplot(2,2,2)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0,'Normalization','probability')
xlabel('Nucleation time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 bin_number 0 .012]);
hold on
histogram('BinCounts',Counts_nuc_edge,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Edge vials')
grid
subplot(2,2,3)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0,'Normalization','probability')
xlabel('Nucleation time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 bin_number 0 .012]);
hold on
histogram('BinCounts',Counts_nuc_face,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Center vials')
grid
subplot(2,2,4)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'EdgeAlpha',0,'Normalization','probability')
xlabel('Nucleation time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([0 bin_number 0 .012]);
hold on
histogram('BinCounts',Counts_nuc_core,'BinEdge',binedges,'EdgeAlpha',0)
legend('All vials','Core vial')
grid






figure(210)
sgtitle('Position dependent solidification time histograms (2000 cycles)','interpreter','latex')
subplot(2,2,1)
histogram(time_solid./60.*Dt, 'BinWidth',0.5,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Solidification time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
hold on
histogram('BinCounts',Counts_sol_corner,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Corner vials')
axis([20 80 0 0.07]);
grid
subplot(2,2,2)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Solidification time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 0.07]);
hold on
histogram('BinCounts',Counts_sol_edge,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Edge vials')
grid
subplot(2,2,3)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Solidification time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 0.07]);
hold on
histogram('BinCounts',Counts_sol_face,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Center vials')
grid
subplot(2,2,4)
histogram(time_solid./60.*Dt, 'BinWidth',.5,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Solidification time [min]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([20 80 0 0.07]);
hold on
histogram('BinCounts',Counts_sol_core,'BinEdge',solidedges,'EdgeAlpha',0)
legend('All vials','Core vial')
grid

figure(220)
sgtitle('Position dependent nucleation temperature histograms (2000 cycles)','interpreter','latex')
subplot(2,2,1)
histogram(temp_nuc, 'BinWidth',0.05,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Nucleation temperature [$^{\circ}$C]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
hold on
histogram('BinCounts',Counts_temp_corner,'BinEdge',tempedges,'EdgeAlpha',0)
legend('All vials','Corner vials')
axis([-15 0 0 0.4]);
grid
subplot(2,2,2)
histogram(temp_nuc, 'BinWidth',0.05,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Nucleation temperature [$^{\circ}$C]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([-15 0 0 0.4]);
hold on
histogram('BinCounts',Counts_temp_edge,'BinEdge',tempedges,'EdgeAlpha',0)
legend('All vials','Edge vials')
grid
subplot(2,2,3)
histogram(temp_nuc, 'BinWidth',0.05,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Nucleation temperature [$^{\circ}$C]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([-15 0 0 0.4]);
hold on
histogram('BinCounts',Counts_temp_face,'BinEdge',tempedges,'EdgeAlpha',0)
legend('All vials','Center vials')
grid
subplot(2,2,4)
histogram(temp_nuc, 'BinWidth',0.05,'EdgeAlpha',0,'Normalization','pdf')
xlabel('Nucleation temperature [$^{\circ}$C]','interpreter','latex');
ylabel('Norm. distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([-15 0 0 0.4]);
hold on
histogram('BinCounts',Counts_temp_core,'BinEdge',tempedges,'EdgeAlpha',0)
legend('All vials','Core vial')
grid

%% Temperature and sigma plots

T_ext = ones(n,1)*(-40);

figure(102)
t = tiledlayout(1,2,'TileSpacing','Compact');
%Tile 1: 0
hn(1) = nexttile();
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_center_all(:,8:10),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('12x20x9 box: Thermal evolution of center vial','interpreter','latex','FontSize',14)
legend('Shelf','First simulation','Second simulation','Third simulation','interpreter','latex')
axis([0 5 -55 25])
grid on;

%Tile 2: 10
hn(2) = nexttile();
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_center_all(:,8:10),'LineWidth',1)
yline(0,'b--','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('(b)','interpreter','latex','FontSize',14)
legend('Shelf','First simulation','Second simulation','Third simulation','interpreter','latex')
axis([1.45 3 -15 5])
grid on;

title(t,'Thermally interacting vials: Stochasticity','interpreter','latex','FontSize',14)
ylabel(t,'Temperature [$$^{\circ}$$C]','FontSize',14,'interpreter','latex')
xlabel(t,'Process time [h]','FontSize',14,'interpreter','latex')

figure(103)
plot(t_time./3600,sigma_center_all(:,10))
axis([0 70 0.95 1])

temp_diff1 = diff(temp_center_all(:,1));
temp_diff2 = diff(diff(temp_center_all(:,1)));

figure(104)
plot(t_time(2:end)./3600,temp_diff1)
axis([2 2.7 -10e-3 10e-3])

figure(105)
plot(t_time(2:end-1)./3600,temp_diff2)
axis([2 2.7 -10e-6 10e-6])

figure(106)
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_corner_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_edge_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_center_all(:,1),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('Pallet freezing (40x36x18 vials): Thermal evolution ','interpreter','latex','FontSize',14)
legend('External temperature','Corner vial','Facial center vial','Core vial','interpreter','latex','FontSize',12)
xlabel('Process time [hr]','interpreter','latex','FontSize',14)
ylabel('Temperature [$$^{\circ}$$C]','interpreter','latex','FontSize',14)
axis([50 130 -41 21])
grid on;

%% Save data
% large box: 500k iterations correspond to 700hr, only 300-400 would have been required 

% temp_nuc_403618_101010_m40 = Temp_nuc_all;
% time_nuc_403618_101010_m40 = Time_nuc_all;
% time_solid_403618_101010_m40 = Time_solid_all;
% 
% temp_vial_cornerm40 = temp_corner_all(:,1:10);
% sigma_vial_cornerm40 = sigma_corner_all(:,1:10);
% sigma_vial_centerm40 = sigma_center_all(:,1:10);
% temp_vial_edgem40 = temp_edge_all(:,1:10);
% temp_vial_centerm40 = temp_center_all(:,1:10);
% 
% % save('3D_403618_xzh_101010_m20_dt5_500000_N1000_2.mat','temp_nuc_403618_101010_m20',...
% % 'time_nuc_403618_101010_m20','time_solid_403618_101010_m20')
% % 
% save('3D_403618_xzh_101010_m40_dt5_200000_N1000_2.mat','temp_nuc_403618_101010_m40',...
% 'time_nuc_403618_101010_m40','time_solid_403618_101010_m40','temp_vial_cornerm40','temp_vial_edgem40','sigma_vial_centerm40','sigma_vial_cornerm40','temp_vial_centerm40')

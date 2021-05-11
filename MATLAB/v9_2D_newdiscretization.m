clear all; close all; clc;

%% Heat conduction model for freezing of packed vials 
%Updated figure annotation scheme in figure18


%% V6: Able to operate with three different values of ht parameters
%currently operates with variable KS value for each iteration. since these
%values are generated within the parfor loop, the runtime is strongly
%affected. predefining kas, 5000 iteration = 29s. variable ks, 5000
%iteration = 143s
%Update: This was wrong. There was an error in the implementation of the internal
%KS variation. Wrong orientation of vector. With a ' it works now, the
%runtime now is comparable to the one with predefined KS. (April 18)

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

%% Next steps: February 17: Add variability of shelf heat transfer
%Use a standard distributed variability in shelf heat transfer coefficient.
%note that there are two ways to implement this. either inside or outside
%the p loop, i.e. the vials may be assigned random k values either in each
%loop or only once. Both of this may be interesting in fact (one may argue
%that the actual behavior here is deterministic)

%% February 22: Added parfor loop to increase speed 
%Introduced temporary variables for characteristic times. This allows
%plotting the CDF and histogram figures, but I lose information on
%temperature profiles. Need to add separately. Also I removed h_freeze

%% Version 7 (April 21)
%Added lots of minor improvements, cleaned up the code
%Completely parallelized and vectorized
%7c: capable of controlled nucleation



%% Material and geometrical properties

%Geometry and materials
lambda = 0.8*0.05; % W/[mK] thermal conductivity of glass / times 0.1 to take into account that the actual contact area is smaller
thickness = 0.002; % [m] thickness of glass (2 mm, usually 2 layers of 1 mm)
leng = 0.01; % [m] side length of cubic vial volume
A = leng*leng;
V = A*leng;

%Sample properties
rho_l = 1000; %kg/m^3 density of water / assumed constant for all phases
k_f = 1.853; %K kg / mol cryoscopic constant
cp_w = 4187; %J/[Kkg] heat capacity liquid water
cp_i = 2108; %J/[Kkg] heat capacity ice 0°C
cp_s = 1240; %J/Kkg, heat capacity sucrose
solid_fraction = 0.05; %mass fraction of solute in solution
Dcp = cp_i - cp_w;
cp_solution = solid_fraction*cp_s + (1-solid_fraction)*cp_w; %heat capacity of solution
T_init = 20; %°C initial temperature
T_ext0 = 20; %°C external temperature in the beginning
T_eq = 0; %°C equilibrium freezing temperature
T_nuc = -10;%°C nucleation temperature
Dh = 333550; %J/kg heat of fusion water 

%Nucleation parameters (Braatz)
kb = 0.000000001; %m−3 s−1 K−b value reduced compared to Braatz paper (he had a value of 10)
b = 12;

%Alternative nucleation parameters to simulate low T and broad
%distribution:
% kb = 0.000000001^0.5; %m−3 s−1 K−b value reduced compared to Braatz paper (he had a value of 10)
% b = 6;

%Sublimation parameters (sublimation is not part of the model anymore)
C_emp = 0.1*10^-6; % Empirical parameter for Kurz and Fisher model
t_factor = 0.225; %tortuosity factor tau^2/epsilon
Mw = 0.01802; %kg/mol molar mass of water
M_s = 342.3/1000; %molar mass sucrose
Rg = 8.3144; %gas constant
%T_sub0 = -45 + 273.15;% Initial temperature of sublimating vials
T_sub0 = 225;% Initial temperature of sublimating vials
T_shelf = -10 + 273.15;
pwc = 13; %Chamber pressure in Pa
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

%% Key HT parameters
kint = 1*20; %Heat transfer coefficient for thermal interaction
kext = 0*20; %external heat transfer coefficient
ks0 = 20;% W/m^2K heat transfer via shelf
ks_sigma = 0*ks0;

cont = 0; %Controlled nucleation: 1 = yes, everything else = no


%% Shortcut definitions
mass = rho_l*V;
hl = mass.*cp_solution;
hs = mass.*cp_i;
depression = k_f/M_s *(solid_fraction/(1-solid_fraction));

sigma_equivalent = Dh/cp_w; % temperature rise required for complete solidification
sigma_solid = sigma_equivalent.*cp_w./cp_i; %similar temperature rise for a frozen sample, i.e. with solid heat capacity

%% Define model parameters

% Number of vials in each dimension
x = 7;
y = x;
z = 1; %only one layer

%number of cycles
N = 5000; 

%Cooling rate
CR_min = 0.5; %K/min
CR = CR_min/60;

Dt = 2; %s timestep
n = 20000; % number of timesteps to carry out 
t_tot = Dt*n/60; %total time in min


t_time = linspace(Dt,Dt*n,n);
T_ext = T_ext0 - t_time.*CR;

%Consider a constant minimum temperature for the shelf 
%otherwise it would go down to -inf which is unreasonable
T_min = -50; %°C
T_ext(T_ext < T_min) = T_min;

%Include first holding step
h_time = 1/6; % duration of holding step in min 
h_length = h_time*60/Dt; %Not ideal here, since it must be an integer.
h_temp = -12; %holding temperature in °C
H_temp = ones(h_length,1).*h_temp;

T_insert = find(T_ext < h_temp,1);
T_insert_end = T_insert + h_length; %This is the timepoint of controlled nucleation, if activated
T_ext = [T_ext(1:T_insert-1) H_temp' T_ext(T_insert:end)];

n = length(T_ext);
t_time = linspace(Dt,Dt*n,n);


%% Include second holding step (for controlled nucleation)
%If no controlled nucleation, just comment this part of the code out

% h_time2 = 360; % duration of holding step in min 
% h_length2 = h_time2*60/Dt; %Not ideal here, since it must be an integer.
% h_temp2 = -10; %holding temperature in °C
% H_temp2 = ones(h_length2,1).*h_temp2;
% 
% T_insert2 = find(T_ext < h_temp2,1);
% T_ext = [T_ext(1:T_insert2-1) H_temp2' T_ext(T_insert2:end)];
% 
% n = length(T_ext);
% t_time = linspace(Dt,Dt*n,n);




%% Model for heat transfer during packed freezing
% Asumptions: Vials are cubic, experience conductive HT with immediate
% neighbors only, nucleation occurs when vials reach a certain temperature,
% nucleation is followed by recalescence to T_eq, stay at T_eq till
% completely frozen (assume pure compound), no temperature gradient within
% vials

kbDtV = kb*Dt*V;


parpool(16);

xyz = zeros(x,y,z);
xyzt = zeros(x,y,z,n+1);
%xyztp = zeros(x,y,z,n+1,N);
xyzp = zeros(x,y,z,N);



%h_freeze = xyztp; %initial height of freezed part

time_frozen = xyzp;
time_nuc = xyzp;


Time_nuc = zeros(N,4);
Time_frozen = zeros(N,4);

time_max_nuc = zeros(N,1);
time_min_nuc = zeros(N,1);
time_max_frozen = zeros(N,1);
time_min_frozen = zeros(N,1);

tic

%% Vectorization

q = zeros(x*y*z,1);
%Define interaction matrix

d0 = -4*ones(x*y*z,1);

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

%Final matrix

INT = sparse(diag(d0) + diag(dx,1) + diag(dx,-1) + diag(dy,x) + diag(dy,-x)); 
% + diag(dz,x*y) + diag(dz,-x*y));
clear d0 dx dy dz;
%Boundaries

EXT = - sum(INT);

INT = INT + diag(EXT);

HTEXT = A.*kext.*Dt;
HTINT = A.*kint.*Dt;

HINT = sparse(INT.*HTINT); %Only works this way if kint = kext

HEXT = EXT.*HTEXT; %Add shelf heat transfer as external heat transfer


parfor p=1:N
    
rng(p) 

T_time = (T_init*ones(x*y*z,1));
sigma_time = zeros(x*y*z,1);
sigma_init = zeros(x*y*z,1); %initial ice fraction at recalescence 

%HT shortcuts
B = -mass*Dh*(1-solid_fraction);
cp_0 = solid_fraction*cp_s+ (1-solid_fraction)*cp_w;
nonsolid_fraction = 1 - solid_fraction;
depression_mass = depression*mass;

%Output characteristic quantities    
time_frozen_temp = zeros(x*y*z,1);
time_nuc_temp = zeros(x*y*z,1);
temp_nuc_temp = zeros(x*y*z,1);

%Output for holding
hold_temp_init = zeros(x*y*z,1);
hold_temp_final = zeros(x*y*z,1);
hold_count_temp = 0;
time_final_temp = 1;

%Output for temperature
temp_vial_corner = zeros(n,1);
sigma_vial_corner = zeros(n,1);
sigma_vial_center = zeros(n,1);
temp_vial_edge = zeros(n,1);
temp_vial_center = zeros(n,1);

T_i = T_ext;

KS = zeros(x*y*z,1);

for vial = 1:(x*y*z)
    KS(vial) = ks0 + ks_sigma*randn(1,1);                
end   

HS =  sparse(KS'.*A.*Dt);
HSEXT = (HEXT+HS)';
HEINT = HINT - diag(HEXT) - diag(HS);

counter = 0; %counts timesteps at hold temperature
solid_counter = zeros(x*y*z,1);

for i=1:n
    
    %Need to think about automating definition of vial positions 
    temp_vial_corner(i) = T_time(1);
    sigma_vial_corner(i) = sigma_time(1);
    sigma_vial_center(i) = sigma_time(25);
    temp_vial_edge(i) = T_time(4);
    temp_vial_center(i) = T_time(25);

    q = HEINT*T_time + HSEXT*T_i(i);
    
    for vial =1:(x*y*z)

              
            %% decision tree: Two states, liquid and ice growth
            %Note that there is currently no implementation for the reverse shift
            %from ice growth state to liquid.
            %Both controlled and stochastic nucleation can be used to shift
            %to ice growth state
            
            
            %ice growth state
            if sigma_time(vial) > 0 
               % cp_sigma = cp_0 + nonsolid_fraction*sigma_time(vial)*Dcp;                
                C = depression_mass*(cp_0 + nonsolid_fraction*sigma_time(vial)*Dcp);
               sigma_time(vial) = sigma_time(vial) + q(vial)/(B - C/((1-sigma_time(vial))^2));
                T_time(vial) = T_eq - depression.*(1/(1-sigma_time(vial)));
                
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
             
             %Condition for controlled nucleation 
             if cont == 1 %cont yes / no in parameter list
             if i == T_insert_end %Time at end of holding step
                 Yy = 0; %Ensure nucleation
             end
             end
            
                if Yy < Xx %then we have nucleation
                time_nuc_temp(vial) = i+1;
                temp_nuc_temp(vial) = T_time(vial);
                q_0 = -(T_eq - T_time(vial))*cp_solution*mass; %needs to be negative otherwise we get negative sigma
                CC = cp_solution*depression_mass;
                sigma_time(vial) = q_0./(B-CC);
             
                T_time(vial) = T_eq - depression*(1/(1-sigma_time(vial)));
                   
                    
                    if T_i(i) == h_temp
                        hold_count_temp = hold_count_temp + 1;
                    end
                 end
                 
                 
            end
            end
            
            
            
            
    end
    
    if T_i(i) == h_temp
        hold_temp_final = T_time;
        time_final_temp = i;
    end
    if T_i(i) > h_temp
        hold_temp_init = T_time;
    end
        
    
end

time_frozen_vial(:,p) = time_frozen_temp;
time_nuc_vial(:,p) = time_nuc_temp;
temp_nuc_vial(:,p) = temp_nuc_temp;
hold_count(p) = hold_count_temp;
hold_final_vial(:,p) = hold_temp_final;
hold_init_vial(:,p) =  hold_temp_init;
time_final(p) = time_final_temp;
KS_all(:,p) = KS;

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
    
    hold_init(xv,yv,zv,:) = hold_init_vial(vial,:);
    hold_final(xv,yv,zv,:) = hold_final_vial(vial,:);
    
end


delete(gcp('nocreate')); %close parallel computing

time_solid = time_frozen - time_nuc;

toc

Time_nuc_all = time_nuc.*Dt./60;
Time_frozen_all = time_frozen.*Dt./60;
Time_solid_all = Time_frozen_all - Time_nuc_all;
Temp_nuc_all = temp_nuc;

%% Data extraction from loop
for p=1:N

    
Time_nuc(p,1) = time_nuc(1,1,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,2) = time_nuc(1,4,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,3) = time_nuc(2,2,1,p).*Dt./60; %[min] nucleation times in min
Time_nuc(p,4) = time_nuc(4,4,1,p).*Dt./60; %[min] nucleation times in min

Time_frozen(p,1) = time_frozen(1,1,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,2) = time_frozen(1,4,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,3) = time_frozen(2,2,1,p).*Dt./60; %[min] nucleation times in min
Time_frozen(p,4) = time_frozen(4,4,1,p).*Dt./60; %[min] nucleation times in min
end

Time_solid = Time_frozen - Time_nuc;

for p=1:N
C = max(max(max(time_nuc(:,:,:,p))));
time_max_nuc(p) = C.*Dt./60;

D = min(min(min(time_nuc(:,:,:,p))));
time_min_nuc(p) = D.*Dt./60;

C_frozen = max(max(max(time_frozen(:,:,:,p))));
time_max_frozen(p) = C_frozen.*Dt./60;

D_frozen = min(min(min(time_frozen(:,:,:,p))));
time_min_frozen(p) = D_frozen.*Dt./60;    

C_solid = max(max(max(time_solid(:,:,:,p))));
time_max_solid(p) = C_solid.*Dt./60;

D_solid = min(min(min(time_solid(:,:,:,p))));
time_min_solid(p) = D_solid.*Dt./60;    
    
end





%% Plot temperature vs time 

times = linspace(0,Dt*n,n+1)./60;
iterations = linspace(1,N,N);

TimeN_sort = sort(Time_nuc,1,'ascend');
TimeF_sort = sort(Time_frozen,1,'ascend');
TimeS_sort = sort(Time_solid,1,'ascend');
fraction = [0 linspace(0,1,(N+1)) 1];

for g=1:4
TimeN_sort_plot(g,:) = [0 TimeN_sort(1,g) TimeN_sort(:,g)' (10*TimeN_sort(end,g))];
TimeF_sort_plot(g,:) = [0 TimeF_sort(1,g) TimeF_sort(:,g)' (10*TimeF_sort(end,g))];
TimeS_sort_plot(g,:) = [0 TimeS_sort(1,g) TimeS_sort(:,g)' (10*TimeS_sort(end,g))];

end

% figure(1)
% hold on;
% plot(times,Temp(:,1:5))
% xlabel('Time [min]');
% ylabel('Temperature [°C]');
% title('Temperature Evolution of Core Vial (2,2,2) in 5 Cycles');
% axis([0 150 -40 21]);
% 
% % 
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

time_max_sort_solid = sort(time_max_solid,'ascend');
time_max_sort_plot_solid = [0 time_max_sort_solid(1) time_max_sort_solid(:)' (10*time_max_sort_solid(end))];

time_min_sort_solid = sort(time_min_solid,'ascend');
time_min_sort_plot_solid = [0 time_min_sort_solid(1) time_min_sort_solid(:)' (10*time_min_sort_solid(end))];


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
title('Simulation of Nucleation Times (1000 Cycles)');
legend('(1,1)','(1,4)','(2,2)','(4,4)','First vial','Last vial');
%set(gca,'Xdir','reverse')
axis([0 185 -0.05 1.05]);
dim = [0.15 0.5 0.4 0.4];


figure(44)
plot(TimeF_sort_plot,fraction)
hold on;
plot(time_min_sort_plot_frozen,fraction)
hold on;
plot(time_max_sort_plot_frozen,fraction)
ylabel('Simulated CDF');
xlabel('Freezing time [min]');
title('Simulation of Freezing Times (1000 Cycles)');
legend('(1,1)','(1,4)','(2,2)','(4,4)','First vial','Last vial');
%set(gca,'Xdir','reverse')
axis([0 185 -0.05 1.05]);
dim = [0.15 0.5 0.4 0.4];

figure(444)
plot(TimeS_sort_plot,fraction)
hold on;
plot(time_min_sort_plot_solid,fraction)
hold on;
plot(time_max_sort_plot_solid,fraction)
ylabel('Simulated CDF');
xlabel('Solidification time [min]');
title('Simulation of Solidification Times (1000 Cycles)');
legend('(1,1)','(1,4)','(2,2)','(4,4)','First vial','Last vial');
%set(gca,'Xdir','reverse')
axis([0 185 -0.05 1.05]);
dim = [0.15 0.5 0.4 0.4];


% %%
% figure(5)
% hold on;
% plot(times,Temp_corner(:,1:5))
% xlabel('Time [min]');
% ylabel('Temperature [°C]');
% title('Temperature Evolution of Corner Vial (1,1) in 5 Cycles');
% axis([0 150 -40 21]);

%% Nucleation times for first run

time_firstrun = time_nuc(:,:,:,1)./60.*Dt; %min

time_nuc_mean = xyz;
time_nuc_std = xyz;

time_frozen_mean = xyz;
time_frozen_std = xyz;

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

figure(11)
histogram(time_nuc./60.*Dt, 'BinWidth',1,'Normalization','probability')
xlabel('Nucleation time [min]');
ylabel('Normalized Distribution [min$$^{-1}$$]','interpreter','latex');
%title('Histogram of all nucleation times for 1000 cycles');
axis([80 160 0 0.2]);
dim = [0.4 0.5 0.4 0.4];

%%
figure(12)
histogram(time_solid./60.*Dt, 'BinWidth',1,'Normalization','probability')
xlabel('Solidification time [min]');
ylabel('Counts');
title('Histogram of all solidification times for 1000 cycles');
axis([0 60 0 1]);
dim = [0.63 0.5 0.4 0.4];


%%
figure(13)
histogram(time_frozen./60.*Dt, 'BinWidth',1)
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 5000]);
title('Histogram of all freezing times for 1000 cycles');

%% 
figure(14)
histogram(temp_nuc, 'BinWidth',0.2)
xlabel('Nucleation Temperature [°C]');
ylabel('Counts');
axis([-20 0 0 5000]);
title('Histogram of all nucleation temperatures for 1000 cycles');


figure(17)
histogram(time_nuc(1,1,1,:)./60.*Dt, 'BinWidth',1)
xlabel('Nucleation time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of corner vial (1,1) nucleation times for 1000 cycles');

%%
figure(18)
histogram(time_solid(1,1,1,:)./60.*Dt,'BinWidth',1)
xlabel('Solidification time [min]');
ylabel('Counts');
axis([0 60 0 500]);


figure(19)
histogram(time_frozen(1,1,1,:)./60.*Dt, 'BinWidth',1)
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of corner vial (1,1) freezing times for 1000 cycles');


%histograms for single corner vials for all runs

figure(20)
histogram(time_nuc(4,4,1,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Nucleation time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of center vial (4,4) nucleation times for 1000 cycles');


figure(21)
histogram(time_solid(4,4,1,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Solidification time [min]');
ylabel('Counts');
axis([0 60 0 500]);
title('Histogram of center vial (4,4) solidification times for 1000 cycles');


figure(22)
histogram(time_frozen(4,4,1,:)./60.*Dt, 'BinWidth',1,'FaceColor', 'magenta')
xlabel('Freezing time [min]');
ylabel('Counts');
axis([60 180 0 500]);
title('Histogram of center vial (4,4) freezing times for 1000 cycles');


%% More figures: Bivariate histograms

Nucedges = linspace(0,200,401);
Soledges = linspace(0,100,201);
Tempedges = linspace(-30,-5,501);

figure(1000)
histogram2(Time_nuc_all(:,:,1,:), Time_solid_all(:,:,1,:),Nucedges,Soledges,'DisplayStyle','tile','Normalization','probability')
title('Nucleation times vs. Solidification times');
xlabel('Nucleation time [min]');
ylabel('Solidification time [min]');
cb = colorbar;
cb.Label.String = 'Normalized probability';


figure(1001)
histogram2(Time_nuc_all(:,:,1,:), Temp_nuc_all(:,:,1,:),Nucedges,Tempedges,'DisplayStyle','tile','Normalization','probability')
title('Nucleation times vs. Nucleation temperatures');
xlabel('Nucleation time [min]');
ylabel('Nucleation temperatures [°C]');
cb = colorbar;
cb.Label.String = 'Normalized probability';

figure(1002)
histogram2(Time_solid_all(:,:,1,:), Temp_nuc_all(:,:,1,:),Soledges,Tempedges,'DisplayStyle','tile','Normalization','probability')
title('Solidification times vs. Nucleation temperatures');
xlabel('Solidification time [min]');
ylabel('Nucleation temperatures [°C]');
cb = colorbar;
cb.Label.String = 'Normalized probability';

hfig = figure(1003);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
sgtitle('Correlations between the characteristic quantities')
subplot(1,3,1)
histogram2(Time_nuc_all(:,:,1,:), Time_solid_all(:,:,1,:),Nucedges,Soledges,'DisplayStyle','tile','Normalization','probability')
title('(a)');
xlabel('Nucleation time [min]');
ylabel('Solidification time [min]');
cb = colorbar;
cb.Label.String = 'Normalized probability';
subplot(1,3,2)
histogram2(Time_nuc_all(:,:,1,:), Temp_nuc_all(:,:,1,:),Nucedges,Tempedges,'DisplayStyle','tile','Normalization','probability')
title('(b)');
xlabel('Nucleation time [min]');
ylabel('Nucleation temperatures [°C]');
cb = colorbar;
cb.Label.String = 'Normalized probability';
subplot(1,3,3)
histogram2(Time_solid_all(:,:,1,:), Temp_nuc_all(:,:,1,:),Soledges,Tempedges,'DisplayStyle','tile','Normalization','probability')
title('(c)');
xlabel('Solidification time [min]');
ylabel('Nucleation temperatures [°C]');
cb = colorbar;
cb.Label.String = 'Normalized probability';

%% Plot temperature curves

figure(100)
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
% plot(t_time./3600,temp_corner_all(:,1),'LineWidth',1)
% plot(t_time./3600,temp_edge_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_center_all(:,1:3),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
%title('(b) Thermally interacting vials: Stochasticity','interpreter','latex','FontSize',14)
title('(a) Thermally independent vials: Stochasticity','interpreter','latex','FontSize',14)
ylabel('Temperature [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
xlabel('Process time [h]','interpreter','latex','FontSize',14)
%legend('Shelf','Corner vial','Edge vial','Center vial')
legend('Shelf','First simulation','Second simulation','Third simulation','interpreter','latex')
axis([0 5 -55 25])
%xticks(linspace(0,16,9))

figure(101)
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_corner_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_edge_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_center_all(:,1),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('(b) Thermally interacting vials: Position','interpreter','latex','FontSize',14)
ylabel('Temperature [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
xlabel('Process time [h]','interpreter','latex','FontSize',14)
legend('Shelf','Corner vial','Edge vial','Center vial','interpreter','latex')
axis([0 15 -55 25])
xticks(linspace(0,16,9))


figure(102)
t = tiledlayout(1,2,'TileSpacing','Compact');
%Tile 1: 0
hn(1) = nexttile();
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_center_all(:,8:10),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('(a)','interpreter','latex','FontSize',14)
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

temp_diff1 = diff(temp_center_all(:,1));
temp_diff2 = diff(diff(temp_center_all(:,1)));

figure(104)
plot(t_time(2:end)./3600,temp_diff1)
axis([2 2.7 -10e-3 10e-3])

figure(105)
plot(t_time(2:end-1)./3600,temp_diff2)
axis([2 2.7 -10e-6 10e-6])
%%
figure(106)
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_center_all(:,1:3),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('Controlled Nucleation: Setup I','interpreter','latex','FontSize',14)
legend('Shelf temperature','First simulation','Second simulation','Third simulation','interpreter','latex','FontSize',12)
xlabel('Process time [hr]','interpreter','latex','FontSize',14)
ylabel('Temperature [$$^{\circ}$$C]','interpreter','latex','FontSize',14)
axis([0 14 -51 21])
grid on;

figure(107)
plot(t_time./3600,T_ext,'LineWidth',1)
hold on;
plot(t_time./3600,temp_corner_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_edge_all(:,1),'LineWidth',1)
plot(t_time./3600,temp_center_all(:,1),'LineWidth',1)
yline(0,'b--','T$$_{eq}$$','interpreter','latex','FontSize',12)
yline(-9,'r--','K = 1h$$^{-1}$$','interpreter','latex','FontSize',12) % Temperature where we have one nucleation event per hour
title('Controlled Nucleation: Setup III','interpreter','latex','FontSize',14)
legend('Shelf temperature','Corner vial','Edge vial','Center vial','interpreter','latex','FontSize',12)
xlabel('Process time [hr]','interpreter','latex','FontSize',14)
ylabel('Temperature [$$^{\circ}$$C]','interpreter','latex','FontSize',14)
axis([0 14 -51 21])
grid on;
%% Save variables

% temp_nuc_20022020_c05 = Temp_nuc_all;
% time_nuc_20022020_c05 = Time_nuc_all;
% time_solid_20022020_c05 = Time_solid_all;
% temp_hold_20022020_c05 = hold_final;
% 
% save('v8_20022020_c05.mat','temp_nuc_20022020_c05',...
% 'time_nuc_20022020_c05','time_solid_20022020_c05','temp_hold_20022020_c05')
% 
% 
% temp_nuc_20022000_c05_h1m5 = Temp_nuc_all;
% time_nuc_20022000_c05_h1m5 = Time_nuc_all;
% time_solid_20022000_c05_h1m5 = Time_solid_all;
% temp_hold_20022000_c05_h1m5 = hold_final;
% 
% save('v8_20022000_c05_h1m5.mat','temp_nuc_20022000_c05_h1m5',...
% 'time_nuc_20022000_c05_h1m5','time_solid_20022000_c05_h1m5','temp_hold_20022000_c05_h1m5')
% 
% 
% temp_nuc_20022020_c05_h540m12 = Temp_nuc_all;
% time_nuc_20022020_c05_h540m12 = Time_nuc_all;
% time_solid_20022020_c05_h540m12 = Time_solid_all;
% temp_hold_20022020_c05_h540m12 = hold_final;
% 
% save('v8e_20022020_c05_h540m12.mat','temp_nuc_20022020_c05_h540m12',...
% 'time_nuc_20022020_c05_h540m12','time_solid_20022020_c05_h540m12','temp_hold_20022020_c05_h540m12')

%%
% temp_nuc_20022020_c05_c180m5h360m10 = Temp_nuc_all;
% time_nuc_20022020_c05_c180m5h360m10 = Time_nuc_all;
% time_solid_20022020_c05_c180m5h360m10 = Time_solid_all;
% temp_hold_20022020_c05_c180m5h360m10 = hold_final;
% 
% save('ssie_20022020_c05_c180m5h360m10.mat','temp_nuc_20022020_c05_c180m5h360m10',...
% 'time_nuc_20022020_c05_c180m5h360m10','time_solid_20022020_c05_c180m5h360m10','temp_hold_20022020_c05_c180m5h360m10')

% temp_nuc_20020020_c05_c180m5ramp = Temp_nuc_all;
% time_nuc_20020020_c05_c180m5ramp = Time_nuc_all;
% time_solid_20020020_c05_c180m5ramp = Time_solid_all;
% temp_hold_20020020_c05_c180m5ramp = hold_final;
% 
% save('v8_20020020_c05_c180m5ramp.mat','temp_nuc_20020020_c05_c180m5ramp',...
% 'time_nuc_20020020_c05_c180m5ramp','time_solid_20020020_c05_c180m5ramp','temp_hold_20020020_c05_c180m5ramp')
% 

% temp_nuc_20022020_c05_h640m12 = Temp_nuc_all;
% time_nuc_20022020_c05_h640m12 = Time_nuc_all;
% time_solid_20022020_c05_h640m12 = Time_solid_all;
% temp_hold_20022020_c05_h640m12 = hold_final;
% 
% save('ssie_20022020_c05_h640m12.mat','temp_nuc_20022020_c05_h640m12',...
% 'time_nuc_20022020_c05_h640m12','time_solid_20022020_c05_h640m12','temp_hold_20022020_c05_h640m12')
% 
% temp_nuc_20003000_c05_low = Temp_nuc_all;
% time_nuc_20003000_c05_low = Time_nuc_all;
% time_solid_20003000_c05_low = Time_solid_all;
% temp_hold_20003000_c05_low = hold_final;
% 
% save('ssie_20003000_c05_low.mat','temp_nuc_20003000_c05_low',...
% 'time_nuc_20003000_c05_low','time_solid_20003000_c05_low','temp_hold_20003000_c05_low')
% 

% temp_nuc_20002000_c05_4900s = Temp_nuc_all;
% time_nuc_20002000_c05_4900s = Time_nuc_all;
% time_solid_20002000_c05_4900s = Time_solid_all;
% temp_hold_20002000_c05_4900s = hold_final;
% 
% save('ssie_20002000_c05_4900s.mat','temp_nuc_20002000_c05_4900s',...
% 'time_nuc_20002000_c05_4900s','time_solid_20002000_c05_4900s','temp_hold_20002000_c05_4900s')
% 
% temp_nuc_20000000_c05n_dt100 = Temp_nuc_all;
% time_nuc_20000000_c05n_dt100 = Time_nuc_all;
% time_solid_20000000_c05n_dt100 = Time_solid_all;
% temp_hold_20000000_c05n_dt100 = hold_final;
% 
% save('ssie_20000000_c05n_dt100.mat','temp_nuc_20000000_c05n_dt100',...
% 'time_nuc_20000000_c05n_dt100','time_solid_20000000_c05n_dt100','temp_hold_20000000_c05n_dt100')

% temp_profile05 = temp_center_all(:,1);
% save('tempcurvenonuc05.mat','temp_profile05')
% % % 
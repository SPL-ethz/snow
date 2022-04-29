%% Numerical validation study

clear all; close all; clc;

%% Loading

%Matlab files


load('c3D_122003_xzh_101010_m40_dt5_20000_N1000_20220422.mat')

%%Python output, csv files

boxm40_tempnuc = csvread('box40_tempnuc.csv');
boxm40_timenuc = csvread('box40_timenuc.csv');
boxm40_tsol = csvread('box40_tsol.csv');

%% Figures

Tempedges = linspace(-30,0,91);
Nucedges = linspace(0,1000,1001);
Soledges = linspace(0,200,201);

figure(1)
% t = tiledlayout(2,2,'TileSpacing','Compact');
% hn(1) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20000000_c05(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Tempedges(1:end-1),histcounts(Tnuc_HT1,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% xlabel(t,'$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
% ylabel(t,'$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
% title('(a) Constant shelf HT','interpreter','latex','FontSize',14)
% grid on
% axis([-15 -5 0 0.6])
% xticks([-15 -13 -11 -9 -7 -5])
% hn(2) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20000020_c05(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Tempedges(1:end-1),histcounts(Tnuc_HText,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% grid on
% title('(b) k$$_{ext} = 20~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% Lgnd = legend('MATLAB','python','interpreter','latex','FontSize',11,'Box','off');
% axis([-15 -5 0 0.6])
% xticks([-15 -13 -11 -9 -7 -5])
% hn(3) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20002000_c05(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% grid on
% hold on
% plot(Tempedges(1:end-1),histcounts(Tnuc_HTint,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% title('(c) k$$_{int} = 20~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% axis([-15 -5 0 0.6])
% xticks([-15 -13 -11 -9 -7 -5])
% hn(4) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20020000_c05(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Tempedges(1:end-1),histcounts(Tnuc_HTvar,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% grid on
% title('(d) s$$_{sh} = 2~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% axis([-15 -5 0 0.6])
% xticks([-15 -13 -11 -9 -7 -5])
% 
% figure(2)
% t = tiledlayout(2,2,'TileSpacing','Compact');
% hn(1) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20000000_c05(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Nucedges(1:end-1),histcounts(tnuc_HT1./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% xlabel(t,'$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel(t,'$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% axis([80 110 0 0.2])
% title('(a) Constant shelf HT','interpreter','latex','FontSize',14)
% grid on
% hn(2) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20000020_c05(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Nucedges(1:end-1),histcounts(tnuc_HText./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% xlabel(t,'$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel(t,'$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% axis([60 110 0 0.12])
% title('(b) k$$_{ext} = 20~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% grid on
% hn(3) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20002000_c05(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Nucedges(1:end-1),histcounts(tnuc_HTint./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% xlabel(t,'$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel(t,'$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% axis([80 150 0 0.12])
% title('(c) k$$_{int} = 20~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% grid on
% hn(4) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20020000_c05(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on
% plot(Nucedges(1:end-1),histcounts(tnuc_HTvar./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% xlabel(t,'$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel(t,'$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% axis([80 130 0 0.12])
% title('(d) s$$_{sh} = 2~Wm^{-2}K^{-1}$$','interpreter','latex','FontSize',14)
% grid on
% Lgnd = legend('MATLAB','python','interpreter','latex','FontSize',11,'Box','off');
% 
% 
% figure(3)
% plot(Soledges(1:end-1),histcounts(time_solid_20020000_c05(:),Soledges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Soledges(1:end-1),histcounts(tsol_HTvar./60,Soledges,'Normalization','pdf'),'r','LineWidth',1)
% 
% hold off;
% xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 65 0 0.15])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(e)','interpreter','latex','FontSize',12)
% xticks([0 20 40 60])
% grid on
% 
% Nucedges = linspace(0,1000,251);
% 
% figure(4)
% t = tiledlayout(1,3,'TileSpacing','Compact');
% hn(1) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20002000_c05_h540m10(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Tempedges(1:end-1),histcounts(Tnuc_HTadnuc,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% hold off;
% xlabel('$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([-15 0 0 0.6])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(a)','interpreter','latex','FontSize',12)
% grid on
% hn(2) = nexttile();
% plot(Soledges(1:end-1),histcounts(time_solid_20002000_c05_h540m10(:),Soledges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Soledges(1:end-1),histcounts(tsol_HTadnuc./60,Soledges,'Normalization','pdf'),'r','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 150 0 0.03])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(b)','interpreter','latex','FontSize',12)
% grid on
% hn(3) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20002000_c05_h540m10(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Nucedges(1:end-1),histcounts(tnuc_HTadnuc./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 700 0 0.012])
% %title(Lgnd,'$$k_\mathrm{int} \left[ \frac{\mathrm{W}}{\mathrm{m^2K}} \right]$$','interpreter','latex')
% title('(c)','interpreter','latex','FontSize',12)
% grid on
% 
% Tempedges = linspace(-30,0,301);
% Nucedges = linspace(0,1000,1001);

t = tiledlayout(1,3,'TileSpacing','Compact');
hn(1) = nexttile();
plot(Tempedges(1:end-1),histcounts(temp_nuc_122003_101010_m40(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
hold on;
plot(Tempedges(1:end-1),histcounts(boxm40_tempnuc,Tempedges,'Normalization','pdf'),'r','LineWidth',1)
hold off;
xlabel('$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
ylabel('$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
%Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
%Lgnd = legend('0','10','20','30');
axis([-15 0 0 0.6])
%title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
title('(a)','interpreter','latex','FontSize',12)
grid on
hn(2) = nexttile();
plot(Soledges(1:end-1),histcounts(time_solid_122003_101010_m40(:),Soledges,'Normalization','pdf'),'k','LineWidth',1)
hold on;
plot(Soledges(1:end-1),histcounts(boxm40_tsol./60,Soledges,'Normalization','pdf'),'r','LineWidth',1)
hold off;
xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
axis([0 120 0 0.06])
%title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
title('(b)','interpreter','latex','FontSize',12)
set(gcf,'units','centimeters','position',[5,5,16,6])
grid on
hn(3) = nexttile();
plot(Nucedges(1:end-1),histcounts(time_nuc_122003_101010_m40(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
hold on;
plot(Nucedges(1:end-1),histcounts(boxm40_timenuc./60,Nucedges,'Normalization','pdf'),'r','LineWidth',1)
hold off;
xlabel('$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
ylabel('$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
%Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
%Lgnd = legend('0','10','20','30');
axis([0 350 0 0.03])
%title(Lgnd,'$$k_\mathrm{int} \left[ \frac{\mathrm{W}}{\mathrm{m^2K}} \right]$$','interpreter','latex')
title('(c)','interpreter','latex','FontSize',12)
grid on


%%
% figure(12040)
% t = tiledlayout(1,3,'TileSpacing','Compact');
% hn(1) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20000000_c05(:),Tempedges,'Normalization','pdf'),'c','LineWidth',1)
% hold on;
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20001000_c05(:),Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20002000_c05(:),Tempedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20003000_c05(:),Tempedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([-15 0 0 0.6])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(a)','interpreter','latex','FontSize',12)
% grid on
% hn(2) = nexttile();
% plot(Soledges(1:end-1),histcounts(time_solid_20001000_c05(:),Soledges,'Normalization','pdf'),'r','LineWidth',1)
% hold on;
% plot(Soledges(1:end-1),histcounts(time_solid_20000000_c05(:),Soledges,'Normalization','pdf'),'c','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20002000_c05(:),Soledges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20003000_c05(:),Soledges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 65 0 0.3])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(b)','interpreter','latex','FontSize',12)
% grid on
% xticks([0 20 40 60])
% 
% hn(3) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20000000_c05(:),Nucedges,'Normalization','pdf'),'c','LineWidth',1)
% hold on;
% plot(Nucedges(1:end-1),histcounts(time_nuc_20001000_c05(:),Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20002000_c05(:),Nucedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20003000_c05(:),Nucedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([80 170 0 0.2])
% %title(Lgnd,'$$k_\mathrm{int} \left[ \frac{\mathrm{W}}{\mathrm{m^2K}} \right]$$','interpreter','latex')
% title('(c)','interpreter','latex','FontSize',12)
% grid on
% xticks([80 110 140 170])
% 
% 
% [leg1,icons,plots,txt] = legend({'$$k_\mathrm{int} = 0 \frac{\mathrm{W}}{\mathrm{m^2K}}$$','$$k_\mathrm{int} = 10$$','$$k_\mathrm{int} = 20$$','$$k_\mathrm{int} = 30$$'},'interpreter','latex','FontSize',11);
% p1 = icons(1).Position;
% p2 = icons(2).Position;
% p3 = icons(3).Position;
% p4 = icons(4).Position;
% icons(1).Position = [0.45 p1(2) 0];
% icons(2).Position = [0.45 p2(2) 0];
% icons(3).Position = [0.45 p3(2) 0];
% icons(4).Position = [0.45 p4(2) 0];
% icons(5).XData = [0.35 0.4];
% icons(7).XData = [0.35 0.4];
% icons(9).XData = [0.35 0.4];
% icons(11).XData = [0.35 0.4];
% set(leg1,'Box','off')
% 
% 
% 
% figure(1204)
% t = tiledlayout(2,3,'TileSpacing','Compact');
% hn(1) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20000000_c05(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20001000_c05(:),Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20002000_c05(:),Tempedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20003000_c05(:),Tempedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([-30 -5 0 0.6])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(a)','interpreter','latex','FontSize',12)
% grid on
% hn(2) = nexttile();
% plot(Soledges(1:end-1),histcounts(time_solid_20000000_c05(:),Soledges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Soledges(1:end-1),histcounts(time_solid_20001000_c05(:),Soledges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20002000_c05(:),Soledges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20003000_c05(:),Soledges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 65 0 0.3])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(b)','interpreter','latex','FontSize',12)
% grid on
% xticks([0 20 40 60])
% 
% hn(3) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20000000_c05(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Nucedges(1:end-1),histcounts(time_nuc_20001000_c05(:),Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20002000_c05(:),Nucedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20003000_c05(:),Nucedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([80 190 0 0.2])
% %title(Lgnd,'$$k_\mathrm{int} \left[ \frac{\mathrm{W}}{\mathrm{m^2K}} \right]$$','interpreter','latex')
% title('(c)','interpreter','latex','FontSize',12)
% grid on
% xticks([80 110 140 170])
% 
% [leg1,icons,plots,txt] = legend({'$$k_\mathrm{int} = 0 \frac{\mathrm{W}}{\mathrm{m^2K}}$$','$$k_\mathrm{int} = 10$$','$$k_\mathrm{int} = 20$$','$$k_\mathrm{int} = 30$$'},'interpreter','latex','FontSize',11);
% p1 = icons(1).Position;
% p2 = icons(2).Position;
% p3 = icons(3).Position;
% p4 = icons(4).Position;
% icons(1).Position = [0.45 p1(2) 0];
% icons(2).Position = [0.45 p2(2) 0];
% icons(3).Position = [0.45 p3(2) 0];
% icons(4).Position = [0.45 p4(2) 0];
% icons(5).XData = [0.35 0.4];
% icons(7).XData = [0.35 0.4];
% icons(9).XData = [0.35 0.4];
% icons(11).XData = [0.35 0.4];
% set(leg1,'Box','off')
% 
% 
% 
% hn(4) = nexttile();
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20000000_c05_low(:),Tempedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20001000_c05_low(:),Tempedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20002000_c05_low(:),Tempedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Tempedges(1:end-1),histcounts(temp_nuc_20003000_c05_low(:),Tempedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$T^\mathrm{nuc}$$ [$$^{\circ}$$C]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nT}}$$ [K$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([-30 -5 0 0.3])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(d)','interpreter','latex','FontSize',12)
% grid on
% 
% hn(5) = nexttile();
% plot(Soledges(1:end-1),histcounts(time_solid_20000000_c05_low(:),Soledges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Soledges(1:end-1),histcounts(time_solid_20001000_c05_low(:),Soledges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20002000_c05_low(:),Soledges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Soledges(1:end-1),histcounts(time_solid_20003000_c05_low(:),Soledges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{sol}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{sol}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([0 65 0 0.15])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(e)','interpreter','latex','FontSize',12)
% xticks([0 20 40 60])
% grid on
% 
% hn(6) = nexttile();
% plot(Nucedges(1:end-1),histcounts(time_nuc_20000000_c05_low(:),Nucedges,'Normalization','pdf'),'k','LineWidth',1)
% hold on;
% plot(Nucedges(1:end-1),histcounts(time_nuc_20001000_c05_low(:),Nucedges,'Normalization','pdf'),'r','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20002000_c05_low(:),Nucedges,'Normalization','pdf'),'m','LineWidth',1)
% plot(Nucedges(1:end-1),histcounts(time_nuc_20003000_c05_low(:),Nucedges,'Normalization','pdf'),'b','LineWidth',1)
% hold off;
% xlabel('$$t^\mathrm{nuc}$$ [min]','FontSize',14,'interpreter','latex');
% ylabel('$$f_{\mathrm{nt}}$$ [min$$^{-1}$$]','interpreter','latex','FontSize',14);
% %Lgnd = legend('k$$_{int}$$ = 0 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 10 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 20 $$\frac{W}{m^2K}$$','k$$_{int}$$ = 30 $$\frac{W}{m^2K}$$','interpreter','latex');
% %Lgnd = legend('0','10','20','30');
% axis([80 190 0 0.1])
% %title(Lgnd,'k$$_{int} \left[ \frac{W}{m^2K} \right]$$','interpreter','latex')
% title('(f)','interpreter','latex','FontSize',12)
% grid on
% xticks([80 110 140 170])
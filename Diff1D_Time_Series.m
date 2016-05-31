
%This script makes simple time series plots of shear etc. from the
%Diff1D model.

depth = -75;

%Setup figure:
figure;
set(gcf,'Position',[457           1        1102         973]);
set(gcf,'Position',get(0,'ScreenSize'));

h1 = subplot('Position',[0.1300    0.7139    0.7750    0.1773]);
xlim([135 165]);
ylim([-0.04 0]);
set(h1,'xtick',[135:5:165]);
set(h1,'xticklabel',[]);
set(h1,'ytick',[-0.04:0.01:-0.01]);
ylabel('$\partial u/\partial z\,\,/\,\,$s$^{-1}$');
hold on;
box on;
h2 = subplot('Position',[0.1300    0.5126   0.7750    0.1773]);
xlim([135 165]);
ylim([-0.7e-3 1.1e-3]);
set(h2,'xtick',[135:5:165]);
set(h2,'xticklabel',[]);
set(h2,'ytick',[-5:5:10]*1e-4);
set(h2,'yticklabel',[-5:5:10]);
ylabel('$Sh^2_{red}\,\,/\,\,10^{-4}$s$^{-2}$');
hold on;
box on;
h3 = subplot('Position',[0.1300    0.3113    0.7750    0.1773]);
xlim([135 165]);
ylim([0 2.1e-3]);
set(h3,'xtick',[135:5:165]);
set(h3,'xticklabel',[]);
set(h3,'ytick',[0:0.5:2]*1e-3);
set(h3,'yticklabel',[0:0.5:2]);
ylabel('$\kappa_T\,\,/\,\,10^{-3}$m$^2$s$^{-1}$');
hold on;
box on;
h4 = subplot('Position',[0.1300    0.1100    0.7750    0.1773]);
xlim([135 165]);
ylim([-400 0]);
set(h4,'xtick',[135:5:165]);
% $$$ set(h1,'ytick',[-0.04:0.01:-0.01]);
ylabel('$J_q\,\,/\,\,$Wm$^{-2}$');
xlabel('Day');
hold on;
box on;

load('KPP_nDIR_strun.mat');
[tmp dind] = min(abs(z_w-depth))
depth = z_w(dind)

%derive variables:
dudz = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
dTdz = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
RSh2 = dudz.^2-4*N2;
Ri = N2./(dudz.^2);

plot(t/86400,dudz(dind,:),'-k','LineWidth',2,'Parent',h1);
plot(t/86400,RSh2(dind,:),'-k','LineWidth',2,'Parent',h2);
plot(t/86400,kt(dind,:),'-k','LineWidth',2,'Parent',h3);
plot(t/86400,-rho0*Cp*kt(dind,:).*dTdz(dind,:),'-k','LineWidth',2,'Parent',h4);

load('KPP_nDIR_SYM2_strun.mat');
[tmp dind] = min(abs(z_w-depth))
depth = z_w(dind)

%derive variables:
dudz = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
dTdz = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
RSh2 = dudz.^2-4*N2;
Ri = N2./(dudz.^2);

plot(t/86400,dudz(dind,:),'--k','LineWidth',2,'Parent',h1);
plot(t/86400,RSh2(dind,:),'--k','LineWidth',2,'Parent',h2);
plot(t/86400,kt(dind,:),'--k','LineWidth',2,'Parent',h3);
plot(t/86400,-rho0*Cp*kt(dind,:).*dTdz(dind,:),'--k','LineWidth',2,'Parent',h4);

[0.4941    0.1843    0.5569]
load('PP_nDIR_strun.mat');
[tmp dind] = min(abs(z_w-depth))
depth = z_w(dind)

%derive variables:
dudz = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
dTdz = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
RSh2 = dudz.^2-4*N2;
Ri = N2./(dudz.^2);

plot(t/86400,dudz(dind,:),'-','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h1);
plot(t/86400,RSh2(dind,:),'-','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h2);
plot(t/86400,kt(dind,:),'-','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h3);
plot(t/86400,-rho0*Cp*kt(dind,:).*dTdz(dind,:),'-','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h4);

load('PP_nDIR_SYM2_strun.mat');
[tmp dind] = min(abs(z_w-depth))
depth = z_w(dind)

%derive variables:
dudz = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
dTdz = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end- ...
                                                  1)),[1 length(u(1,:))]);
RSh2 = dudz.^2-4*N2;
Ri = N2./(dudz.^2);

plot(t/86400,dudz(dind,:),'--','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h1);
plot(t/86400,RSh2(dind,:),'--','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h2);
plot(t/86400,kt(dind,:),'--','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h3);
plot(t/86400,-rho0*Cp*kt(dind,:).*dTdz(dind,:),'--','color',[0.4941 0.1843 0.5569],'LineWidth',2,'Parent',h4);

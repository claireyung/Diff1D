%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    DIFF1Dplot                                    %
%                                                                  %
% This script sets up the plotting                                 %
% time step from DIFF1D                                            %
%                                                                  %
% Ryan Holmes June 2014                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup:
figure;
set(gcf,'Position',[76 35 1840 931]);

PosVec = [0.0718    0.5539    0.1898    0.4044; ...
          0.2904    0.5539    0.1898    0.4044; ...
          0.5025    0.5539    0.1898    0.4044; ...
          0.7201    0.5539    0.1898    0.4044; ...
          0.0718    0.0789    0.1898    0.4044; ...
          0.2904    0.0789    0.1898    0.4044; ...
          0.5025    0.0789    0.1898    0.4044; ...
          0.7201    0.0789    0.1898    0.4044];
xlims = [-1.5 1.5;...
         13 25; ...
         34.9 35.2; ...
         -0.04 0.02; ...
         0 3e-4; ...
         -10e-4 5e-4; ...
         0 6e-3; ...
         0 300; ...
         ]
% $$$          -5e-6 5e-6; ...
ylims = [-300 0];
txty = -285;

%u and v:
Pa1 = subplot(2,4,1);
plot(uI,zI,'-k','LineWidth',2);
hold on;
plot(vI,zI,'-r','LineWidth',2);
Pa11 = plot(uI,zI,'-k','LineWidth',2);
Pa12 = plot(vI,zI,'-r','LineWidth',2);
ylim(ylims);
xlim(xlims(1,:));
set(gca,'ytick',[ylims(1):20:ylims(2)]);
ylabel('Depth (m)');
set(gca,'Position',PosVec(1,:));
grid on;
text(-0.9,txty,'Velocity (ms$^{-1}$)');
text(0.5,-40,'$u$','FontSize',25);
text(-0.5,-60,txty,'$v$','FontSize',25,'Color','r');
hold on;
Pa2 = subplot(2,4,2);
plot(TI,zI,'-k','LineWidth',2);
hold on;
Pa21 = plot(TI,zI,'-k','LineWidth',2);
ylim(ylims);
xlim(xlims(2,:));
text(16,txty,'Temperature ($^\circ$C)');
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
grid on;
set(gca,'Position',PosVec(2,:));
hold on;
Pa3 = subplot(2,4,3);
plot(SI,zI,'-k','LineWidth',2);
hold on;
Pa31 = plot(SI,zI,'-k','LineWidth',2);
ylim(ylims);
xlim(xlims(3,:));
text(35.05,txty,'Salinity (psu)');
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
grid on;
set(gca,'Position',PosVec(3,:));
hold on;
Pa4 = subplot(2,4,4);
plot((uI(2:end)-uI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'k'],'LineWidth',2);
hold on;
plot((vI(2:end)-vI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'r'],'LineWidth',2);
Pa41 = plot((uI(2:end)-uI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'k'],'LineWidth',2);
Pa42 = plot((vI(2:end)-vI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'r'],'LineWidth',2);
xlim(xlims(4,:));
ylim(ylims);
text(-0.03,txty,'Vertical Shear ($s^{-1}$)');
text(-0.02,-140,'$du/dz$','FontSize',25);
text(0.005,-40,'$dv/dz$','FontSize',25,'Color','r');
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
grid on;
set(gca,'Position',PosVec(4,:));
hold on;
Pa5 = subplot(2,4,5);
plot((bI(2:end)-bI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'k'],'LineWidth',2);
hold on;
Pa51 = plot((bI(2:end)-bI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))),zwI,['-' ...
                    'k'],'LineWidth',2);
xlim(xlims(5,:));
ylim(ylims);
text(2e-4,txty,'$N^2$ ($s^{-2}$)');
set(gca,'ytick',[ylims(1):20:ylims(2)]);
ylabel('Depth (m)');
grid on;
set(gca,'Position',PosVec(5,:));
hold on;
Pa6 = subplot(2,4,6);
plot(((uI(2:end)-uI(1:(end-1)))./(zI(2:end)-zI(1:(end-1)))+...
     (vI(2:end)-vI(1:(end-1)))./(zI(2:end)-zI(1:(end-1)))).^2-4*...
     (bI(2:end)-bI(1:(end-1)))./(zI(2:end)-zI(1:(end-1)))...
         ,zwI,'-k','LineWidth',2);
hold on;
Pa61 = plot(((uI(2:end)-uI(1:(end-1)))./(zI(2:end)-zI(1:(end-1)))+ ...
            (vI(2:end)-vI(1:(end-1)))./(zI(2:end)-zI(1:(end-1)))).^2-4* ...
           (bI(2:end)-bI(1:(end-1)))./(zI(2:end)-zI(1:(end-1))) ...
           ,zwI,'-k','LineWidth',2);
xlim(xlims(6,:));
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
grid on;
set(gca,'Position',PosVec(6,:));
hold on;
ylim(ylims);
text(-1e-4,txty,'$Sh^2_{red}$ ($s^{-2}$)');
drawnow;

Pa7 = subplot(2,4,7);
plot(kv(:,1),z_w,'-k','LineWidth',2);
hold on;
plot(kt(:,1),z_w,'-b','LineWidth',2);
Pa71=plot(kv(:,1),z_w,'-k','LineWidth',2);
Pa72=plot(kt(:,1),z_w,'-b','LineWidth',2);
ylim(ylims);
xlim(xlims(7,:));
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
set(gca,'Position',PosVec(7,:));
grid on;
text(1e-3,txty,'$\kappa$ (m$^2$s$^{-1}$)');
text(1e-3,-110,'$\kappa_v$','FontSize',25);
text(1e-3,-140,'$\kappa_t$','FontSize',25,'color','b');
hold on;

%Momentum equation terms:
% $$$ Pa8 = subplot(2,4,8);
% $$$ plot(UTND(:,1),z_rho,'-k','LineWidth',2);
% $$$ hold on;
% $$$ plot(URST(:,1),z_rho,'-b','LineWidth',2);
% $$$ plot(UMFX(:,1),z_rho,'-r','LineWidth',2);
% $$$ plot(UDIV(:,1),z_rho,'-','color',[0 0.5 0],'LineWidth',2);
% $$$ plot(UPGF(:,1),z_rho,'-m','LineWidth',2);
% $$$ Pa81=plot(UTND(:,1),z_rho,'-k','LineWidth',2);
% $$$ Pa82=plot(URST(:,1),z_rho,'-b','LineWidth',2);
% $$$ Pa83=plot(UMFX(:,1),z_rho,'-r','LineWidth',2);
% $$$ Pa84=plot(UDIV(:,1),z_rho,'-','color',[0 0.5 0],'LineWidth',2);
% $$$ Pa85=plot(UPGF(:,1),z_rho,'-m','LineWidth',2);
% $$$ ylim(ylims);
% $$$ xlim(xlims(8,:));
% $$$ set(gca,'ytick',[ylims(1):20:ylims(2)]);
% $$$ set(gca,'yticklabel',[]);
% $$$ set(gca,'Position',PosVec(8,:));
% $$$ grid on;
% $$$ text(-4.5e-6,txty,'Force bal. (ms$^{-2}$)');
% $$$ text(-3e-6,-110,'UTND','FontSize',25);
% $$$ text(-3e-6,-140,txty,'URST','FontSize',25,'Color','b');
% $$$ text(-3e-6,-170,txty,'UMFX','FontSize',25,'Color','r');
% $$$ text(-3e-6,-200,txty,'UDIV','FontSize',25,'Color',[0 0.5 0]);
% $$$ text(-3e-6,-230,txty,'UPGF','FontSize',25,'Color','m');
% $$$ hold on;

%heat fluxes:
Pa8 = subplot(2,4,8);
plot(kt(2:(end-1),1).*Cp.*rho0.*(T(2:end,1)-T(1:(end-1),1))./ ...
     (z_rho(2:end)-z_rho(1:(end-1))),z_w(2:(end-1)),'-k','LineWidth',2);
hold on;
Pa81 = plot(kt(2:(end-1),1).*Cp.*rho0.*(T(2:end,1)-T(1:(end-1),1))./(z_rho(2:end)-z_rho(1:(end-1))),z_w(2:(end-1)),'-k','LineWidth',2);
ylim(ylims);
xlim(xlims(8,:));
set(gca,'ytick',[ylims(1):20:ylims(2)]);
set(gca,'yticklabel',[]);
set(gca,'Position',PosVec(8,:));
grid on;
text(150,txty,'Jq (Wm$^{-2}$)');
hold on;



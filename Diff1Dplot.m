%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    DIFF1Dplot                                    %
%                                                                  %
% This script updates the plot for the current time step           %
%                                                                  %
% Ryan Holmes June 2014                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%u and v:
delete(Pa11);
delete(Pa12);
Pa11 = line(mean(u(:,ti),2),z_rho,'LineWidth',2,'color','k','LineStyle','--','Parent',Pa1);
Pa12 = line(mean(v(:,ti),2),z_rho,'LineWidth',2,'color','r','LineStyle','--','Parent',Pa1);

%T:
delete(Pa21);
Pa21 = line(mean(T(:,ti),2),z_rho,'color','k','LineWidth',2,'LineStyle','--','Parent',Pa2);

%S:
delete(Pa31);
Pa31 = line(mean(S(:,ti),2),z_rho,'color','k','LineWidth',2,'LineStyle','--','Parent',Pa3);

%Shears:
delete(Pa41);
delete(Pa42);
Pa41 = line(mean((u(2:end,ti)-u(1:(end-1),ti))./repmat((z_rho(2:end)-z_rho(1:(end-1))),[1 length(ti)]),2),z_w(2:(end-1)),...
           'color','k','LineWidth',2,'LineStyle','--','Parent',Pa4);
Pa42 = line(mean((v(2:end,ti)-v(1:(end-1),ti))./repmat((z_rho(2:end)-z_rho(1:(end-1))),[1 length(ti)]),2),z_w(2:(end-1)),...
           'color','r','LineWidth',2,'LineStyle','--','Parent',Pa4);

%Stratification:
delete(Pa51);
Pa51 = line(mean((b(2:end,ti)-b(1:(end-1),ti))./repmat((z_rho(2:end)-z_rho(1:(end-1))),[1 length(ti)]),2),z_w(2:(end-1)),...
           'color','k','LineWidth',2,'LineStyle','--','Parent', ...
            Pa5);

%Reduced shear squared:
delete(Pa61);
Pa61 = line(mean(((u(2:end,ti)-u(1:(end-1),ti))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(ti)])+...
     (v(2:end,ti)-v(1:(end-1),ti))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(ti)])).^2-4*...
     (b(2:end,ti)-b(1:(end-1),ti))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(ti)]),2)...
         ,z_w(2:(end-1)),'color','k','LineWidth',2,'LineStyle','--','Parent',Pa6);


%Diffusivities:
delete(Pa71);delete(Pa72);
Pa71=line(mean(kv(:,ti),2),z_w,'color','k','LineWidth',2,'LineStyle','--','Parent',Pa7);
Pa72=line(mean(kt(:,ti),2),z_w,'color','b','LineWidth',2,'LineStyle','--','Parent',Pa7);

% $$$ %Force balance:
% $$$ delete(Pa81);delete(Pa82);delete(Pa83);delete(Pa84);
% $$$ Pa81=line(mean(UTND(:,ti),2),z_rho,'color','k','LineWidth',2,'LineStyle','--','Parent',Pa8);
% $$$ Pa82=line(mean(URST(:,ti),2),z_rho,'color','b','LineWidth',2,'LineStyle','--','Parent',Pa8);
% $$$ Pa83=line(mean(UMFX(:,ti),2),z_rho,'color','r','LineWidth',2,'LineStyle','--','Parent',Pa8);
% $$$ Pa84=line(mean(UDIV(:,ti),2),z_rho,'color',[0 0.5 0],'LineWidth',2,'LineStyle','--','Parent',Pa8);
% $$$ Pa85=line(mean(UPGF(:,ti),2),z_rho,'color','m','LineWidth',2,'LineStyle','--','Parent',Pa8);

%Jq:
delete(Pa81);
Pa81 = plot(mean(kt(2:(end-1),ti).*Cp.*rho0.*(T(2:end,ti)-T(1:(end- ...
                                                  1),ti))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(ti)]),2),z_w(2:(end-1)),'--k','LineWidth',2);

drawnow;

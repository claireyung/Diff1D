figure;
set(gcf,'Position',[453 26 1006 650]);
% $$$ set(gcf,'Position',[206 69 1619 905]);
[Tr,Zr] = meshgrid(t/86400,z_rho);
[Tw,Zw] = meshgrid(t/86400,z_w);
axs = [15 t(end)/86400 -250 0];
intp = 0;
lnfilt = 0;
txtx = (t(end)/86400)*0.995;
txty = -200;

%Times to plot:
% $$$ [tmp tII] = min(abs(t/86400-135));
% $$$ [tmp tFF] = min(abs(t/86400-165));
[tmp tII] = min(abs(t/86400));
[tmp tFF] = min(abs(t/86400-120));
tplot = 15;
tvec = tII:tplot:tFF;

%Derive variables:
DUDZ = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
DTDZ = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
DVDZ = (v(2:end,:)-v(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
RSh2 = DUDZ.^2+DVDZ.^2-4*N2;
Ri = 1./(N2./(DUDZ.^2+DVDZ.^2));
Jq = -kt(2:(end-1),:).*DTDZ*Cp*rho0 -kt(2:(end-1),:).*gamt(2:(end-1),:)*Cp*rho0;
Jq_bl = -kt_bl(2:(end-1),:).*DTDZ*Cp*rho0;
Jq_shear = -kt_int(2:(end-1),:).*DTDZ*Cp*rho0;
Jq_nonlocal = -kt(2:(end-1),:).*gamt(2:(end-1),:)*Cp*rho0;
Forcing = srflux +shflux;

%Get EUC:
[zL,tL] = size(u);
EUC = zeros(tL,1);
for x=1:tL
    tmp = u(:,x);%filter_field(U(x,:),3,'-t');
    [tmp2 ind] = max(tmp);
    EUC(x) = z_rho(ind);
end
EUC = filter_field(EUC,5,'-t');

% 7-panels:
% $$$ naxs = 7;
% $$$ % $$$           0.1191    0.8700    0.6861    0.0973;...
% $$$ PosVec = [0.1191    0.87      0.7509    0.12; ...
% $$$           0.1191    0.735      0.7509     0.12; ...
% $$$           0.1191    0.6      0.7509    0.12; ...
% $$$           0.1191    0.465      0.7509    0.12; ...
% $$$           0.1191    0.33      0.7509    0.12; ...
% $$$           0.1191    0.195      0.7509    0.12; ...
% $$$           0.1191    0.06      0.7509    0.12];
% $$$ xlab = [0 0 0 0 0 0 1];
% $$$ xtic = [0 0 0 0 0 0 1];
% $$$ VarOp{1} = {'dvdy','Tr','Zr'};
% $$$ VarOp{2} = {'T','Tr','Zr'};
% $$$ VarOp{3} = {'DUDZ','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ VarOp{4} = {'N2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ VarOp{5} = {'RSh2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ VarOp{6} = {'kt','Tw','Zw'};
% $$$ VarOp{7} = {'Jq','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ caxs = [-3e-6 3e-6;...
% $$$         20 25; ...
% $$$         -0.05 0.03; ...
% $$$         0 3e-4; ...
% $$$         -2e-3 1e-3;...
% $$$         0 2.2e-3;...
% $$$         0 400;];
% $$$ names = {'$\partial$v/$\partial$y$\,\,/\,\,$s$^{-1}$',...
% $$$          'T$\,\,/\,\,^\circ$C', ...
% $$$          '$\partial$u/$\partial$z$\,\,/\,\,$s$^{-1}$', ...
% $$$          'N$^2\,\,/\,\,$s$^{-2}$', ...
% $$$          'Sh$^2_{red}\,\,/\,\,10^{-4}$s$^{-2}$',...
% $$$         '$\kappa_v\,\,/\,\,$m$^2$s$^{-1}$',...
% $$$          '$J_q\,\,/\,\,$Wm$^{-2}$'};

% % 6-panels:
% naxs = 6;
% PosVec = [0.06    0.8475      0.7509     0.14; ...
%           0.06    0.69      0.7509    0.14; ...
%           0.06    0.5325      0.7509    0.14; ...
%           0.06    0.3750      0.7509    0.14; ...
%           0.06    0.2175      0.7509    0.14; ...
%           0.06    0.06      0.7509    0.14];

% 9-panels:
% naxs = 9;
% PosVec = [0.06    0.85      0.7509     0.085; ...
%           0.06    0.75      0.7509    0.085; ...
%           0.06    0.65     0.7509    0.085; ...
%           0.06    0.55      0.7509    0.085; ...
%           0.06    0.45      0.7509    0.085; ...
%           0.06    0.35      0.7509    0.085; ...
%           0.06    0.25      0.7509    0.085; ...
%           0.06    0.15      0.7509    0.085; ...
%           0.06    0.05      0.7509    0.085];

% 5-panels:
naxs = 5;
PosVec = [0.06    0.77      0.705     0.145; ...
          0.06    0.59      0.7509    0.145; ...
          0.06    0.41     0.7509    0.145; ...
          0.06    0.23      0.7509    0.145; ...
          0.06    0.05      0.7509    0.145];


%Combined PosVec's:
%left panels:
% $$$ PosVec = [0.06    0.8225      0.389     0.135; ...
% $$$           0.06    0.67      0.389    0.135; ...
% $$$           0.06    0.5175      0.389    0.135; ...
% $$$           0.06    0.3650      0.389    0.135; ...
% $$$           0.06    0.2125      0.389    0.135; ...
% $$$           0.06    0.06      0.389    0.135];
%right panels:
% $$$ PosVec = [0.463    0.8225      0.43     0.135; ...
% $$$           0.463    0.67      0.43    0.135; ...
% $$$           0.463    0.5175      0.43    0.135; ...
% $$$           0.463    0.3650      0.43    0.135; ...
% $$$           0.463    0.2125      0.43    0.135; ...
% $$$           0.463    0.06      0.43    0.135];
xlab = [0 0 0 0 0 0 0 0 1];
xtic = [0 0 0 0 0 0 0 0 1];
%VarOp{1} = {'dvdy','Tr','Zr'};
%VarOp{2} = {'DUDZ','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
%VarOp{3} = {'N2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ VarOp{4} = {'RSh2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
%VarOp{4} = {'Ri','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{1} = {'Forcing','Tw'};
VarOp{2} = {'Jq','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{3} = {'Jq_bl','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{4} = {'Jq_nonlocal','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{5} = {'Jq_shear','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};

caxs = [-500 500;...
        -500 500;...
        -500 500;...
        -500 500;...
        -500 500;        ];
names = {'$Q$ total surface forcing $\,\,/\,\,$Wm$^{-2}$',...
         '$J_q$ total mixing $\,\,/\,\,$Wm$^{-2}$',...
         '$J_q$ boundary layer diffusion $\,\,/\,\,$Wm$^{-2}$',...
         '$J_q$ boundary layer nonlocal flux $\,\,/\,\,$Wm$^{-2}$',...
         '$J_q$ internal diffusion $\,\,/\,\,$Wm$^{-2}$'};

% $$$ subplot('Position',PosVec(1,:));
% $$$ plot(t/86400,dvdy,'-k','LineWidth',2);
% $$$ hold on;
% $$$ grid on;
% $$$ ylabel('$dv/dy$ s$^{-1}$');
% $$$ set(gca,'xticklabel',[]);
for sp = 2:naxs

subplot('Position',PosVec(sp,:));
eval(['Var = ' VarOp{sp}{1} ';']);
eval(['X = ' VarOp{sp}{2} ';']);
eval(['Y = ' VarOp{sp}{3} ';']);
if (lnfilt ~= 0)
    Var = filter_field(Var',lnfilt,'-t')';
    Bft = filter_field(b',lnfilt,'-t')';
else
    Bft = b;
end
if (intp ~= 0)
    pcolPlot(interp2(X(:,tvec),intp),interp2(Y(:,tvec),intp),interp2(Var(:,tvec),intp));
else
    pcolPlot(X(:,tvec),Y(:,tvec),Var(:,tvec));
end
hold on;
contour(Tr(:,tvec),Zr(:,tvec),Bft(:,(tvec)),-g/rho0*[0:0.1:40],'-k'); %0:0.2:40
caxis(caxs(sp,:));
plot(t(tvec)/86400,Hsbl(tvec),'-','Color',[1 1 1],'LineWidth',2);
plot(t(tvec)/86400,EUC(tvec),'-','Color',[0.6 0.6 0.6],'LineWidth',2);
axis(axs);
ylabel('Depth (m)','FontSize',15);
% $$$ set(gca,'ytick',[]);
set(gca,'FontSize',15);
if (~xlab(sp))
    xlabel('');
else
    xlabel('Time (days)','FontSize',15);
end
if (~xtic(sp))
    set(gca,'xtick',[]);
end
title('');
% $$$ colorbar('off');
cb = colorbar;
set(cb,'FontSize',15);
% $$$ text(txtx,txty,names(sp),'FontSize',15,'BackgroundColor','w', ...
% $$$      'HorizontalAlignment','right');
text(X(1,tvec(1))+1,txty,names(sp),'FontSize',13,'BackgroundColor','w', ...
     'HorizontalAlignment','left');
xlim([X(1,tvec(1)) X(1,tvec(end))]);

meanplot = 0;
if (meanplot)
    %Plot average curves inside on left:
xlim([-7 t(end)/86400]);
set(gca,'xtick',[-6 -3.5 -1 0:5:(t(end)/86400)]);
if (~xtic(sp))
    set(gca,'xticklabel',[]);
else
    set(gca,'xticklabel',{[],[],[],0,5,10,15,20,25,30});
end
% $$$ set(gca,'xticklabel',{caxs(sp,1),mean(caxs(sp,:)),caxs(sp,2),[],[],[],[],[],[],[]});
hold on;
% $$$ plot([-6 -6],[axs(3) axs(4)],'-k');
% $$$ plot([-1 -1],[axs(3) axs(4)],'-k');
plot((mean(Var,2)-mean(caxs(sp,:)))/(caxs(sp,2)-caxs(sp,1))*5-3.5,mean(Y,2),'-k','LineWidth',2);
plot((zeros(size(mean(Var,2)))-mean(caxs(sp,:)))/(caxs(sp,2)-caxs(sp,1))*5-3.5,mean(Y,2),'-k');
end
end

sp = 1;
subplot('Position',PosVec(sp,:));
eval(['Var = ' VarOp{sp}{1} ';']);
eval(['X = ' VarOp{sp}{2} ';']);
plot(X(:,tvec),Var(:,tvec));
hold on;
yline(0)
text(X(1,tvec(1))+1,txty+170,names(sp),'FontSize',13, ...
     'HorizontalAlignment','left');
xlim([X(1,tvec(1)) X(1,tvec(end))]);
xlabel('');
set(gca,'xtick',[]);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',15);

title('');
ylabel('');


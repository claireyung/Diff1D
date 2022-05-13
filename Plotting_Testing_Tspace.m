
% Plot TIW and nTIW standard runs
fnTIWrun = 'C:/Users/holme/data/Diff1D/KPP_nDIR_nTIW_strun.mat';
fTIWrun = 'C:/Users/holme/data/Diff1D/KPP_nDIR_strun.mat';
load(fTIWrun);

%%% Calculate temperature binned heat fluxes:

% Restrict times:
% $$$ d1 = 105;
% $$$ d2 = 165;
d1 = 0;
d2 = 45;
[tmp tII] = min(abs(t/86400-d1));
[tmp tFF] = min(abs(t/86400-d2));
Nt = tFF-tII+1;

T = T(:,tII:tFF);
kt = kt(:,tII:tFF);

% Calculate J and dJdz:
DTDZ = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
J = zeros(Nz+1,Nt);
J(2:end-1,:) = -kt(2:(end-1),:).*DTDZ*Cp*rho0;
dJdz = (J(2:end,:)-J(1:(end-1),:)); % Wm-2 within each grid cell

% Define temperature grid:
dTi = 0.5;
Tie = 10:dTi:30;
Ti = avg(Tie);
NTi = length(Ti);

% Do temperature binning:
JT = zeros(NTi,Nt);
for ii = NTi:-1:1
    ii    
    for zi = 1:Nz
        inds = T(zi,:)>Tie(ii) & T(zi,:)<=Tie(ii+1);
        JT(ii,:) = JT(ii,:)+dJdz(zi,:).*(inds*1);
    end
end
JT = cat(1,zeros(1,Nt),cumsum(JT));

% Do temperature binning of mean:
Tmean = mean(T,2);
dJdzmean = mean(dJdz,2);
JTmean = zeros(NTi,1);
for ii = NTi:-1:1
    for zi = 1:Nz
        inds = Tmean(zi)>Tie(ii) & Tmean(zi)<=Tie(ii+1);
        JTmean(ii) = JTmean(ii)+dJdzmean(zi).*(inds*1);
    end
end
JTmean = cat(1,zeros(1,1),cumsum(JTmean));

% Do some plotting in T and z space:
figure;
tind = 50;
subplot(1,2,1);
plot(J(:,tind),z_w,'--k');
hold on;
plot(mean(J,2),z_w,'-k','linewidth',2);

subplot(1,2,2);
plot(JT(:,tind),Tie,'--k');
hold on;
plot(mean(JT,2),Tie,'-k','linewidth',2);
plot(JTmean,Tie,':k','linewidth',2);











%%%% Plotting time-depth series:

[Tr,Zr] = meshgrid(t/86400,z_rho);
[Tw,Zw] = meshgrid(t/86400,z_w);

%Times to plot:
d1 = 105;
d2 = 165;
[tmp tII] = min(abs(t/86400-d1));
[tmp tFF] = min(abs(t/86400-d2));
tplot = 15;
tvec = tII:tplot:tFF;

axs = [d1 d2 -250 0];
intp = 0;
lnfilt = 0;
txtx = d1 + (d2-d1)*0.995;
txty = -220;

%Derive variables:
DUDZ = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
DVDZ = (v(2:end,:)-v(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
RSh2 = DUDZ.^2+DVDZ.^2-4*N2;
Ri = 1./(N2./(DUDZ.^2+DVDZ.^2));
DTDZ = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
Jq = -kt(2:(end-1),:).*DTDZ*Cp*rho0;

figure;
set(gcf,'Position',[453 26 1006 947]);

% 6-panels:
naxs = 6;
PosVec = [0.06    0.8475      0.7509     0.14; ...
          0.06    0.69      0.7509    0.14; ...
          0.06    0.5325      0.7509    0.14; ...
          0.06    0.3750      0.7509    0.14; ...
          0.06    0.2175      0.7509    0.14; ...
          0.06    0.06      0.7509    0.14];
xlab = [0 0 0 0 0 1];
xtic = [0 0 0 0 0 1];
VarOp{1} = {'T','Tr','Zr'};
VarOp{2} = {'DTDZ','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{3} = {'N2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{4} = {'RSh2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
% $$$ VarOp{4} = {'Ri','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
VarOp{5} = {'kt','Tw','Zw'};
VarOp{6} = {'Jq','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
caxs = [0 23;...
        0 1e-1; ...
        0 3e-4; ...
        -2e-3 1e-3;...
        0 2.2e-3;...
        -500 0;];
names = {'Temperature', ...
         '$\partial$T/$\partial$z$\,\,/\,\,$s$^{-1}$', ...
         'N$^2\,\,/\,\,$s$^{-2}$', ...
         'Sh$^2_{red}\,\,/\,\,10^{-4}$s$^{-2}$',...
        '$\kappa_T\,\,/\,\,$m$^2$s$^{-1}$',...
         '$J_q\,\,/\,\,$Wm$^{-2}$'};

for sp = 1:naxs

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
contour(Tr(:,tvec),Zr(:,tvec),Bft(:,(tvec)),-g/rho0*[0:0.2:40],'-k');
caxis(caxs(sp,:));
plot(t(tvec)/86400,Hsbl(tvec),'-','Color',[1 1 1],'LineWidth',2);
% $$$ plot(t(tvec)/86400,EUC(tvec),'-','Color',[0.6 0.6 0.6],'LineWidth',2);
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
text(txtx,txty,names(sp),'FontSize',15,'BackgroundColor','w', ...
     'HorizontalAlignment','right');
xlim([X(1,tvec(1)) X(1,tvec(end))]);

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

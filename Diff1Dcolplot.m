figure;
set(gcf,'Position',[453 26 1006 947]);
[Tr,Zr] = meshgrid(t/86400,z_rho);
[Tw,Zw] = meshgrid(t/86400,z_w);
axs = [0 t(end)/86400 -200 0];
intp = 0;
lnfilt = 0;
txtx = axs(2)*0.995;
txty = axs(3)+2;

meanplot = 0; %plot mean curves on left of pcolor plots.

%Times to plot:
[tmp tII] = min(abs(t/86400-0));
[tmp tFF] = min(abs(t/86400-50000));
ntpts = 400;
if (ntpts<(tFF-tII+1))
    tvec = round(linspace(tII,tFF,ntpts));
else
    tvec = tII:tFF;
end

% $$$ %Derive variables:
DUDZ = (u(2:end,:)-u(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    length(t)]);
DVDZ = (v(2:end,:)-v(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    length(t)]);
N2 = (b(2:end,:)-b(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    length(t)]);
RiInv = 1./(N2./(DUDZ.^2+DVDZ.^2));

% 6-panels:
naxs = 6;
PosVec = [0.12    0.8475      0.7509     0.14; ...
          0.12    0.69      0.7509    0.14; ...
          0.12    0.5325      0.7509    0.14; ...
          0.12    0.3750      0.7509    0.14; ...
          0.12    0.2175      0.7509    0.14; ...
          0.12    0.06      0.7509    0.14];

xlab = [0 0 0 0 0 1];
xtic = [0 0 0 0 0 1];
caxs = [];
names = {};
%Variable 1:
VarOp{1} = {'T','Tr','Zr'};
caxs(1,:) = [21 22.25];
names{1} = '$T\,\,/\,\,^\circ$C';
%Variable 2:
VarOp{2} = {'u','Tr','Zr'};
caxs(2,:) = [-0.1 0.1];
names{2} = '$u\,\,/\,\,$ms$^{-1}$';
%Variable 3:
VarOp{3} = {'v','Tr','Zr'};
caxs(3,:) = [-0.1 0.1];
names{3} = '$v\,\,/\,\,$ms$^{-1}$';
%Variable 4:
VarOp{4} = {'N2','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
caxs(4,:) = [0 4.5e-5];
names{4} = 'N$^2\,\,/\,\,$s$^{-2}$';
%Variable 5:
VarOp{5} = {'RiInv','Tw(2:(end-1),:)','Zw(2:(end-1),:)'};
caxs(5,:) = [0 10];
names{5} = 'Ri$^{-1}$';
%Variable 6:
VarOp{6} = {'log10(kt)','Tw','Zw'};
caxs(6,:) = [-5 -1];
names{6} = '$log_{10}(\kappa_T\,\,/\,\,$m$^2$s$^{-1})$';

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
contour(Tr(:,tvec),Zr(:,tvec),Bft(:,tvec),-g/rho0*[0:0.05:40],'-k');
caxis(caxs(sp,:));
plot(t(tvec)/86400,Hsbl(tvec),'-','Color',[1 1 1],'LineWidth',2);
axis(axs);
ylabel('Depth (m)','FontSize',15);
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
cb = colorbar;
set(cb,'FontSize',15);
text(txtx,txty,names(sp),'FontSize',15,'BackgroundColor','w', ...
     'HorizontalAlignment','right');

if (meanplot)
    %Plot average curves inside on left:
xlim([-7 t(end)/86400]);
set(gca,'xtick',[-6 -3.5 -1 0:5:(t(end)/86400)]);
if (~xtic(sp))
    set(gca,'xticklabel',[]);
else
    set(gca,'xticklabel',{[],[],[],0,5,10,15,20,25,30});
end
hold on;
plot((mean(Var,2)-mean(caxs(sp,:)))/(caxs(sp,2)-caxs(sp,1))*5-3.5,mean(Y,2),'-k','LineWidth',2);
plot((zeros(size(mean(Var,2)))-mean(caxs(sp,:)))/(caxs(sp,2)-caxs(sp,1))*5-3.5,mean(Y,2),'-k');
end
end

%%%%% Scatter PLOT:
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

load('KPP_nDIR_strun.mat');

%Setup figure:
figure;
set(gcf,'Position',[71 6 1845 999]);
mnK = 0;
% $$$ mxK = 2.7e-3;
mxK = 4.1e-3;
mnDT = 0.015;
mxDT = 0.12;
xres = 50;
yres = 50;
[K,DT] = ndgrid(mnK:((mxK-mnK)/xres):mxK,mnDT:((mxDT-mnDT)/yres):mxDT);

PhaseDiag = subplot('Position',[0.1658    0.2963    0.5722    0.6457]);
contourf(K,DT,-Cp*rho0*K.*DT,[-500:50:0]);
axis([mnK mxK mnDT mxDT]);
caxis([-500 0]);
cb = colorbar('Location','South');
xlabel(cb,'$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
set(cb,'Position',[0.4439    0.3098    0.2820    0.0250]);
colormap(dark);%flipud(gray));
set(gca,'xticklabel',[]);
ylabel('$dT/dz\,\,/\,\,^\circ$Cm$^{-1}$','FontSize',25);
set(gca,'YDir','reverse');
hold on;

KDist = subplot('Position',[0.1658    0.0969    0.5724    0.1658]);
xlim([mnK mxK]);
ylim([0 4000]);
xlabel('$\kappa_T\,\,/\,\,$m$^{2}$s$^{-1}$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;

Jqplot = subplot('Position',[0.7512    0.3423    0.1854    0.5495]);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
ylabel('Depth (m)','FontSize',25);
set(gca,'yaxislocation','right');
ylim([-150 0]);
xlim([-240 0]);
box on;
grid on;
hold on;

%Setup limits:
[tmp tII] = min(abs(t/86400-135));
[tmp tFF] = min(abs(t/86400-150));
z = z_w(2:(end-1));
[tmp indtop] = min(abs(z+50));
[tmp indbot] = min(abs(z+150));

tvec = tII:(tFF-1);
% $$$ zvec = indbot:indtop;
zvec = [9 12 14 16];
% $$$ zvec = [9 13 16];
% $$$ zvec = (indbot+1):2:indtop;
linetype = {'.','+','s','x'};
dots = 720/4;
linesize = [20 20 7 20];
colors = {'k',[0.4941    0.1843    0.5569]};

hold on;
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);

%%Plot histograms:
xcount = mnK:5e-5:mxK;
hcount1 = histc(reshape(kt(zvec,:),1,[]),xcount);

%%Plot scatter plot:
ktS = kt(:,1:dots:end);dTdzS=dTdz(:,1:dots:end);
for zi = 1:length(zvec)
plot(ktS(zvec(zi),:),dTdzS(zvec(zi),:), ...
     linetype{1},'color',colors{1},'MarkerSize',linesize(1), ...
     'Parent',PhaseDiag,'LineWidth',2);
end
plot(mean(kt(zvec,:),2),mean(dTdz(zvec,:),2),...
     linetype{2},'color',colors{1},'MarkerSize',linesize(2), ...
     'LineWidth',3,'Parent',PhaseDiag);

%%Plot Jq:
plot(mean(-kt.*dTdz*Cp*rho0,2),z,'-','LineWidth',2,'color', ...
     colors{1},'Parent',Jqplot);


load('KPP_nDIR_K04em3.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);

%%Plot Histograms:
hcount2 = histc(reshape(kt(zvec,:),1,[]),xcount);
bar(xcount,hcount1,0.85,'Parent',KDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
bar(xcount,hcount2,0.5,'Parent',KDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});

%%Plot scatter plot:
ktS = kt(:,1:dots:end);dTdzS=dTdz(:,1:dots:end);
for zi = 1:length(zvec)
plot(ktS(zvec(zi),:),dTdzS(zvec(zi),:), ...
     linetype{3},'color',colors{2},'MarkerSize',linesize(3),'Parent',PhaseDiag,'LineWidth',2);
end
plot(mean(kt(zvec,:),2),mean(dTdz(zvec,:),2),...
     linetype{4},'color',colors{2},'MarkerSize',linesize(4),'LineWidth',3,'Parent',PhaseDiag);

%%Plot Jq:
plot(mean(-kt.*dTdz*Cp*rho0,2),z,'-','LineWidth',2,'color', ...
     colors{2},'Parent',Jqplot);


load('KPP_nDIR_nTIW_strun.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);

%Plot Jq:
plot(mean(-kt.*dTdz*Cp*rho0,2),z,'--','LineWidth',2,'color', ...
     colors{1},'Parent',Jqplot);

load('KPP_nDIR_nTIW_K04em3.mat');
kt = kt(2:(end-1),tvec);
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);

%Plot Jq:
plot(mean(-kt.*dTdz*Cp*rho0,2),z,'--','LineWidth',2,'color', ...
     colors{2},'Parent',Jqplot);

% $$$ %Depth levels on Jq:
% $$$ for zi = 1:length(zvec)
% $$$ plot([-400 0],z(zvec(zi))*[1 1],'--k','Parent',Jqplot);
% $$$ end

%Depth levels short on Jq:
for zi = 2:length(zvec)
plot([-50 -25],z(zvec(zi))*[1 1],'-k','LineWidth',2,'Parent', ...
     Jqplot);
text(-50+5,z(zvec(zi))+3,[num2str(round(10*z(zvec(zi)))/10) 'm'],'FontSize', ...
     15,'BackgroundColor','w','Parent',Jqplot);
end
plot([-200 -175],z(zvec(1))*[1 1],'-k','LineWidth',2,'Parent', ...
     Jqplot);
text(-200+5,z(zvec(1))+3,[num2str(round(10*z(zvec(1)))/10) 'm'],'FontSize', ...
     15,'BackgroundColor','w','Parent',Jqplot);


%%Text labels:
if (length(zvec) == 4)
PosVec = [0.0004 0.0941;...
          0.0014 0.0607;...
          0.0020 0.0423;...
          0.0020 0.0198];
else
PosVec = [0.0004 0.0941;...
          0.0016 0.0559;...
          0.0021 0.0320];
end

for zi=1:length(zvec)
text(PosVec(zi,1),PosVec(zi,2),[num2str(round(10*z(zvec(zi)))/10) 'm'],'FontSize', ...
     20,'BackgroundColor','w','Parent',PhaseDiag);
end


%%Legends:

lg = subplot('Position',[0.01 0.01 0.01 0.01]);
hold on;
plot(0,0,[linetype{1} '-'],'color',colors{1},'MarkerSize', ...
     linesize(1),'LineWidth',2,'visible','off');
plot(0,0,[linetype{3} '-'],'color',colors{2},'MarkerSize', ...
     linesize(3),'LineWidth',2,'visible','off');

plot(0,0,[linetype{2}],'color',colors{1},'MarkerSize', ...
     linesize(2),'LineWidth',3,'visible','off');
plot(0,0,[linetype{4}],'color',colors{2},'MarkerSize', ...
     linesize(4),'LineWidth',3,'visible','off');
set(lg,'visible','off');
lg = legend('KPP 50\%','KPP','KPP 50\% mean','KPP mean');
set(lg,'FontSize',25);
set(lg,'Position',[0.6188    0.4438    0.0967    0.1278]);

lg2 = subplot('Position',[0.99 0.01 0.01 0.01]);
hold on;
plot(0,0,'-k','LineWidth',2,'visible','off');
plot(0,0,'--k','LineWidth',2,'visible','off');
set(lg2,'visible','off');
lg2 = legend('TIWs','no TIW');
set(lg2,'FontSize',25);
set(lg2,'Position',[0.7594    0.3581    0.0908    0.0703]);


set(findobj(gcf,'Position',[0.4439    0.3098    0.2820    0.0250]),'xtick',[-500:100:0]);

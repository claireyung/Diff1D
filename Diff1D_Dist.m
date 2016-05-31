%%%%% Distributions PLOT:
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

fnames = {'KPP_nDIR_strun.mat', ...
          'KPP_nDIR_SYM2_strun.mat', ...
          'KPP_nDIR_nTIW_strun.mat',...
           'PP_nDIR_strun.mat',...
           'PP_nDIR_SYM2_strun.mat',...
           'PP_nDIR_nTIW_strun.mat',...
          };
load(fnames{1});

% $$$ PosVec = [0.08    0.1    0.43    0.1658;...
% $$$           0.08    0.34    0.43    0.1658;...
% $$$           0.08    0.58    0.43    0.1658;...
% $$$           0.08    0.82    0.43    0.1658;...
% $$$           0.53    0.1    0.43    0.1658;...
% $$$           0.53    0.34    0.43    0.1658;...
% $$$           0.53    0.58    0.43    0.1658;...
% $$$           0.53    0.82    0.43    0.1658];
PosVec = [0.1    0.12    0.40    0.1520;...
          0.1    0.35    0.40    0.1520;...
          0.1    0.58    0.40    0.1520;...
          0.1    0.81    0.40    0.1520;...
          0.52    0.12    0.40    0.1520;...
          0.52    0.35    0.40    0.1520;...
          0.52    0.58    0.40    0.1520;...
          0.52    0.81    0.40    0.1520];
          
%Setup figure:
figure;
set(gcf,'Position',[71 6 1845 999]);

zvec = [13];
z = z_w(2:(end-1));
z(zvec)
colors = {[0 0 1],[0 0.6 0],'k',[0 0 1],[0 0.6 0],'k'};
mnK = 0;
mxK = 2e-3;
Khlims = [0 3000];
Kxcount = mnK:((mxK-mnK)/54):mxK;
mnSh = -0.04;
mxSh = -0.005;
Shhlims = [0 1500];
Shxcount = mnSh:((mxSh-mnSh)/54):mxSh;
mnN2 = 1e-4;
mxN2 = 1.6e-4;
N2hlims = [0 4000];
N2xcount = mnN2:((mxN2-mnN2)/54):mxN2;
mnRi = 0;
mxRi = 1.25;
Rihlims = [0 2000];
Rixcount = mnRi:((mxRi-mnRi)/54):mxRi;
mnRSh2 = -7e-4;
mxRSh2 = 8e-4;
RSh2hlims = [0 2000];
RSh2xcount = mnRSh2:((mxRSh2-mnRSh2)/54):mxRSh2;
mnJq = 0;
mxJq = 400;
Jqhlims = [0 3000];
Jqxcount = mnJq:((mxJq-mnJq)/54):mxJq;

KDist = subplot('Position',PosVec(2,:));
xlim([mnK-0.03e-3 mxK]);
ylim(Khlims);
xlabel('$\kappa_T\,\,/\,\,10^{-3}$m$^{2}$s$^{-1}$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;
set(gca,'ytick',[0:1000:3000]);
set(gca,'xtick',[0:0.5:2]*1e-3);
set(gca,'xticklabel',[0:0.5:2]);

ShDist = subplot('Position',PosVec(4,:));
xlim([mnSh mxSh]);
ylim(Shhlims);
xlabel('$\partial u/\partial z\,\,/\,\,$s$^{-1}$','FontSize',25);
ylabel('Histogram');
title('KPP Interior Mixing');
box on;
grid on;
hold on;

% $$$ N2Dist = subplot('Position',[0.08    0.58   0.43    0.1658]);
% $$$ xlim([mnN2 mxN2]);
% $$$ ylim(N2hlims);
% $$$ xlabel('$N^2\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ ylabel('Histogram');
% $$$ box on;
% $$$ grid on;
% $$$ hold on;
RiDist = subplot('Position',PosVec(3,:));
xlim([mnRi mxRi]);
ylim(Rihlims);
xlabel('$Ri$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;
% $$$ RSh2Dist = subplot('Position',PosVec(3,:));
% $$$ xlim([mnRSh2 mxRSh2]);
% $$$ ylim(RSh2hlims);
% $$$ xlabel('$Sh^2_{red}\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ ylabel('Histogram');
% $$$ box on;
% $$$ grid on;
% $$$ hold on;

JqDist = subplot('Position',PosVec(1,:));
xlim([mnJq-5 mxJq-20]);
ylim(Jqhlims);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;
set(gca,'ytick',[0:1000:3000]);
set(gca,'xtick',[0:100:300]);

load(fnames{1});
%Setup limits:
[tmp tII] = min(abs(t/86400-135));
[tmp tFF] = min(abs(t/86400-150));
tvec = tII:(tFF-1);

hold on;
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount1 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount1 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount1 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount1 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount1 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount1 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew1 = skewness(reshape(kt(zvec,:),1,[]));
Shskew1 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew1 = skewness(reshape(N2(zvec,:),1,[]));
Riskew1 = skewness(reshape(Ri(zvec,:),1,[]));
RSh2skew1 = skewness(reshape(RSh2(zvec,:),1,[]));
Jqskew1 = skewness(reshape(Jq(zvec,:),1,[]));
Khmean1 = mean(reshape(kt(zvec,:),1,[]));
Shmean1 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean1 = mean(reshape(N2(zvec,:),1,[]));
Rimean1 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean1 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean1 = mean(reshape(Jq(zvec,:),1,[]));

load(fnames{2});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount2 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount2 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount2 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount2 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount2 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount2 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew2 = skewness(reshape(kt(zvec,:),1,[]));
Shskew2 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew2 = skewness(reshape(N2(zvec,:),1,[]));
Riskew2 = skewness(reshape(Ri(zvec,:),1,[]));
RSh2skew2 = skewness(reshape(RSh2(zvec,:),1,[]));
Jqskew2 = skewness(reshape(Jq(zvec,:),1,[]));
Khmean2 = mean(reshape(kt(zvec,:),1,[]));
Shmean2 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean2 = mean(reshape(N2(zvec,:),1,[]));
Rimean2 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean2 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean2 = mean(reshape(Jq(zvec,:),1,[]));

load(fnames{3});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount3 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount3 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount3 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount3 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount3 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount3 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khmean3 = mean(reshape(kt(zvec,:),1,[]));
Shmean3 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean3 = mean(reshape(N2(zvec,:),1,[]));
Rimean3 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean3 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean3 = mean(reshape(Jq(zvec,:),1,[]));

bar(Kxcount,Khcount1,0.85,'Parent',KDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Shxcount,Shhcount1,0.85,'Parent',ShDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
% $$$ bar(N2xcount,N2hcount1,0.85,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Rixcount,Rihcount1,0.85,'Parent',RiDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
% $$$ bar(RSh2xcount,RSh2hcount1,0.85,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Jqxcount,Jqhcount1,0.85,'Parent',JqDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Kxcount,Khcount2,0.5,'Parent',KDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});
bar(Shxcount,Shhcount2,0.5,'Parent',ShDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});
% $$$ bar(N2xcount,N2hcount2,0.5,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{2},'EdgeColor',colors{2});
bar(Rixcount,Rihcount2,0.5,'Parent',RiDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});
% $$$ bar(RSh2xcount,RSh2hcount2,0.5,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{2},'EdgeColor',colors{2});
bar(Jqxcount,Jqhcount2,0.5,'Parent',JqDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});
% $$$ bar(Kxcount,Khcount3,0.15,'Parent',KDist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(Shxcount,Shhcount3,0.15,'Parent',ShDist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(N2xcount,N2hcount3,0.15,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(Rixcount,Rihcount3,0.15,'Parent',RiDist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(RSh2xcount,RSh2hcount3,0.15,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(Jqxcount,Jqhcount3,0.15,'Parent',JqDist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});

%Skewnesses:
xlims = get(KDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(KDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Khskew1)/100),'color',colors{1},'Parent',KDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Khskew2)/100),'color',colors{2},...
     'Parent',KDist,'HorizontalAlignment','center');
xlims = get(ShDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(ShDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Shskew1)/100),'color',colors{1},'Parent',ShDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Shskew2)/100),'color',colors{2},...
     'Parent',ShDist,'HorizontalAlignment','center');
% $$$ xlims = get(N2Dist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
% $$$ ylims = get(N2Dist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
% $$$ text(txtx,txty,num2str(round(100*N2skew1)/100),'color',colors{1},'Parent',N2Dist,'HorizontalAlignment','center');
% $$$ text(txtx,txty-txty*0.2,num2str(round(100*N2skew2)/100),'color',colors{2},...
% $$$      'Parent',N2Dist,'HorizontalAlignment','center');
xlims = get(RiDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(RiDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Riskew1)/100),'color',colors{1},'Parent',RiDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Riskew2)/100),'color',colors{2},...
     'Parent',RiDist,'HorizontalAlignment','center');
% $$$ xlims = get(RSh2Dist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
% $$$ ylims = get(RSh2Dist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
% $$$ text(txtx,txty,num2str(round(100*RSh2skew1)/100),'color',colors{1},'Parent',RSh2Dist,'HorizontalAlignment','center');
% $$$ text(txtx,txty-txty*0.2,num2str(round(100*RSh2skew2)/100),'color',colors{2},...
% $$$      'Parent',RSh2Dist,'HorizontalAlignment','center');
xlims = get(JqDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(JqDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Jqskew1)/100),'color',colors{1},'Parent',JqDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Jqskew2)/100),'color',colors{2},...
     'Parent',JqDist,'HorizontalAlignment','center');

%Means:
xlims = get(KDist,'xlim');ylims = get(KDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Khmean1-width Khmean1 Khmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',KDist,'LineWidth',3);
line([Khmean2-width Khmean2 Khmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',KDist,'LineWidth',3);
line([Khmean3-width Khmean3 Khmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',KDist,'LineWidth',3);
xlims = get(ShDist,'xlim');ylims = get(ShDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Shmean1-width Shmean1 Shmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',ShDist,'LineWidth',3);
line([Shmean2-width Shmean2 Shmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',ShDist,'LineWidth',3);
line([Shmean3-width Shmean3 Shmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',ShDist,'LineWidth',3);
xlims = get(RiDist,'xlim');ylims = get(RiDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Rimean1-width Rimean1 Rimean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',RiDist,'LineWidth',3);
line([Rimean2-width Rimean2 Rimean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',RiDist,'LineWidth',3);
line([Rimean3-width Rimean3 Rimean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',RiDist,'LineWidth',3);
% $$$ xlims = get(RSh2Dist,'xlim');ylims = get(RSh2Dist,'ylim');
% $$$ basey = ylims(2)-0.04*(ylims(2)-ylims(1));
% $$$ width = (xlims(2)-xlims(1))*0.0175;
% $$$ height = (ylims(2)-ylims(1))*0.0525;
% $$$ line([RSh2mean1-width RSh2mean1 RSh2mean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',RSh2Dist,'LineWidth',3);
% $$$ line([RSh2mean2-width RSh2mean2 RSh2mean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',RSh2Dist,'LineWidth',3);
% $$$ line([RSh2mean3-width RSh2mean3 RSh2mean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',RSh2Dist,'LineWidth',3);
xlims = get(JqDist,'xlim');ylims = get(JqDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Jqmean1-width Jqmean1 Jqmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',JqDist,'LineWidth',3);
line([Jqmean2-width Jqmean2 Jqmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',JqDist,'LineWidth',3);
line([Jqmean3-width Jqmean3 Jqmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',JqDist,'LineWidth',3);

%Legends:
% $$$ axes(KDist);
% $$$ lg = legend('KPP regular TIWs','KPP symmetric TIWs','KPP no TIWs');
% $$$ set(lg,'Position',[0.1804 0.4338 0.0604 0.0510]);

%%% Sym simulation:
KDist = subplot('Position',PosVec(6,:));
xlim([mnK-0.03e-3 mxK]);
ylim(Khlims);
set(gca,'yticklabel',[]);
xlabel('$\kappa_T\,\,/\,\,$m$^{2}$s$^{-1}$','FontSize',25);
% $$$ ylabel('Histogram');
box on;
grid on;
hold on;
set(gca,'ytick',[0:1000:3000]);
set(gca,'xtick',[0:0.5:2]*1e-3);
set(gca,'xticklabel',[0:0.5:2]);

ShDist = subplot('Position',PosVec(8,:));
xlim([mnSh mxSh]);
ylim(Shhlims);
set(gca,'yticklabel',[]);
xlabel('$\partial u/\partial z\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ ylabel('Histogram');
title('P\&P Interior Mixing');
box on;
grid on;
hold on;

% $$$ N2Dist = subplot('Position',[0.53    0.58   0.43    0.1658]);
% $$$ xlim([mnN2 mxN2]);
% $$$ ylim(N2hlims);
% $$$ set(gca,'yticklabel',[]);
% $$$ xlabel('$N^2\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ box on;
% $$$ grid on;
% $$$ hold on;
RiDist = subplot('Position',PosVec(7,:));
xlim([mnRi mxRi]);
ylim(Rihlims);
set(gca,'yticklabel',[]);
xlabel('$Ri$','FontSize',25);
box on;
grid on;
hold on;
% $$$ RSh2Dist = subplot('Position',PosVec(7,:));
% $$$ xlim([mnRSh2 mxRSh2]);
% $$$ ylim(RSh2hlims);
% $$$ set(gca,'yticklabel',[]);
% $$$ xlabel('$Sh^2_{red}\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ box on;
% $$$ grid on;
% $$$ hold on;

JqDist = subplot('Position',PosVec(5,:));
xlim([mnJq-5 mxJq-20]);
ylim(Jqhlims);
set(gca,'yticklabel',[]);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
box on;
grid on;
hold on;
set(gca,'ytick',[0:1000:3000]);
set(gca,'xtick',[0:100:300]);

load(fnames{4});
hold on;
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount1 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount1 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount1 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount1 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount1 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount1 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew1 = skewness(reshape(kt(zvec,:),1,[]));
Shskew1 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew1 = skewness(reshape(N2(zvec,:),1,[]));
Riskew1 = skewness(reshape(Ri(zvec,:),1,[]));
RSh2skew1 = skewness(reshape(RSh2(zvec,:),1,[]));
Jqskew1 = skewness(reshape(Jq(zvec,:),1,[]));
Khmean1 = mean(reshape(kt(zvec,:),1,[]));
Shmean1 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean1 = mean(reshape(N2(zvec,:),1,[]));
Rimean1 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean1 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean1 = mean(reshape(Jq(zvec,:),1,[]));

load(fnames{5});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount2 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount2 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount2 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount2 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount2 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount2 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew2 = skewness(reshape(kt(zvec,:),1,[]));
Shskew2 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew2 = skewness(reshape(N2(zvec,:),1,[]));
Riskew2 = skewness(reshape(Ri(zvec,:),1,[]));
RSh2skew2 = skewness(reshape(RSh2(zvec,:),1,[]));
Jqskew2 = skewness(reshape(Jq(zvec,:),1,[]));
Khmean2 = mean(reshape(kt(zvec,:),1,[]));
Shmean2 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean2 = mean(reshape(N2(zvec,:),1,[]));
Rimean2 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean2 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean2 = mean(reshape(Jq(zvec,:),1,[]));

load(fnames{6});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = kt.*dTdz*Cp*rho0;
%%Testing Jq/mean(dTdz):
%kt = Jq./repmat(mean(dTdz,2),[1 length(dTdz(1,:))])/Cp/rho0;
Ri = N2./(Sh2.^2);
RSh2 = Sh2.^2-4*N2;

Khcount3 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount3 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount3 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount3 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
RSh2hcount3 = histc(reshape(RSh2(zvec,:),1,[]),RSh2xcount);
size(find(Ri<0.2))
Jqhcount3 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khmean3 = mean(reshape(kt(zvec,:),1,[]));
Shmean3 = mean(reshape(Sh2(zvec,:),1,[]));
N2mean3 = mean(reshape(N2(zvec,:),1,[]));
Rimean3 = mean(reshape(Ri(zvec,:),1,[]));
RSh2mean3 = mean(reshape(RSh2(zvec,:),1,[]));
Jqmean3 = mean(reshape(Jq(zvec,:),1,[]));

bar(Kxcount,Khcount1,0.85,'Parent',KDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Shxcount,Shhcount1,0.85,'Parent',ShDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
% $$$ bar(N2xcount,N2hcount1,0.85,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Rixcount,Rihcount1,0.85,'Parent',RiDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
% $$$ bar(RSh2xcount,RSh2hcount1,0.85,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Jqxcount,Jqhcount1,0.85,'Parent',JqDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Kxcount,Khcount2,0.5,'Parent',KDist, ...
    'FaceColor',colors{5},'EdgeColor',colors{5});
bar(Shxcount,Shhcount2,0.5,'Parent',ShDist, ...
    'FaceColor',colors{5},'EdgeColor',colors{5});
% $$$ bar(N2xcount,N2hcount2,0.5,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{5},'EdgeColor',colors{5});
bar(Rixcount,Rihcount2,0.5,'Parent',RiDist, ...
    'FaceColor',colors{5},'EdgeColor',colors{5});
% $$$ bar(RSh2xcount,RSh2hcount2,0.5,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{5},'EdgeColor',colors{5});
bar(Jqxcount,Jqhcount2,0.5,'Parent',JqDist, ...
    'FaceColor',colors{5},'EdgeColor',colors{5});

% $$$ bar(Kxcount,Khcount3,0.15,'Parent',KDist, ...
% $$$     'FaceColor',colors{6},'EdgeColor',colors{6});
% $$$ bar(Shxcount,Shhcount3,0.15,'Parent',ShDist, ...
% $$$     'FaceColor',colors{6},'EdgeColor',colors{6});
% $$$ % $$$ bar(N2xcount,N2hcount3,0.15,'Parent',N2Dist, ...
% $$$ % $$$     'FaceColor',colors{6},'EdgeColor',colors{6});
% $$$ % $$$ bar(Rixcount,Rihcount3,0.15,'Parent',RiDist, ...
% $$$ % $$$     'FaceColor',colors{6},'EdgeColor',colors{6});
% $$$ bar(RSh2xcount,RSh2hcount3,0.15,'Parent',RSh2Dist, ...
% $$$     'FaceColor',colors{6},'EdgeColor',colors{6});
% $$$ bar(Jqxcount,Jqhcount3,0.15,'Parent',JqDist, ...
% $$$     'FaceColor',colors{6},'EdgeColor',colors{6});

%Skewnesses:
xlims = get(KDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(KDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Khskew1)/100),'color',colors{4},'Parent',KDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Khskew2)/100),'color',colors{5},...
     'Parent',KDist,'HorizontalAlignment','center');
xlims = get(ShDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(ShDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Shskew1)/100),'color',colors{4},'Parent',ShDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Shskew2)/100),'color',colors{5},...
     'Parent',ShDist,'HorizontalAlignment','center');
% $$$ xlims = get(N2Dist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
% $$$ ylims = get(N2Dist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
% $$$ text(txtx,txty,num2str(round(100*N2skew1)/100),'color',colors{4},'Parent',N2Dist,'HorizontalAlignment','center');
% $$$ text(txtx,txty-txty*0.2,num2str(round(100*N2skew2)/100),'color',colors{5},...
% $$$      'Parent',N2Dist,'HorizontalAlignment','center');
xlims = get(RiDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(RiDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Riskew1)/100),'color',colors{4},'Parent',RiDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Riskew2)/100),'color',colors{5},...
     'Parent',RiDist,'HorizontalAlignment','center');
% $$$ xlims = get(RSh2Dist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
% $$$ ylims = get(RSh2Dist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
% $$$ text(txtx,txty,num2str(round(100*RSh2skew1)/100),'color',colors{4},'Parent',RSh2Dist,'HorizontalAlignment','center');
% $$$ text(txtx,txty-txty*0.2,num2str(round(100*RSh2skew2)/100),'color',colors{5},...
% $$$      'Parent',RSh2Dist,'HorizontalAlignment','center');
xlims = get(JqDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(JqDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Jqskew1)/100),'color',colors{4},'Parent',JqDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Jqskew2)/100),'color',colors{5},...
     'Parent',JqDist,'HorizontalAlignment','center');

%Means:
xlims = get(KDist,'xlim');ylims = get(KDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Khmean1-width Khmean1 Khmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',KDist,'LineWidth',3);
line([Khmean2-width Khmean2 Khmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',KDist,'LineWidth',3);
line([Khmean3-width Khmean3 Khmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',KDist,'LineWidth',3);
xlims = get(ShDist,'xlim');ylims = get(ShDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Shmean1-width Shmean1 Shmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',ShDist,'LineWidth',3);
line([Shmean2-width Shmean2 Shmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',ShDist,'LineWidth',3);
line([Shmean3-width Shmean3 Shmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',ShDist,'LineWidth',3);
xlims = get(RiDist,'xlim');ylims = get(RiDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Rimean1-width Rimean1 Rimean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',RiDist,'LineWidth',3);
line([Rimean2-width Rimean2 Rimean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',RiDist,'LineWidth',3);
line([Rimean3-width Rimean3 Rimean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',RiDist,'LineWidth',3);
% $$$ xlims = get(RSh2Dist,'xlim');ylims = get(RSh2Dist,'ylim');
% $$$ basey = ylims(2)-0.04*(ylims(2)-ylims(1));
% $$$ width = (xlims(2)-xlims(1))*0.0175;
% $$$ height = (ylims(2)-ylims(1))*0.0525;
% $$$ line([RSh2mean1-width RSh2mean1 RSh2mean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',RSh2Dist,'LineWidth',3);
% $$$ line([RSh2mean2-width RSh2mean2 RSh2mean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',RSh2Dist,'LineWidth',3);
% $$$ line([RSh2mean3-width RSh2mean3 RSh2mean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',RSh2Dist,'LineWidth',3);
xlims = get(JqDist,'xlim');ylims = get(JqDist,'ylim');
basey = ylims(2)-0.04*(ylims(2)-ylims(1));
width = (xlims(2)-xlims(1))*0.0175;
height = (ylims(2)-ylims(1))*0.0525;
line([Jqmean1-width Jqmean1 Jqmean1+width],[basey-height basey basey-height],'color',colors{1},'Parent',JqDist,'LineWidth',3);
line([Jqmean2-width Jqmean2 Jqmean2+width],[basey-height basey basey-height],'color',colors{2},'Parent',JqDist,'LineWidth',3);
line([Jqmean3-width Jqmean3 Jqmean3+width],[basey-height basey basey-height],'color',colors{3},'Parent',JqDist,'LineWidth',3);

%Legends:
axes(ShDist);
% $$$ lg = legend('KPP SYM','P\&P SYM');
% $$$ lg = legend('P\&P regular TIWs','P\&P symmetric TIWs','P\&P no TIWs');
lg = legend('non-linear TIW stretching','linear TIW stretching',['no ' ...
                    'TIW stretching']);
set(lg,'Position',[0.5341    0.8919+0.0    0.1314    0.0866]);

%%Time series:
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ mnSh = -0.07;
% $$$ mxSh = -0.01;
% $$$ Shxcount = mnSh:((mxSh-mnSh)/54):mxSh;
% $$$ 
% $$$ subplot('Position',[0.1300    0.1100    0.5589    0.8150]);
% $$$ load KPP_nDIR_strun.mat;
% $$$ Sh2 = ((u(2:end,:)-u(1:(end-1),:))./...
% $$$     repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(u(1,:))]));
% $$$ plot(t(1:100:end),Sh2(16,1:100:end),'-k','LineWidth',2);
% $$$ Sh2Id = Sh2(16,1)*exp(-2.8e-6*(cos(2*pi/period* ...
% $$$                                                   t(1:100:end))-1)*period/2/pi);
% $$$ hold on; plot(t(1:100:end),Sh2Id,'-r','LineWidth',2);
% $$$ Sh2IdnA = -2.8e-6*(cos(2*pi/period*t(1:100:end))- ...
% $$$                                     1)*period/2/pi*Sh2(16,1)+Sh2(16,1);
% $$$ hold on; plot(t(1:100:end),Sh2IdnA,'-b','LineWidth',2);
% $$$ Sh2KPP = Sh2(16,1:100:end);
% $$$ load PP_nDIR_strun.mat;
% $$$ Sh2 = ((u(2:end,:)-u(1:(end-1),:))./...
% $$$     repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(u(1,:))]));
% $$$ plot(t(1:100:end),Sh2(16,1:100:end),'--k','LineWidth',2);
% $$$ load KPP_nDIR_nPha_strun.mat;
% $$$ Sh2 = ((u(2:end,:)-u(1:(end-1),:))./...
% $$$     repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(u(1,:))]));
% $$$ plot(t(1:100:end),Sh2(16,1:100:end),':k','LineWidth',2);
% $$$ legend('Diff1D KPP nDIR','Ideal','Ideal noAsym','Diff1D PP nDIR',...
% $$$        'Diff1D KPP nDIR nPha');
% $$$ xlabel('time');
% $$$ ylabel('$\partial u/\partial z$');
% $$$ 
% $$$ subplot('Position',[0.75    0.1100    0.1    0.8150]);
% $$$ Shcount = histc(reshape(Sh2KPP,1,[]),Shxcount);
% $$$ barh(Shxcount,Shcount,0.85, ...
% $$$     'FaceColor','k','EdgeColor','k');
% $$$ hold on;
% $$$ Shcount = histc(reshape(Sh2Id,1,[]),Shxcount);
% $$$ barh(Shxcount,Shcount,0.5, ...
% $$$     'FaceColor','r','EdgeColor','r');
% $$$ Shcount = histc(reshape(Sh2IdnA,1,[]),Shxcount);
% $$$ barh(Shxcount,Shcount,0.2, ...
% $$$     'FaceColor','b','EdgeColor','b');
% $$$ 

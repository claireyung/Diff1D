%%%%% Distributions PLOT:
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% $$$ fnames = {'KPP_nDIR_strun.mat', ...
% $$$          'PP_nDIR_strun.mat', ...
% $$$          'KPP_nDIR_SYM_strun.mat',...
% $$$          'PP_nDIR_SYM_strun.mat'};
fnames = {'KPP_nDIR_strun.mat', ...
         'KPP_nDIR_nTIW_strun.mat', ...
         'PP_nDIR_strun.mat',...
         'PP_nDIR_nTIW_strun.mat'};
% $$$          'PP1_nDIR_strun.mat',...
% $$$          'PP1.2_nDIR_strun.mat'};
fnamesnTIW = fnames;
for f=1:length(fnamesnTIW)
    fnamesnTIW{f} = strrep(fnamesnTIW{f},'strun','nTIW_strun');
end

%Setup figure:
figure;
set(gcf,'Position',[71 6 1845 999]);

zvec = [10:16];
z = z_w(2:(end-1));
z(zvec)
colors = {'k',[0.4941    0.1843    0.5569],[0 0 0.6],[0 0.6 0]};
mnK = 0;
mxK = 4e-3;
Khlims = [0 8000];
Kxcount = mnK:((mxK-mnK)/54):mxK;
mnSh = -0.045;
mxSh = -0.005;
Shhlims = [0 6000];
Shxcount = mnSh:((mxSh-mnSh)/54):mxSh;
mnN2 = 1e-4;
mxN2 = 1.6e-4;
N2hlims = [0 4000];
N2xcount = mnN2:((mxN2-mnN2)/54):mxN2;
mnRi = 0;
mxRi = 2;
Rihlims = [0 11000];
Rixcount = mnRi:((mxRi-mnRi)/54):mxRi;
mnJq = -600;
mxJq = 0;
Jqhlims = [0 5000];
Jqxcount = mnJq:((mxJq-mnJq)/54):mxJq;

KDist = subplot('Position',[0.08    0.1    0.43    0.1658]);
xlim([mnK mxK]);
ylim(Khlims);
xlabel('$\kappa_T\,\,/\,\,$m$^{2}$s$^{-1}$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;

ShDist = subplot('Position',[0.08    0.34    0.43    0.1658]);
xlim([mnSh mxSh]);
ylim(Shhlims);
xlabel('$\partial u/\partial z\,\,/\,\,$s$^{-2}$','FontSize',25);
ylabel('Histogram');
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
RiDist = subplot('Position',[0.08    0.58   0.43    0.1658]);
xlim([mnRi mxRi]);
ylim(Rihlims);
xlabel('$Ri$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;

JqDist = subplot('Position',[0.08    0.82  0.43    0.1658]);
xlim([mnJq mxJq]);
ylim(Jqhlims);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
ylabel('Histogram');
box on;
grid on;
hold on;

load(fnames{1});
%Setup limits:
[tmp tII] = min(abs(t/86400-135));
[tmp tFF] = min(abs(t/86400-150));
z = z_w(2:(end-1));
[tmp indtop] = min(abs(z+50));
[tmp indbot] = min(abs(z+150));

tvec = tII:(tFF-1);

hold on;
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = -kt.*dTdz*Cp*rho0;
Ri = N2./(Sh2.^2);

Khcount1 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount1 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount1 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount1 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
size(find(Ri<0.2))
Jqhcount1 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew1 = skewness(reshape(kt(zvec,:),1,[]));
Shskew1 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew1 = skewness(reshape(N2(zvec,:),1,[]));
Riskew1 = skewness(reshape(Ri(zvec,:),1,[]));
Jqskew1 = skewness(reshape(Jq(zvec,:),1,[]));

load(fnames{2});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = -kt.*dTdz*Cp*rho0;
Ri = N2./(Sh2.^2);

Khcount2 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount2 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount2 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount2 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
size(find(Ri<0.2))
Jqhcount2 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew2 = skewness(reshape(kt(zvec,:),1,[]));
Shskew2 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew2 = skewness(reshape(N2(zvec,:),1,[]));
Riskew2 = skewness(reshape(Ri(zvec,:),1,[]));
Jqskew2 = skewness(reshape(Jq(zvec,:),1,[]));

bar(Kxcount,Khcount1,0.85,'Parent',KDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Shxcount,Shhcount1,0.85,'Parent',ShDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
% $$$ bar(N2xcount,N2hcount1,0.85,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{1},'EdgeColor',colors{1});
bar(Rixcount,Rihcount1,0.85,'Parent',RiDist, ...
    'FaceColor',colors{1},'EdgeColor',colors{1});
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
bar(Jqxcount,Jqhcount2,0.5,'Parent',JqDist, ...
    'FaceColor',colors{2},'EdgeColor',colors{2});

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
xlims = get(JqDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(JqDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Jqskew1)/100),'color',colors{1},'Parent',JqDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Jqskew2)/100),'color',colors{2},...
     'Parent',JqDist,'HorizontalAlignment','center');

%Legends:
axes(RiDist);
lg = legend('KPP','P\&P');
set(lg,'Position',[0.1804 0.6738 0.0604 0.0510]);

%%% Sym simulation:
KDist = subplot('Position',[0.53    0.1    0.43    0.1658]);
xlim([mnK mxK]);
ylim(Khlims);
set(gca,'yticklabel',[]);
xlabel('$\kappa_T\,\,/\,\,$m$^{2}$s$^{-1}$','FontSize',25);
% $$$ ylabel('Histogram');
box on;
grid on;
hold on;

ShDist = subplot('Position',[0.53    0.34    0.43    0.1658]);
xlim([mnSh mxSh]);
ylim(Shhlims);
set(gca,'yticklabel',[]);
xlabel('$\partial u/\partial z\,\,/\,\,$s$^{-2}$','FontSize',25);
% $$$ ylabel('Histogram');
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
RiDist = subplot('Position',[0.53    0.58   0.43    0.1658]);
xlim([mnRi mxRi]);
ylim(Rihlims);
set(gca,'yticklabel',[]);
xlabel('$Ri$','FontSize',25);
box on;
grid on;
hold on;

JqDist = subplot('Position',[0.53    0.82  0.43    0.1658]);
xlim([mnJq mxJq]);
ylim(Jqhlims);
set(gca,'yticklabel',[]);
xlabel('$J_q\,\,/\,\,$Wm$^{-2}$','FontSize',25);
box on;
grid on;
hold on;

load(fnames{3});
hold on;
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = -kt.*dTdz*Cp*rho0;
Ri = N2./(Sh2.^2);

Khcount1 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount1 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount1 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount1 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
size(find(Ri<0.2))
Jqhcount1 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew1 = skewness(reshape(kt(zvec,:),1,[]));
Shskew1 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew1 = skewness(reshape(N2(zvec,:),1,[]));
Riskew1 = skewness(reshape(Ri(zvec,:),1,[]));
Jqskew1 = skewness(reshape(Jq(zvec,:),1,[]));

load(fnames{4});
kt = kt(2:(end-1),tvec);
N2 = (b(2:end,tvec)-b(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Sh2 = ((u(2:end,tvec)-u(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]));
dTdz = (T(2:end,tvec)-T(1:(end-1),tvec))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(tvec)]);
Jq = -kt.*dTdz*Cp*rho0;
Ri = N2./(Sh2.^2);

Khcount2 = histc(reshape(kt(zvec,:),1,[]),Kxcount);
Shhcount2 = histc(reshape(Sh2(zvec,:),1,[]),Shxcount);
N2hcount2 = histc(reshape(N2(zvec,:),1,[]),N2xcount);
Rihcount2 = histc(reshape(Ri(zvec,:),1,[]),Rixcount);
size(find(Ri<0.2))
Jqhcount2 = histc(reshape(Jq(zvec,:),1,[]),Jqxcount);
Khskew2 = skewness(reshape(kt(zvec,:),1,[]));
Shskew2 = skewness(reshape(Sh2(zvec,:),1,[]));
N2skew2 = skewness(reshape(N2(zvec,:),1,[]));
Riskew2 = skewness(reshape(Ri(zvec,:),1,[]));
Jqskew2 = skewness(reshape(Jq(zvec,:),1,[]));

bar(Kxcount,Khcount1,0.85,'Parent',KDist, ...
    'FaceColor',colors{3},'EdgeColor',colors{3});
bar(Shxcount,Shhcount1,0.85,'Parent',ShDist, ...
    'FaceColor',colors{3},'EdgeColor',colors{3});
% $$$ bar(N2xcount,N2hcount1,0.85,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{3},'EdgeColor',colors{3});
bar(Rixcount,Rihcount1,0.85,'Parent',RiDist, ...
    'FaceColor',colors{3},'EdgeColor',colors{3});
bar(Jqxcount,Jqhcount1,0.85,'Parent',JqDist, ...
    'FaceColor',colors{3},'EdgeColor',colors{3});
bar(Kxcount,Khcount2,0.5,'Parent',KDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Shxcount,Shhcount2,0.5,'Parent',ShDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
% $$$ bar(N2xcount,N2hcount2,0.5,'Parent',N2Dist, ...
% $$$     'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Rixcount,Rihcount2,0.5,'Parent',RiDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});
bar(Jqxcount,Jqhcount2,0.5,'Parent',JqDist, ...
    'FaceColor',colors{4},'EdgeColor',colors{4});

%Skewnesses:
xlims = get(KDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(KDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Khskew1)/100),'color',colors{3},'Parent',KDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Khskew2)/100),'color',colors{4},...
     'Parent',KDist,'HorizontalAlignment','center');
xlims = get(ShDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(ShDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Shskew1)/100),'color',colors{3},'Parent',ShDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Shskew2)/100),'color',colors{4},...
     'Parent',ShDist,'HorizontalAlignment','center');
% $$$ xlims = get(N2Dist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
% $$$ ylims = get(N2Dist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
% $$$ text(txtx,txty,num2str(round(100*N2skew1)/100),'color',colors{3},'Parent',N2Dist,'HorizontalAlignment','center');
% $$$ text(txtx,txty-txty*0.2,num2str(round(100*N2skew2)/100),'color',colors{4},...
% $$$      'Parent',N2Dist,'HorizontalAlignment','center');
xlims = get(RiDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(RiDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Riskew1)/100),'color',colors{3},'Parent',RiDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Riskew2)/100),'color',colors{4},...
     'Parent',RiDist,'HorizontalAlignment','center');
xlims = get(JqDist,'xlim');txtx=xlims(1)+(xlims(2)-xlims(1))*0.9;
ylims = get(JqDist,'ylim');txty=ylims(1)+(ylims(2)-ylims(1))*0.9;
text(txtx,txty,num2str(round(100*Jqskew1)/100),'color',colors{3},'Parent',JqDist,'HorizontalAlignment','center');
text(txtx,txty-txty*0.2,num2str(round(100*Jqskew2)/100),'color',colors{4},...
     'Parent',JqDist,'HorizontalAlignment','center');

%Legends:
axes(RiDist);
% $$$ lg = legend('KPP SYM','P\&P SYM');
lg = legend('P\&P1','P\&P1.2');
set(lg,'Position',[0.6304 0.6738 0.0604 0.0510]);

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



%%% Check convergence of Jq plots:
load('KPP_nDIR_nPha_lng_strun.mat');
z = z_w(2:(end-1));
kt = kt(2:(end-1),:);
dTdz = (T(2:end,:)-T(1:(end-1),:))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(T(1,:))]);
Jq_nPha = -kt.*dTdz*Cp*rho0;

%divide into 15-day intervals and average:
[tmp tII] = min(abs(t/86400));
[tmp tFF] = min(abs(t/86400-15));
d15i = tFF-tII;
Jq_nPha15 = zeros(length(z),floor(length(t)/d15i));
ind = 1;
for ti = 1:length(Jq_nPha15(1,:))
    Jq_nPha15(:,ti) = mean(Jq_nPha(:,ind:(ind+d15i-1)),2);
    ind = ind+d15i;
end


load('KPP_nDIR_lng_strun.mat');
kt = kt(2:(end-1),:);
dTdz = (T(2:end,:)-T(1:(end-1),:))./...
    repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 length(T(1,:))]);
Jq = -kt.*dTdz*Cp*rho0;
Jq15 = zeros(length(z),floor(length(t)/d15i));
ind = 1;
for ti = 1:length(Jq_nPha15(1,:))
    Jq15(:,ti) = mean(Jq(:,ind:(ind+d15i-1)),2);
    ind = ind+d15i;
end

thres = 50;
pdiff = abs(((Jq15-Jq_nPha15)./Jq15))*100;
pdiff(abs(Jq_nPha15)<thres)=0;

thres = 10;
diff = abs(Jq15-Jq_nPha15);

[X,Y] = ndgrid(15:15:max(t/86400),z);
figure;
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);
pcolPlot(X,Y,Jq15');
caxis([-250 0]);
ylim([-250 0]);
colorbar;
xlabel('Time (days)');
ylabel('Depth (m)');
title('$J_q$ 15-day means $\alpha>0$');
subplot(2,2,2);
pcolPlot(X,Y,Jq_nPha15');
caxis([-250 0]);
ylim([-250 0]);
colorbar;
xlabel('Time (days)');
ylabel('Depth (m)');
title('$J_q$ 15-day means $\alpha<0$');
subplot(2,2,3);
% $$$ pcolPlot(X,Y,pdiff');
% $$$ hold on;
% $$$ contour(X,Y,pdiff',[5 5],'-m');
% $$$ caxis([0 25]);
% $$$ title(['Percentage error values $>' num2str(thres) '$Wm$^{-2}$']);
pcolPlot(X,Y,diff');
hold on;
contour(X,Y,diff',[thres thres],'-w','LineWidth',2);
ylim([-250 0]);
caxis([0 100]);
colorbar;
title('Difference between $\alpha>0$ and $\alpha<0$ runs');
xlabel('Time (days)');
ylabel('Depth (m)');
subplot(2,2,4);
ylim([-250 0]);

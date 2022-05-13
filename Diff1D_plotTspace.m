%%% Calculate temperature binned heat fluxes and plot:

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

% Calculate J non-local:
Jnl = zeros(Nz+1,Nt);
Jnl = -gamt*Cp*rho0;
dJnldz = (Jnl(2:end,:)-Jnl(1:(end-1),:)); % Wm-2 within each grid cell

% Define temperature grid:
dTi = 0.5;
Tie = 10:dTi:30;
Ti = avg(Tie);
NTi = length(Ti);

% Do temperature binning:
JT = zeros(NTi,Nt);
JnlT = zeros(NTi,Nt);
for ii = NTi:-1:1
    ii    
    for zi = 1:Nz
        inds = T(zi,:)>Tie(ii) & T(zi,:)<=Tie(ii+1);
        JT(ii,:) = JT(ii,:)+dJdz(zi,:).*(inds*1);
        JnlT(ii,:) = JnlT(ii,:)+dJnldz(zi,:).*(inds*1);
    end
end
JT = cat(1,zeros(1,Nt),cumsum(JT));
JnlT = cat(1,zeros(1,Nt),cumsum(JnlT));

% Do temperature binning of mean:
Tmean = mean(T,2);
dJdzmean = mean(dJdz,2);
JTmean = zeros(NTi,1);
dJnldzmean = mean(dJnldz,2);
JnlTmean = zeros(NTi,1);
for ii = NTi:-1:1
    for zi = 1:Nz
        inds = Tmean(zi)>Tie(ii) & Tmean(zi)<=Tie(ii+1);
        JTmean(ii) = JTmean(ii)+dJdzmean(zi).*(inds*1);
        JnlTmean(ii) = JnlTmean(ii)+dJnldzmean(zi).*(inds*1);
    end
end
JTmean = cat(1,zeros(1,1),cumsum(JTmean));
JnlTmean = cat(1,zeros(1,1),cumsum(JnlTmean));

% Do some plotting in T and z space:
figure;
tind = 50;
subplot(1,2,1);
plot(J(:,tind),z_w,'--k');
hold on;
plot(mean(J,2),z_w,'-k','linewidth',2);

plot(Jnl(:,tind),z_w,'--r');
plot(mean(Jnl,2),z_w,'-r','linewidth',2);

subplot(1,2,2);
plot(JT(:,tind),Tie,'--k');
hold on;
plot(mean(JT,2),Tie,'-k','linewidth',2);
plot(JTmean,Tie,':k','linewidth',2);
plot(JnlT(:,tind),Tie,'--r');
plot(mean(JnlT,2),Tie,'-r','linewidth',2);
plot(JnlTmean,Tie,':r','linewidth',2);


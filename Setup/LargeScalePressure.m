
%This script calculates large-scale SOND averaged free surface and
%buoyancy gradients for use in GOTM.

load MIXnames;
fname = BMIX_avg;
lonmn = -150;
lonmx = -130;
latmn = -0.5;
latmx = 0.5;
tmin = 240;
tmax = 360;
zbox = [1 50];
lon = ncread(fname,'lon_rho');
lat = ncread(fname,'lat_rho');
[tmp latmnI] = min(abs(lat(1,:)-latmn));
[tmp latmxI] = min(abs(lat(1,:)-latmx));
lon = lon(:,round((latmnI+latmxI)/2));
lat = lat(:,round((latmnI+latmxI)/2));
[tmp lonmnI] = min(abs(lon-lonmn));
[tmp lonmxI] = min(abs(lon-lonmx));
lon = lon(lonmnI:lonmxI);
lat = lat(lonmnI:lonmxI);
[tmp lon140] = min(abs(lon+140));
time = ncread(fname,'ocean_time')/86400;
load('ZVEC.mat');
z_rho = ZRHO;
[tmp tminI] = min(abs(time-tmin));
[tmp tmaxI] = min(abs(time-tmax));

slice = {[lonmnI lonmxI],[latmnI latmxI],...
           zbox,[tminI tmaxI]};
slice2D = {[lonmnI lonmxI],[latmnI latmxI],...
           [tminI tmaxI]};

CF = pi*6371000/180; %conversion factor from lat/lon to distance.
dx = CF*(lon(end)-lon(1))*cos(pi/180*lat(1));

zeta = squeeze(nanmean(GetVar(fname,0,{'zeta'},slice2D),2));
zeta1 = mean(zeta(1,:));
zeta2 = mean(zeta(end,:));
dzetadx = (zeta2-zeta1)/dx;

T = squeeze(nanmean(GetVar(fname,0,{'temp'},slice),2));
S = squeeze(nanmean(GetVar(fname,0,{'salt'},slice),2));
T1 = mean(T(1,:,:),3);
T2 = mean(T(end,:,:),3);
S1 = mean(S(1,:,:),3);
S2 = mean(S(end,:,:),3);
dTdx = (T2-T1)/dx;
dSdx = (S2-S1)/dx;
dSdy = 0*dSdx;
dTdy = 0*dTdx;

%Calculate PGF = -1/rho0 dPdx:
alpha = 2.489e-4; %1/degC
beta  = 7.453e-4; %1/psu
g = 9.81; %ms-2
rho0 = 1025; %kg m-3

dbdx = g*alpha*dTdx-g*beta*dSdx;
PGF = zeros(size(dbdx));
for zi=(length(PGF)-1):-1:1
    PGF(zi) = PGF(zi+1)+(dbdx(zi+1)+dbdx(zi))/2*(z_rho(zi+1)- ...
                                                 z_rho(zi));
end
PGF=PGF-g*dzetadx;
figure;
plot(PGF,z_rho,'-k');
hold on;
plot(2.9e-7*cos(2*pi*z_rho/560).*exp(z_rho/200),z_rho,'-r');
ylim([-1000 0]);
ylabel('Depth (m)');
xlabel('PGF (ms$^{-2}$)');
text(1e-7,-300,'$2.9\times10^{-7}\cos(2\pi z/560)e^{z/200}$', ...
     'color','r');
title('large-scale PGF from BMIX simulations');



%From MMIX diagnostics:
load('MMIX_140W_TSDIAG.mat');
PGF = mean(U_PRSGRD,2);
Z = mean(Z,2);
[tmp minPGFind] = min(Z+300);
minPGFind = minPGFind-1;
[minPGF tmp] = min(PGF);
PGFgauss = (-log(PGF((minPGFind+1):end)-2*minPGF)).^(1/3);
Zgauss = Z((minPGFind+1):end);
figure;
plot(PGFgauss,Zgauss);
hold on;
plot(polyval(polyfit(Zgauss,PGFgauss,1),Zgauss),Zgauss,'-r');
polyfit(Zgauss,PGFgauss,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DIFF1D - 1D Diffusion equation solver                             %
%                                                                    %
% This script solves the 1D vertical diffusion equation on a ROMS    %
% stretched vertical coordinates grid for the temperature, salinity, %
% zonal and meridional velocities.                                   %
%                                                                    %
% Boundary layer mixing is used. The interior mixing scheme can      %
% be chosen as KPP or Pacanowski & Philander or Peters 88.           %
%                                                                    %
% This script relies on several additional scripts;                  %
%                                                                    %
% Diff1Dconst.m -> contains most constants                          %
%                                                                    %
% Diff1Dplot.m -> plots the standard image for movie as running.     %
%                                                                    %
% Diff1Dmix.m -> contains the code that determines the mixing        %
% scheme diffusivity.                                                %
%                                                                    %
% Diff1Dstep.m -> step forward the fields in time using a flux       %
% formulation.                                                       %
%                                                                    %
% Ryan Holmes June 2014                                              %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup:
close all;
clear all;
Diff1Dconst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME stepping 
dt = 480;
tfin = 180*86400;%30*86400; %sec
tplot = 160;%720/4; %plot every tplot slices.
meanplot = 0; %plot mean curves on left of pcolor plots.
% $$$ pplot = 0; %plot profiles as code runs
% $$$ % $$$ colplot = 1; %plot pcolor at end.
t = 0:dt:tfin;Nt = length(t); %time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE forcing 
tau_x = -0.08; %Nm-2 = kg m-1 s-2 % 0.0766 - avg. sep-dec from BMIX_avg.
tau_y = 0; %Nm-2

ssflux = 0; %psu kg-1 m-2 s-1
srflux = 275; %radiative heat flux. near averages (277.83
              %avg. sep-dec from BMIX_avg)
shflux = -180; %surface heat flux (use srflux = 0 and shflux = total
              %for all at surface). Wm-2 = J s-1 m-2 = kg s-3
              % (181.88 avg. sep-dec from BMIX_avg).

DIR = 0; %Include a diurnal cycle?

if (DIR)
    hr = mod(t/3600,24);
    srflux = srflux*4*(cos((2*(hr/24)-1)*pi)).^2;
    srflux(hr/24<0.25)=0;
    srflux(hr/24>0.75)=0;
else
    srflux = srflux*ones(size(t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BODY forces 
% This is currently setup to restore T-S to the initial profile
% with damping coefficient TS_RST s-1 and to apply a zonal body force
% proportional to the zonal velocity itself.

%TS restoring:
TS_RST = 1/(15*86400); %nudging coefficient (s-1)

%u restoring:
u_RST  = 1/(200000000*86400); %nudging coefficient (s-1)

%v restoring:
v_RST  = 1/(200000000*86400); %nudging coefficient (s-1)

%PGF:
PGFscale = 120; %depth scale of cubic gaussian.
%set equal to wind stress:
PGFamp = -tau_x/rho0/trapz(-10000:0.1:0,exp(-(-(-10000:0.1:0)/PGFscale).^3));
%PGFamp = 3.2e-7; %Value of PGF at surface
PGF_X = [num2str(PGFamp) '*exp(-(-z_rho/' num2str(PGFscale) ').^3)'];

%TIW forcing:
SYM = 0; %0 -> body force is dvdy*u
         %1 -> body force is dvdy*u_initial
         %2 -> body force is dvdy*uF, where uF is taken from the
         %file uSYM.
rTIW = 0; %0 -> dvdy is ideal (see next few lines)
          %1 -> dvdy taken from file dvdyFILE and interpolated onto
          %grid, repeating as many times as neccessary.
dvdyFILE = 'ROMS_dvdy.mat';
uSYM = 'uSYM_PP.mat';%'uSYM_PP.mat';
period = 15*86400; %peroid of oscillation (s)
amplitude = 2.8e-6; %amplitude of stretching dv/dy.
dvdy = amplitude*sin(2*pi/period*t); %dv/dy time change
dvdy_v = '(5.2e-9/2.8e-6)*z_rho+1'; %Vertical form of dvdy

%Kelvin wave forcing:
Kperiod = 60*86400; %period of oscillaiton (s)
Kamp = 10; %amplitude (ms-1)
Kpha = pi/2; %initial phase
Kmod = 1; %Vertical mode number.
load('TSvadv.mat');
Vadv_only = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICAL mixing

%Interior:
%0 = no interior, 1 = KPP, 2 = PP, 3 = PP88.
INT = 2;

%Background:
kv0 = 1e-4; %m2s-1 interior background
kt0 = 1e-5; %m2s-1 interior background
ks0 = 1e-5; %m2s-1 interior background

% KPP 
Ri0 = 0.7; %Critical Richardson number
% $$$ Ri0 = 1; %Critical Richardson number
K0 = 2e-3; %Interior diffusivity maximum
% $$$ K0 = 4e-3; %Interior diffusivity maximum

%PP parameters:
Ri0_PP = 0.2; %Decay scale
K0_PP = 0.01; %max int. diff.
% $$$ K0_PP = 0.005; %max int. diff.
PR = 0;       %PR = 1; PP1 parameterization with Pr=1.
              %PR = 0; normal PP parametrization
              %PR = 2; PP1.2 parametrization, where kv is set to
              %kt.

%P88 parameters:
%contained in Diff1Dconst.m
P88_Kmax = 0.1; %maximum diffusivity for P88 scheme.

%Boundary layer:
KPPBL = 1; %use boundary layer mixing.
EKMO = 1; %use Ekman/MO length restrictions.
nl = 1; %use non-local fluxes?
KPPMLD = -15; %Initial guess mixed layer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL conditions

% Initial profiles from mean BMIX profiles:
load('BMIX_140W_TS.mat');
zI = mean(Z,2);
uI = mean(U,2);
vI = mean(V,2)*0; %SET TO ZERO!!!
TI = mean(T,2);
SI = mean(S,2);
bI = g*alpha*TI-g*beta*SI;
zwI = avg(zI);

%Kelvin wave B modes:
P = sw_pres(-zI,0);
KN2 = sw_bfrq(SI,sw_temp(SI,TI,P,0),P,lat);
Kz = [-4000; zwI; 0;];
KN2 = [KN2(1); KN2; KN2(end)];
[Kc,KS,KC] = Bmodes(Kz,KN2,Kmod);
Kc = Kc(Kmod);KS = KS(:,Kmod);KC=KC(:,Kmod);
KS = interp1(Kz,KS,z_rho,'spline');
KC = interp1(avg(Kz),KC,z_rho,'spline');
KN2 = interp1(Kz,KN2,z_rho,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKINGCODE 

% SETUP 
%Distributed forcing fluxes:
TAU_X = [zeros(Nz,1); tau_x];
TAU_Y = [zeros(Nz,1); tau_y];
TAU_S = [zeros(Nz,1); ssflux];

%Initial body-forces:
BF_X = zeros(Nz,1);
BF_Y = zeros(Nz,1);
BF_T = zeros(Nz,1);
BF_S = zeros(Nz,1);

%TIW forcing:
if (rTIW==0) %idealized dv/dy
eval(['dvdy = repmat(dvdy,[Nz 1]).*repmat(' dvdy_v ...
      ',[1 Nt]);']);
elseif (rTIW==1) %read dvdy from file dvdyFILE:
    load(dvdyFILE);
    ROMSNt = length(ROMSdvdy(1,:));
    dvdyDL = zeros(Nz,ROMSNt);
    for tiii=1:ROMSNt
        dvdyDL(:,tiii) = interp1(ROMSz,ROMSdvdy(:,tiii),z_rho,'spline');
    end
    dvdyDL = repmat(dvdyDL,[1 ceil(t(end)/86400/ROMSNt)]);
    dvdy = zeros(Nz,Nt);
    for ziii=1:Nz
        dvdy(ziii,:) = interp1((0:1:(length(dvdyDL(1,:))-1))*86400,dvdyDL(ziii,:),...
                               t,'linear');
    end
end
    
%Wind-stress balancing PGF:
eval(['PGF_X = ' num2str(PGF_X) ';']);

%Kelvin wave forcing:
Komega = 2*pi/Kperiod;
PGF_K = -Kamp*Komega*real(i*exp(-i*Komega*repmat(t,[Nz 1])+i*Kpha)).*repmat(KC,[1 ...
                    Nt]);
KW = -Kamp*Komega*Kc/g*real(i*exp(-i*Komega*repmat(t,[Nz 1])+i* ...
                                  Kpha)).*repmat(KS,[1 Nt]);

%Intialize arrays:
u = zeros(Nz,Nt);v = zeros(Nz,Nt);T = zeros(Nz,Nt);S = zeros(Nz,Nt);b = zeros(Nz,Nt);%z_rho variables
kv = zeros(Nz+1,Nt);kt = zeros(Nz+1,Nt);ks = zeros(Nz+1,Nt);%z_w variables
Hsbl = zeros(Nt,1);gamv=zeros(Nz+1,Nt);gamt=zeros(Nz+1,Nt);gams= ...
       zeros(Nz+1,Nt);

%SYM = 2 symmetric stretching base velocity profile:
if (SYM == 2)
    load(uSYM);
end

%Initialize mom. balance arrays:
UDIV=zeros(Nz,Nt);URST=zeros(Nz,Nt);UTND=zeros(Nz,Nt);UMFX=zeros(Nz,Nt);
UPGF=zeros(Nz,Nt);
%Temp. eqn:
TRST=zeros(Nz,Nt);THFX=zeros(Nz,Nt);TTND=zeros(Nz,Nt);
%SST budget:
% $$$ SST=zeros(Nt,1);SSTTND = zeros(Nt,1);SSTENT = zeros(Nt,1);SSTSHX=zeros(Nt,1);
% $$$ SSTTUB=zeros(Nt,1);SSTRST=zeros(Nt,1);

%Interp from initial profiles:
u(:,1) = interp1(zI,uI,z_rho,'spline');v(:,1) = interp1(zI,vI,z_rho,'spline');
T(:,1) = interp1(zI,TI,z_rho,'spline');S(:,1) = interp1(zI,SI,z_rho,'spline');
kv(:,:) = kv0;kt(:,:) = kt0;ks(:,:) = ks0;b(:,1) = g*alpha*T(:,1)- ...
          g*beta*S(:,1);Hsbl(1) = KPPMLD;

%Calculate constant mixing parameters:
Ustar = sqrt(sqrt(tau_x.^2+tau_y.^2)/rho0); %Friction velocity (const).
hekman = 0.7*Ustar/max(abs([f 1e-10])); %Ekman depth (const).
wt0 = shflux/Cp/rho0;
ws0 = ssflux/rho0;

% TIME stepping start
for ti = 1:(length(t)-1)
    ['Doing step ' num2str(ti) ' of ' num2str(length(t)-1)]
    %Calculate diffusivities:
    Diff1Dmix; 
    
    %Calculate heat flux down:
    TAU_T = [zeros(Nz,1); shflux/Cp];
    for zi = 1:length(z_w)
        TAU_T(zi) = TAU_T(zi)+srflux(ti)/Cp*swdk(z_w(zi)); %degC kg m-2 s-1 = degC s-1 m * rho0
    end
    
    %SST budget:
% $$$     SST(ti)  = (sum(T(z_rho>=KPPMLD,ti).*Hz(z_rho>=KPPMLD))+ ...
% $$$         (min(z_rho(z_rho>=KPPMLD))-KPPMLD)*(T(min(find(z_rho>= ...
% $$$                                                       KPPMLD)),ti)))/(-KPPMLD);
% $$$     SSTSHX(ti) = (shflux+(1-swdk(KPPMLD))*srflux(ti))/(-KPPMLD*rho0*Cp);
% $$$     SSTRST(ti) = -TS_RST*(SST(ti)-SST(1));
% $$$     SSTENT(ti) = -1/(-KPPMLD)*(KPPMLD-KPPMLDp)/dt*(SST(ti)-interp1(z_rho,T(:,ti),KPPMLD,'spline'));
% $$$     SSTTUB(ti) = -interp1(z_w(2:(end-1)),kt(2:(end-1),ti).*((T(2:end,ti)-...
% $$$                            T(1:(end-1),ti))./Hzw+gamt(2:(end-1),ti)), ...
% $$$                           KPPMLD,'spline')/(-KPPMLD);
    
    %Calculate body forces:
    if (SYM==1)
        UDIV(:,ti) = dvdy(:,ti).*u(:,1);
    elseif (SYM==2)
        UDIV(:,ti) = dvdy(:,ti).*uF(:,1);
    elseif (SYM == 0)
        UDIV(:,ti) = dvdy(:,ti).*u(:,ti);
    end
    URST(:,ti) = -u_RST*(u(:,ti)-u(:,1));
    UPGF(:,ti) = PGF_X+PGF_K(:,ti);
    TRST(:,ti) = -TS_RST*(T(:,ti)-Tvad(:,ti));
    TVAD(:,ti) = -interp1(avg(z_rho),diff(T(:,ti))./diff(z_rho), ...
                          z_rho,'spline').*KW(:,ti);
    SVAD(:,ti) = -interp1(avg(z_rho),diff(S(:,ti))./diff(z_rho), ...
                          z_rho,'spline').*KW(:,ti);
    
    if (~Vadv_only)
    BF_X = URST(:,ti)+UDIV(:,ti)+UPGF(:,ti);
    BF_T = TRST(:,ti)+TVAD(:,ti);
    BF_Y = -v_RST*(v(:,ti)-v(:,1));
    BF_S = -TS_RST*(S(:,ti)-Svad(:,ti))+SVAD(:,ti);
    else %Apply only vertical advection
    BF_X = 0*URST(:,ti);
    BF_T = TVAD(:,ti);
    BF_Y = 0*-v_RST*(v(:,ti)-v(:,1));
    BF_S = SVAD(:,ti);
    kv(:,ti) = 0;kt(:,ti) = 0;ks(:,ti) = 0;
    gamv(:,ti) = 0;gamt(:,ti) = 0;gams(:,ti) = 0;
    TAU_X(:) = 0;TAU_T(:) = 0;TAU_S(:) = 0;
    end

    %Calculate step ti+1:
    [u(:,ti+1),UMFX(:,ti)] = Diff1Dstep(u(:,ti),kv(:,ti),gamv(:,ti),Hz,Hzw,-TAU_X/rho0,BF_X,Nz,dt);
    [v(:,ti+1),tmp] = Diff1Dstep(v(:,ti),kv(:,ti),gamv(:,ti),Hz,Hzw,-TAU_Y/rho0,BF_Y,Nz,dt);
    [T(:,ti+1),THFX(:,ti)] = Diff1Dstep(T(:,ti),kt(:,ti),gamt(:,ti),Hz,Hzw,-TAU_T/rho0,BF_T,Nz,dt);
    [S(:,ti+1),tmp] = Diff1Dstep(S(:,ti),ks(:,ti),gams(:,ti),Hz,Hzw,-TAU_S/rho0,BF_S,Nz,dt);
    b(:,ti+1) = g*alpha*T(:,ti+1)-g*beta*S(:,ti+1);
    UTND(:,ti) = (u(:,ti+1)-u(:,ti))/dt;
    TTND(:,ti) = (T(:,ti+1)-T(:,ti))/dt;
end

%save('TIW_PP.mat');
%save('Kelvin_vadvNUD.mat');
% $$$ save('RestoringTS_Sensitivity/PP_nDIR_nTIW_rst15.mat');

% $$$ Diff1Dcolplot;
% $$$ ti = 1:(length(t)-1);
% $$$ Diff1Dsetplot;
% $$$ Diff1Dplot;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DIFF1D - 1D Diffusion equation solver                             %
%                                                                    %
% This script solves the 1D vertical diffusion equation on a ROMS    %
% stretched vertical coordinates grid for the temperature, salinity, %
% zonal and meridional velocities.                                   %
%                                                                    %
% Boundary layer mixing is used. The interior mixing scheme can      %
% be chosen as KPP or Pacanowski & Philander or Peters 88.           %
%                                                                    %
% To do:                                                             %
%       - Add in double diffusive contribution to interior           %
%       KPP?                                                         %
%       - Add in a PGF to counter the wind-stress?                   %
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
clear all;Diff1Dconst;

%%%% RUN NUMBER %%%%%%%
run = 12;
out_folder = './data/'; %folder to save data to (will be labeled with
                      %run number). 

%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME stepping 
dt = 120;
tfin = 45*86400; %45*86400; %sec
t = 0:dt:tfin;Nt = length(t); %time
NOUT = 1; %output averaged every NOUT time steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE forcing 
period = 15*86400; %peroid of oscillation (s) 15*86400
amplitude = 0.04; %amplitude of wind stress oscillation (Nm-2)
tau_x = -0.08*ones(size(t)) ;%+ amplitude*sin(2*pi/period*t); %Nm-2 = kg m-1 s-2
tau_y = 0*ones(size(t)); %Nm-2

ssflux = 0*ones(size(t)); %psu kg-1 m-2 s-1
period = 15*86400; %peroid of oscillation (s)
amplitude = 150; %2.8e-6; %amplitude of heat flux variation
shflux = -180*ones(size(t)); %surface heat flux (use srflux = 0 and shflux = total
              %for all at surface). Wm-2 = J s-1 m-2 = kg s-3

DIR = 0; %Include a diurnal cycle?
srflux = 275; %275; %radiative heat flux.
if (DIR)
    hr = mod(t/3600,24);
    srflux = srflux*4*(cos((2*(hr/24)-1)*pi)).^2;
    srflux(hr/24<0.25)=0;
    srflux(hr/24>0.75)=0;
else
    srflux = srflux*ones(size(t)) + amplitude*sin(2*pi/period*t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BODY forces 
% This is currently setup to restore T-S to the initial profile
% with damping coefficient TS_RST s-1 and to apply a zonal body force
% proportional to the zonal velocity itself.

%Depth independent nudging:
TS_RST = 1/(15*86400); %nudging coefficient (s-1)
u_RST  = 1/(200000000*86400); %nudging coefficient (s-1)
v_RST  = 1/(200000000*86400); %nudging coefficient (s-1)

%Pressure-gradient force for EUC:
PGFscale = 120; %depth scale of cubic gaussian.
%set equal to wind stress:
PGFamp = -tau_x(1)/rho0/trapz(-10000:0.1:0,exp(-(-(-10000:0.1:0)/PGFscale).^3));
%PGFamp = 3.2e-7; %Value of PGF at surface
PGF_X = PGFamp.*exp(-(-z_rho/PGFscale).^3);

%Vertical advection:
w = zeros(Nz+1,Nt);

%TIW forcing:
SYM = 0; %0 -> body force is dvdy*u
         %1 -> body force is dvdy*u_initial
period = 15*86400; %peroid of oscillation (s)
amplitude = 0; %2.8e-6; %amplitude of stretching dv/dy.
dvdy = amplitude*sin(2*pi/period*t); %dv/dy time change
dvdy_v = '(5.2e-9/2.8e-6)*z_rho+1'; %Vertical form of dvdy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICAL mixing

%Interior:
%0 = no interior, 1 = KPP, 2 = PP, 3 = PP88.
INT = 1;

%Background:
kv0 = 1e-4; %m2s-1 interior background
kt0 = 1e-5; %m2s-1 interior background
ks0 = 1e-5; %m2s-1 interior background

if (INT==1)     % KPP
    Ri0 = 0.7; %Critical Richardson number
    K0 = 2e-3; %Interior diffusivity maximum

elseif (INT == 2) %PP
    Ri0_PP = 0.2; %Decay scale
    K0_PP = 0.01; %max int. diff.
    PR = 0;      %PR = 1; PP1 parameterization with Pr=1.
                 %PR = 0; normal PP parametrization
                 %PR = 2; PP1.2 parametrization, where kv is set to
                 %kt.

elseif (INT == 3)
    %P88 parameters: %contained in Diff1Dconst.m
    P88_Kmax = 0.1;
end

%Boundary layer:
KPPBL = 1; %use boundary layer mixing.
EKMO = 1; %use Ekman/MO length restrictions.
nl = 1; %use non-local fluxes?
KPPMLD = -15; %Initial guess mixed layer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL conditions

% Initial profiles from mean BMIX profiles:
load('data/BMIX_140W_TS.mat');
zI = mean(Z,2);
uI = mean(U,2);
vI = mean(V,2)*0; %SET TO ZERO!!!
TI = mean(T,2);
SI = mean(S,2);
bI = g*alpha*TI-g*beta*SI;
zwI = (zI(2:end)+zI(1:(end-1)))/2;

% $$$ % Initial profiles from TPOS20 Sep-Dec:
% $$$ load('TPOS20_140W_TS.mat');
% $$$ dvec = datevec(time);
% $$$ inds = dvec(:,2)>=9;
% $$$ zI = Z;
% $$$ uI = mean(U(:,inds),2);
% $$$ vI = 0*mean(V(:,inds),2);
% $$$ TI = mean(T(:,inds),2);
% $$$ SI = mean(S(:,inds),2);
% $$$ bI = g*alpha*TI-g*beta*SI;
% $$$ zwI = (zI(2:end)+zI(1:(end-1)))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKINGCODE 

% SETUP 

%Initial body-forces:
BF_X = zeros(Nz,1);
BF_Y = zeros(Nz,1);
BF_T = zeros(Nz,1);
BF_S = zeros(Nz,1);

eval(['dvdy = repmat(dvdy,[Nz 1]).*repmat(' dvdy_v ...
      ',[1 Nt]);']);

%Intialize arrays:
u = zeros(Nz,Nt);v = zeros(Nz,Nt);T = zeros(Nz,Nt);
S = zeros(Nz,Nt);b = zeros(Nz,Nt);%z_rho variables
kv = zeros(Nz+1,Nt);kt = zeros(Nz+1,Nt);ks = zeros(Nz+1,Nt);
kt_int = zeros(Nz+1,Nt);kt_bl = zeros(Nz+1,Nt);
gamv=zeros(Nz+1,Nt);gamt=zeros(Nz+1,Nt);gams= zeros(Nz+1,Nt);
bulkRiN = zeros(Nz+1,Nt);bulkRiD = zeros(Nz+1,Nt);%z_w variables
Hsbl = zeros(1,Nt);

%Interp from initial profiles:
u(:,1) = interp1(zI,uI,z_rho,'spline');v(:,1) = interp1(zI,vI,z_rho,'spline');
T(:,1) = interp1(zI,TI,z_rho,'spline');S(:,1) = interp1(zI,SI,z_rho,'spline');
kv(:,:) = kv0;kt(:,:) = kt0;ks(:,:) = ks0;
b(:,1) = g*alpha*T(:,1)-g*beta*S(:,1);Hsbl(1) = KPPMLD;

% TIME stepping start
for ti = 1:(length(t)-1)
    if (mod(ti,50)==0)
    ['Doing step ' num2str(ti) ' of ' num2str(length(t)-1)]
    end
    
    %Calculate mixing parameters from time dependent forcing:
    Ustar = sqrt(sqrt(tau_x(ti).^2+tau_y(ti).^2)/rho0); %Friction velocity (const).
    hekman = 0.7*Ustar/max(abs([f 1e-10])); %Ekman depth
    wt0 = shflux(ti)/Cp/rho0;
    ws0 = ssflux(ti)/rho0;
            
    %Calculate diffusivities:
    Diff1Dmix; 
    bulkRiN(:,ti) = RiKPP_Numer;
    bulkRiD(:,ti) = RiKPP_Denom;
    
    %Calculate heat flux down:
    TAU_T = [zeros(Nz,1); shflux(ti)/Cp];
    for zi = 1:length(z_w)
        TAU_T(zi) = TAU_T(zi)+srflux(ti)/Cp*swdk(z_w(zi)); %degC kg m-2 s-1 = degC s-1 m * rho0
    end
    
    %Calculate momentum and salinity fluxes:
    TAU_X = [zeros(Nz,1); tau_x(ti)];
    TAU_Y = [zeros(Nz,1); tau_y(ti)];
    TAU_S = [zeros(Nz,1); ssflux(ti)];

    %Calculate body forces:
    %depth-independent restoring:
    URST = -u_RST*(u(:,ti)-u(:,1));
    VRST = -v_RST*(v(:,ti)-v(:,1));
    TRST = -TS_RST*(T(:,ti)-T(:,1));
    SRST = -TS_RST*(S(:,ti)-S(:,1));

    %Vertical advection:
    UVAD = zeros(Nz+1,1);
    UVAD(2:(end-1),:) = -diff(u(:,ti))./diff(z_rho).*w(2:(end-1),ti);
    UVAD = avg(UVAD);
    VVAD = zeros(Nz+1,1);
    VVAD(2:(end-1),:) = -diff(v(:,ti))./diff(z_rho).*w(2:(end-1),ti);
    VVAD = avg(VVAD);
    TVAD = zeros(Nz+1,1);
    TVAD(2:(end-1),:) = -diff(T(:,ti))./diff(z_rho).*w(2:(end-1),ti);
    TVAD = avg(TVAD);
    
    %Zonal pressure gradient:
    UPGF = PGF_X;
    
    %TIW Stretching:
    UDIV = dvdy(:,ti).*u(:,ti);
  
    BF_X = URST+UPGF+UDIV;
    BF_Y = VRST;
    BF_T = TRST;
    BF_S = SRST;
    
    %Calculate step ti+1:
    [u(:,ti+1),tmp] = Diff1Dstep(u(:,ti),kv(:,ti),gamv(:,ti),Hz,Hzw,-TAU_X/rho0,BF_X,Nz,dt);
    [v(:,ti+1),tmp] = Diff1Dstep(v(:,ti),kv(:,ti),gamv(:,ti),Hz,Hzw,-TAU_Y/rho0,BF_Y,Nz,dt);
    [T(:,ti+1),tmp] = Diff1Dstep(T(:,ti),kt(:,ti),gamt(:,ti),Hz,Hzw,-TAU_T/rho0,BF_T,Nz,dt);
    [S(:,ti+1),tmp] = Diff1Dstep(S(:,ti),ks(:,ti),gams(:,ti),Hz,Hzw,-TAU_S/rho0,BF_S,Nz,dt);
    b(:,ti+1) = g*alpha*T(:,ti+1)-g*beta*S(:,ti+1);

end

Diff1Dredout;
if (exist(out_folder)==7 | exist(out_folder)==5)
    ['Saving run ' num2str(run) ' to folder ' out_folder]
    save(sprintf([out_folder 'run_%03d.mat'],run));
end


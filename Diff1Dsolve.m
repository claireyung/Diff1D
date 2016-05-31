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
%       - add in nonlocal transport part -> only nonzero in          %
%       unstable situations! -> Almost always stable as long as ML   %
%       is at least as deep as first grid point (to absorb some SW   %
%       radiation inside the mixed layer depth), which is set as a   %
%       minimum depth.                                               %
%       - Add in double diffusive contribution to interior           %
%       KPP?                                                         %
%       - Add in diurnal cycle? -> Needs nonlocal fluxes             %
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
clear all;
Diff1Dconst;

%%%% RUN NUMBER %%%%%%%
run = 15;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME stepping 
dt = 30;
tfin = 16*86400; %sec
t = 0:dt:tfin;Nt = length(t); %time
NOUT = 1; %output averaged every NOUT time steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE forcing 
AL = 0.12;
AH = 0;
omgL = 0.09*1e-4;
omgH = 0.5*1e-4;
tau_x = AL*cos(omgL*t)+AH*cos(omgH*t); %Nm-2 = kg m-1 s-2
tau_y = 0*ones(size(t)); %Nm-2

ssflux = 0*ones(size(t)); %psu kg-1 m-2 s-1
shflux = 0*ones(size(t)); %surface heat flux (use srflux = 0 and shflux = total
              %for all at surface). Wm-2 = J s-1 m-2 = kg s-3

DIR = 0; %Include a diurnal cycle?
srflux = 0; %radiative heat flux.
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

%Depth independent nudging:
TS_RST = 1/(15*86400); %nudging coefficient (s-1)
u_RST  = 1/(2*86400); %nudging coefficient (s-1)
v_RST  = 1/(2*86400); %nudging coefficient (s-1)

%Coriolis terms:
% $$$ load('testvars_rmholmes.mat','dudytest1');
dudytest1 = 5e-5;
ucor = (f-dudytest1)*ones(Nz,1); %Coriolis coefficient in u-mom eqn.
vcor = (-f*ones(Nz,1)); %Coriolis coefficient in v-mom eqn.

%Advection from ROMS output:
dtROMS = 10800;
tROMS = 0:dtROMS:(599*dtROMS);

%Vertical advecting velocity:
% $$$ load('rmholmes_w_rho_buoy65mdia.mat','w_zt','z_wROMS');
w = zeros(Nz+1,Nt); %Vertical advection velocity
% $$$ wROMS = filter_field(w_zt,21,'-t');
% $$$ wROMS(isnan(wROMS)) = 0;
% $$$ for ii = 1:length(z_w)
% $$$     [tmp ind] = min(abs(z_wROMS-z_w(ii)));
% $$$     w(ii,:) = interp1(tROMS',wROMS(ind,:)',t,'linear');
% $$$ end
% $$$ w = 0*w;

%Temperature advection:
load('rmholmes_temp_adv_220.mat');
temp_adv = zeros(length(z_rhoROMS),Nt);
for ii=1:length(z_rhoROMS)
    temp_adv(ii,:) = interp1(tROMS',temp_hadv220(ii,:)+...
        temp_vadv220(ii,:),t,'linear');
end
Tadv = zeros(Nz,Nt);
for ii=1:length(t)
    Tadv(:,ii) = interp1(z_rhoROMS,temp_adv(:,ii),...
        z_rho,'linear');
end
twin = 9*(dtROMS/dt);
Tadv = filter_field(Tadv,twin + (1-mod(twin,2)),'-t');
Tadv(isnan(Tadv)) = 0;
Tadv = filter_field(Tadv',5,'-t')';
Tadv(isnan(Tadv)) = 0;
clear temp_adv;

% %Horizontal advection terms:
% %dTdy:
%       % Constant dTdy:
%      load('testvars_rmholmes.mat','dTdytest1');
%      tyad = (-dTdytest1*ones(Nz,Nt));
% %      % Ideal Ekman varying dTdy:
% %      dTdy = diff(temp(:,end))./diff(y_rho(1,:)');
% %      y0 = y_rho(1,214);
% %      ypos = y0-cumsum(-tau_x./0.7.*sqrt(rho0./abs(tau_x))*...
% %                       dt*(f/(f-dudytest1))/rho0);
% %     tyad = -repmat(interp1(avg(y_rho(1,:)'),dTdy,ypos,'linear'),[Nz 1]);
% %      % v-advected dTdy:
% %      dTdy = diff(temp(:,end))./diff(y_rho(1,:)'); %dTdy vector
% %      y0 = y_rho(1,214); %initial position
% %      yp = zeros(Nz,Nt);yp(:,1) = y0;
% %      Y_RST = 1/(8*86400); %y-restoring coefficient.
% %      tyad = zeros(Nz,Nt);
% %      
% %PV<0 depth-dependent restoring:
%       % Constant dbdy:
%       dbdy = g*alpha*dTdytest1*ones(1,Nt);
% %      % Ideal Ekman varying dTdy:
% % % $$$      dbdy = g*alpha*(-tyad(1,:)); 
%       TS_PVth = dbdy.^2./(f-dudytest1)./f; %Threshold N2 for PV restoring
%       TS_PVRST = (f-dudytest1); %Coefficient for PV restoring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICAL mixing

%Interior:
%0 = no interior, 1 = KPP, 2 = PP, 3 = PP88.
INT = 1;

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
%Ri0_PP = 0.2; %Decay scale
%K0_PP = 0.01; %max int. diff.
% $$$ K0_PP = 0.005; %max int. diff.
%PR = 0;       %PR = 1; PP1 parameterization with Pr=1.
              %PR = 0; normal PP parametrization
              %PR = 2; PP1.2 parametrization, where kv is set to
              %kt.

%P88 parameters:
%contained in Diff1Dconst.m
%P88_Kmax = 0.1; %maximum diffusivity for P88 scheme.

%Boundary layer:
KPPBL = 1; %use boundary layer mixing.
EKMO = 1; %use Ekman/MO length restrictions.
nl = 1; %use non-local fluxes?
KPPMLD = -20; %Initial guess mixed layer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL conditions

% Initial temperature profile from Dan's Sims:
load('testvars_rmholmes.mat','t0test1','z_in');
zI = z_in';
uI = 0*zI;
vI = 0*zI;
TI = t0test1';
SI = 35*ones(size(zI));
bI = g*alpha*TI-g*beta*SI;
zwI = (zI(2:end)+zI(1:(end-1)))/2;

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

%Intialize arrays:
u = zeros(Nz,Nt);v = zeros(Nz,Nt);T = zeros(Nz,Nt);
S = zeros(Nz,Nt);b = zeros(Nz,Nt);%z_rho variables
kv = zeros(Nz+1,Nt);kt = zeros(Nz+1,Nt);ks = zeros(Nz+1,Nt);
gamv=zeros(Nz+1,Nt);gamt=zeros(Nz+1,Nt);gams= zeros(Nz+1,Nt);
bulkRiN = zeros(Nz+1,Nt);bulkRiD = zeros(Nz+1,Nt);%z_w variables
Hsbl = zeros(1,Nt);

%Interp from initial profiles:
u(:,1) = interp1(zI,uI,z_rho,'spline');v(:,1) = interp1(zI,vI,z_rho,'spline');
T(:,1) = interp1(zI,TI,z_rho,'spline');S(:,1) = interp1(zI,SI,z_rho,'spline');
kv(:,:) = kv0;kt(:,:) = kt0;ks(:,:) = ks0;
b(:,1) = g*alpha*T(:,1)-g*beta*S(:,1);Hsbl(1) = KPPMLD;

S(:,:) = 35;

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
    
% $$$     %Calculate heat flux down:
% $$$     TAU_T = [zeros(Nz,1); shflux(ti)/Cp];
% $$$     for zi = 1:length(z_w)
% $$$         TAU_T(zi) = TAU_T(zi)+srflux(ti)/Cp*swdk(z_w(zi)); %degC kg m-2 s-1 = degC s-1 m * rho0
% $$$     end

    %Distributed forcing fluxes:
    TAU_X = [zeros(Nz,1); tau_x(ti)];
    TAU_Y = [zeros(Nz,1); tau_y(ti)];
% $$$     TAU_S = [zeros(Nz,1); ssflux(ti)];
    TAU_T = zeros(Nz+1,1);

    %T Restoring restoring:
    TS_RSTwork = TS_RST*ones(Nz,1);
%    % PV Restoring:
%     ind = find(N2<TS_PVth(ti),1,'first');
%     TS_RSTwork(ind:end) = TS_RSTwork(ind:end)+TS_PVRST;

    %Calculate body forces:
    %depth-indepent restoring:
    URST = -u_RST*(u(:,ti)-u(:,1));
    VRST = -v_RST*(v(:,ti)-v(:,1));
    TRST = -TS_RSTwork.*(T(:,ti)-T(:,1));
    %Coriolis:
    UCOR = ucor.*v(:,ti);
    VCOR = vcor.*u(:,ti);

    %Horizontal advection:
%    tyad(:,ti) = -interp1(avg(y_rho(1,:)'),dTdy,yp(:,ti),'linear');
%    TYAD = tyad(:,ti).*v(:,ti);
    TYAD = Tadv(:,ti);
    %Only non-zero above Hsbl:
    % TYAD(z_rho>Hsbl(ti)) = 0;
    
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
  
    BF_X = UCOR+URST+UVAD;
    BF_Y = VCOR+VRST+VVAD;
    BF_T = TYAD+TRST+TVAD;
% $$$     BF_S = -TS_RST*(S(:,ti)-S(:,1));
    
    %Calculate step ti+1:
    [u(:,ti+1),tmp] = Diff1Dstep(u(:,ti),kv(:,ti), ...
                                 zeros(Nz+1,1),Hz,Hzw,-TAU_X/rho0,BF_X,Nz,dt);
    [v(:,ti+1),tmp] = Diff1Dstep(v(:,ti),kv(:,ti), ...
                                 zeros(Nz+1,1),Hz,Hzw,-TAU_Y/rho0,BF_Y,Nz,dt);
    [T(:,ti+1),tmp] = Diff1Dstep(T(:,ti),kt(:,ti), ...
                                 zeros(Nz+1,1),Hz,Hzw,-TAU_T/rho0,BF_T,Nz,dt);
% $$$     [S(:,ti+1),tmp] = Diff1Dstep(S(:,ti),ks(:,ti),gams(:,ti),Hz,Hzw,-TAU_S/rho0,BF_S,Nz,dt);
    b(:,ti+1) = g*alpha*T(:,ti+1)-g*beta*S(:,ti+1);
% $$$     UTND(:,ti) = (u(:,ti+1)-u(:,ti))/dt;
% $$$     TTND(:,ti) = (T(:,ti+1)-T(:,ti))/dt;
    %y-advected dTdy:
  %  yp(:,ti+1) = yp(:,ti)-v(:,ti)*dt-Y_RST*(yp(:,ti)-y0);
end

Diff1Dredout;
save(sprintf('data/run_%03d.mat',run));


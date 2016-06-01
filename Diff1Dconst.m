%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DIFF1Dconst               %
% This script contains constants used in Diff1D %
%                                               %
% Ryan Holmes June 2014                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICAL coordinates ROMS
N = 50;
h = 4000;
Vtransform = 1;
Vstretching = 1;
theta_s = 5;
theta_b = 0;
hc = 75;
zbot = -300;

%Generate coordinates
[z_rho,s_rho,Cs_rho] = scoord(h,0,0,Vtransform,Vstretching,theta_s,...
                              theta_b,hc,N,0,1,1,0,1);
z_rho = z_rho(z_rho>zbot)';
[z_w,s_w,Cs_w] = scoord(h,0,0,Vtransform,Vstretching,theta_s,...
                              theta_b,hc,N,1,1,1,0,1);
Nz = length(z_rho);
z_w = z_w((N+1-Nz):(N+1))';

%High res regular coodinates:
% $$$ z_rho = (-299:2:-1)';
% $$$ z_w = (-300:2:0)';

Hz = z_w(2:end)-z_w(1:(end-1));
Hzw = z_rho(2:end)-z_rho(1:(end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL constants
rho0 = 1025; %kg m-3
alpha = 2.489e-4; %1/degC
beta  = 7.453e-4; %1/psu
Cp = 4000; %J/kg/degC = m2 s-2 1/degC
g = 9.81; %ms-2
lat = 0; %Latitude
f = 2*7.29e-5*sin(lat*pi/180); %Coriolis parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHORTWAVE absorption curve:
lmd_mu1 = 0.35;
lmd_mu2 = 23;
lmd_r1 = 0.58;
cff1=1/lmd_mu1
cff2=1/lmd_mu2
swdk = @(z)(exp(z(1)*cff1)*lmd_r1+...
       exp(z(1)*cff2)*(1-lmd_r1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY layer constants KPP
Ric = 0.3; %critical Ri for MLD.
vonKar = 0.41;  %VonKarmen constant
epsl = 0.1; %Surface layer fraction (epsilon)
Cstar = 10; %C star parameter
small = 1e-10; %small parameter

%Dimensionless flux profile constants (appendix B Large 1994):
zetam = -0.2; %Maximum zeta value for mom
zetas = -1; %Maximum zeta value for scalars
am = 1.257;
cm = 8.36;
as = -28.86;
cs = 98.96;

%Turbulent velicity scale constant:
Vtc = 1.25*sqrt(0.2)/(sqrt(cs*epsl)*Ric*vonKar*vonKar);

%Diffusivity shape function coefficients:
% $$$ a0 = 0.001; a1 = 1; 
a0 = 0; a1 = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Peters 88 constants:
% K0_P88_L_m = 5.6e-8; %low Ri K0 for mom.
% K0_P88_L_s = 3.0e-9; %low Ri K0 for scl.
% EX_P88_L_m = -8.2; %low Ri exponent for mom.
% EX_P88_L_s = -9.6; %low Ri exponent for scl.
% 
% K0_P88_U_m = 5e-4; %high Ri K0 for mom.
% K0_P88_U_s = 5e-4; %high Ri K0 for scls.
% EX_P88_U_m = -1.5; %high Ri exponent for mom.
% EX_P88_U_s = -2.5; %high Ri exponent for scl.
% Ri0_P88_m = 0.2; %decay scale high Ri for mom.
% Ri0_P88_s = 0.2; %decay scale high Ri for scls.

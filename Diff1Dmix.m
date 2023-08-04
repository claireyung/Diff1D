%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        DIFF1Dmix                                %
%                                                                 %
% This script calculates the parametrized diffusivity profiles at %
% time ti given the forcing and fields u,v,T,S(:,ti).             %
%                                                                 %
% Ryan Holmes June 2014                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2 = zeros(size(z_w)); %stratification
N2(2:(end-1)) = (b(2:end,ti)-b(1:(end-1),ti))./Hzw;
N2(1) = N2(2);N2(end) = N2(end-1);


if (INT == 1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KPPINTERIOR Mixing

kint = zeros(size(kv(:,ti)));

%Calculate Ri:
Sh2 = zeros(size(z_w)); %Shear squared
Sh2(2:(end-1)) = ((u(2:end,ti)-u(1:(end-1),ti))./Hzw).^2+ ...
    ((v(2:end,ti)-v(1:(end-1),ti))./Hzw).^2;
Sh2(1) = Sh2(2);Sh2(end) = Sh2(end-1);

Ri = N2./Sh2; %Richardson number

%Set diffusivity:
kint(Ri<0) = K0; 
kint(Ri>0 & Ri<Ri0) = K0*(1-(Ri(Ri>0 & Ri<Ri0)/Ri0).^2).^3;

%Add on these diffusivities:
kv(:,ti) = kv(:,ti)+kint;
kt(:,ti) = kt(:,ti)+kint;
ks(:,ti) = ks(:,ti)+kint;      

%Add to kt_int
kt_int(:,ti) = kt_int(:,ti)+kint;

elseif (INT == 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PPINTERIOR Mixing
kint = zeros(size(kv(:,ti)));

%Calculate Ri:
Sh2 = zeros(size(z_w)); %Shear squared
Sh2(2:(end-1)) = ((u(2:end,ti)-u(1:(end-1),ti))./Hzw).^2+ ...
    ((v(2:end,ti)-v(1:(end-1),ti))./Hzw).^2;
Sh2(1) = Sh2(2);Sh2(end) = Sh2(end-1);

Ri = N2./Sh2; %Richardson number

%Set diffusivity:
kint(Ri<0) = K0_PP; 
kint(Ri>=0) = K0_PP*(1+Ri(Ri>=0)/Ri0_PP).^(-2);

%Add on these diffusivities:
kv(:,ti) = kv(:,ti)+kint;
kint(Ri<0)  = K0_PP;
kint(Ri>=0) = kv(Ri>=0,ti)./(1+Ri(Ri>=0)/Ri0_PP);
kt(:,ti) = kt(:,ti)+kint;
ks(:,ti) = ks(:,ti)+kint;

if (PR==1)% For PP1 parametrization, set Prandtl = 1.
    kint = (kv(:,ti)+kt(:,ti))/2;
    kv(:,ti) = kint;
    kt(:,ti) = kint;
    ks(:,ti) = kint;
elseif (PR==2)% For PP1.2 parametrization, set kv=kt.
    kv(:,ti) = kt(:,ti);
end

elseif (INT == 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PP88INTERIOR Mixing
kint = zeros(size(kv(:,ti)));

%Calculate Ri:
Sh2 = zeros(size(z_w)); %Shear squared
Sh2(2:(end-1)) = ((u(2:end,ti)-u(1:(end-1),ti))./Hzw).^2+ ...
    ((v(2:end,ti)-v(1:(end-1),ti))./Hzw).^2;
Sh2(1) = Sh2(2);Sh2(end) = Sh2(end-1);

Ri = N2./Sh2; %Richardson number

%Set kv diffusivity:
kint(Ri<0) = K0_P88_U_m; 
kint(Ri>0) = K0_P88_L_m*Ri(Ri>0).^EX_P88_L_m+ ...
         K0_P88_U_m*(1+Ri(Ri>0)/Ri0_P88_m).^EX_P88_U_m;
kv(:,ti) = kv(:,ti)+kint;

%Set kt & ks diffusivity:
kint(Ri<0) = K0_P88_U_s; 
kint(Ri>0) = K0_P88_L_s*Ri(Ri>0).^EX_P88_L_s+ ...
         K0_P88_U_s*(1+Ri(Ri>0)/Ri0_P88_s).^EX_P88_U_s;
kt(:,ti) = kt(:,ti)+kint;
ks(:,ti) = ks(:,ti)+kint;

%Restrict maximum diffusivity:
kt(kt(:,ti)>P88_Kmax,ti) = P88_Kmax;
ks(ks(:,ti)>P88_Kmax,ti) = P88_Kmax;
kv(kv(:,ti)>P88_Kmax,ti) = P88_Kmax;

end 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY layer KPP Mixing

if (KPPBL)
    
%radiative flux at depth z:
wr = zeros(size(z_w));
for zi=1:length(z_w)
wr(zi)  = srflux(ti)*(1-swdk(z_w(zi)))/Cp/rho0;
end

%Buoyancy flux:
Bf = g*alpha*(wr+wt0)-g*beta*ws0;

%Monin-Obukhov length:
L = Ustar.^3/vonKar./Bf;
L(abs(L)<small) = small*sign(L(abs(L)<small));

%stability parameter:
zeta = -z_w./L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate velocity scales: %%%%%%
wm = zeros(size(z_w));
ws = zeros(size(z_w));
%stable:
cond = zeta>=0;
wm(cond) = vonKar*Ustar./(1+5*zeta(cond));
ws(cond) = wm(cond);
%unstable:
zeta_t = max([zeta epsl*(-KPPMLD)./L],[],2);
cond = zeta_t<0 & zeta_t >= zetam;
wm(cond) = vonKar*Ustar*(1-16*zeta_t(cond)).^(-1/4);
cond = zeta_t < zetam;
wm(cond) = vonKar*Ustar*(am-cm*zeta_t(cond)).^(-1/3);
cond = zeta_t<0 & zeta_t >= zetas;
ws(cond) = vonKar*Ustar*(1-16*zeta_t(cond)).^(-1/2);
cond = zeta_t < zetas;
ws(cond) = vonKar*Ustar*(as-cs*zeta_t(cond)).^(-1/3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate bulk Richardson number:
bw = interp1(z_rho,b(:,ti),z_w,'spline'); %w-points
uw = interp1(z_rho,u(:,ti),z_w,'spline');
vw = interp1(z_rho,v(:,ti),z_w,'spline');
RiKPP_Numer = -(bw(end)-bw).*z_w;
RiKPP_Denom = (uw(end)-uw).^2+(vw(end)-vw).^2+Vtc*(-z_w).*ws.* ...
                             sqrt(abs(N2));
RiKPP = RiKPP_Numer./RiKPP_Denom;
RiKPP(end) = 0;

%Interpolate to find MLD:
KPPMLD = interp1(RiKPP,z_w,Ric,'linear');

%Interpolate to find MO length with new BLD:
LB = interp1(z_w,L,KPPMLD);

%Restrict to be less than ekman and MO depths:
if (EKMO & LB>0)
    KPPMLD = -min([-KPPMLD hekman LB]);
end

%Restrict to be greater than the minimum z_rho:
KPPMLD = min([KPPMLD max(z_rho)]);

%Output:
Hsbl(ti) = KPPMLD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recalculate velocity scales: %%%%%%
wm = zeros(size(z_w));
ws = zeros(size(z_w));
%stable:
cond = zeta>=0;
wm(cond) = vonKar*Ustar./(1+5*zeta(cond));
ws(cond) = wm(cond);
%unstable:
zeta_t = max([zeta epsl*(-KPPMLD)./L],[],2);
cond = zeta_t<0 & zeta_t >= zetam;
wm(cond) = vonKar*Ustar*(1-16*zeta_t(cond)).^(-1/4);
cond = zeta_t < zetam;
wm(cond) = vonKar*Ustar*(am-cm*zeta_t(cond)).^(-1/3);
cond = zeta_t<0 & zeta_t >= zetas;
ws(cond) = vonKar*Ustar*(1-16*zeta_t(cond)).^(-1/2);
cond = zeta_t < zetas;
ws(cond) = vonKar*Ustar*(as-cs*zeta_t(cond)).^(-1/3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc non-dimensional lengths:
sig = z_w/KPPMLD;

%Calculate vertical derivatives of velocity scales:
wmDZ = ((wm(2:end)-wm(1:(end-1)))./Hz);
wmDZ = interp1(z_rho,wmDZ,z_w,'spline');
wsDZ = ((ws(2:end)-ws(1:(end-1)))./Hz);
wsDZ = interp1(z_rho,wsDZ,z_w,'spline');

%Calculate interior diffusivity derivatives:
kvDZ = ((kv(2:end,ti)-kv(1:(end-1),ti))./Hz);
kvDZ = interp1(z_rho,kvDZ,z_w,'spline');
ktDZ = ((kt(2:end,ti)-kt(1:(end-1),ti))./Hz);
ktDZ = interp1(z_rho,ktDZ,z_w,'spline');
ksDZ = ((ks(2:end,ti)-ks(1:(end-1),ti))./Hz);
ksDZ = interp1(z_rho,ksDZ,z_w,'spline');

%Interpolate to MLD:
kvN = interp1(z_w,kv(:,ti),KPPMLD,'linear');
ktN = interp1(z_w,kt(:,ti),KPPMLD,'linear');
ksN = interp1(z_w,ks(:,ti),KPPMLD,'linear');
wmN = interp1(z_w,wm,KPPMLD,'linear');
wsN = interp1(z_w,ws,KPPMLD,'linear');
kvDZN = interp1(z_w,kvDZ,KPPMLD,'linear');
ktDZN = interp1(z_w,ktDZ,KPPMLD,'linear');
ksDZN = interp1(z_w,ksDZ,KPPMLD,'linear');
wmDZN = interp1(z_w,wmDZ,KPPMLD,'linear');
wsDZN = interp1(z_w,wsDZ,KPPMLD,'linear');

%Calculate intermediate shape function coefficients:
G1v = kvN/(-KPPMLD*wmN);
G1t = ktN/(-KPPMLD*wsN);
G1s = ksN/(-KPPMLD*wsN);
G1vDZ = -kvDZN/wmN-kvN*wmDZN/(-KPPMLD*wmN^2);
G1tDZ = -ktDZN/wsN-ktN*wsDZN/(-KPPMLD*wsN^2);
G1sDZ = -ksDZN/wsN-ksN*wsDZN/(-KPPMLD*wsN^2);

%Calculate shape function coefficients:
a2v = -2+3*G1v-G1vDZ;
a2t = -2+3*G1t-G1tDZ;
a2s = -2+3*G1s-G1sDZ;
a3v = 1-2*G1v+G1vDZ;
a3t = 1-2*G1t+G1tDZ;
a3s = 1-2*G1s+G1sDZ;

%Calculate shape function:
Gv = a0+a1*sig+a2v*sig.*sig+a3v.*sig.*sig;
Gt = a0+a1*sig+a2t*sig.*sig+a3t.*sig.*sig;
Gs = a0+a1*sig+a2s*sig.*sig+a3s.*sig.*sig;

%Calculate diffusivities:
kv(sig<1 & sig>0,ti) = (-KPPMLD)*wm(sig<1 & sig>0).*Gv(sig<1 & sig>0);
kt(sig<1 & sig>0,ti) = (-KPPMLD)*ws(sig<1 & sig>0).*Gt(sig<1 & sig>0);
ks(sig<1 & sig>0,ti) = (-KPPMLD)*ws(sig<1 & sig>0).*Gs(sig<1 & sig>0); 

%Add to kt_bl
kt_bl(sig<1 & sig>0,ti) = (-KPPMLD)*ws(sig<1 & sig>0).*Gt(sig<1 & sig>0);


%Calculate non-local transports:
if (nl)
    gams(zeta<0,ti) = Cstar*vonKar*(cs*vonKar*epsl)^(1/3)*ws0./(ws(zeta<0)*(-KPPMLD)+small);
    gamt(zeta<0,ti) = Cstar*vonKar*(cs*vonKar*epsl)^(1/3)*(wt0+wr(zeta< ...
                                                  0))./(ws(zeta<0)*(-KPPMLD)+small);
    gams(sig>1 | sig<0,ti) = 0;
    gamt(sig>1 | sig<0,ti) = 0;

end

end

function dxdt = ZonePcmStateFcnCT_3x(x,u)

% 4 states
%==============
% x(1)=Tin [K]
% x(2)=Tinterior [K]
% x(3)=Tpcm [K]
% x(4)=Tsup [K]

% 4 input
%===============
% u(1)=solrad [W/m2]
% u(2)=Tout [degC]
% u(3)=occ [-]
% u(4)=Vair [m3/h]

m_pcm=1; %kg

solrad=u(1);
Tout=u(2);
occ=u(3);
Vair=u(4);

Tin=x(1);
Tinterior=x(2);
Tpcm=x(3);
Tsup=x(4);

Vi=486.5/4;
Rext=3.0;
Rint=0.05;
tmass=25/4;
imass=17.89/4;
shgc=7.94;
occ_gain=100;
Cp_air= 1005; %J/(kg*K)
rho_air= 1.225; %kg/m3

Re=Rext/(3*(Vi.^(2/3)));
Ri=Rint/(3*(Vi.^(2/3)));
Ce=tmass*1.2*1005*Vi;
Ci=imass*1.2*1005*Vi;

a11 = -1/Re/Ce-1/Ri/Ce;
a12 = 1/Ri/Ce;
a21 = 1/Ri/Ci;
a22 =-1/Ri/Ci;

b11 = shgc/Ce/4;
b12 = 1/Re/Ce;
b13 = occ_gain/Ce/4;
b14 = 1/Ce;

dxdt=zeros(4,1);

if (Tpcm>288.15)&(Tpcm<298.15) % Tpcm>15degC, Tpcm<25degC
    Cp_pcm= 17000; %J/(kg*K)
else
    Cp_pcm= 2000; %J/(kg*K)
end


if (Tin>297.15)&(Tpcm>=Tout) % cooling with ambient air
    q= Cp_air*rho_air*Vair*(Tout-Tin)/3600;
    dxdt(3)= 0;
elseif (Tin<=297.15)&(Tpcm>=Tout) % regeneration
    q=0;
    dxdt(3)= Cp_air*rho_air*Vair*(Tout-Tpcm)/3600/Cp_pcm/m_pcm/96;
else % cooling with PCM or no air supply
    q= Cp_air*rho_air*Vair*(Tsup-Tin)/3600;
    dxdt(3)= Cp_air*rho_air*Vair*(Tout-Tpcm)/3600/Cp_pcm/m_pcm/96;
end  
dxdt(1)= a11*Tin+a12*Tinterior+b11*solrad+b12*Tout+b13*occ+b14*q;
dxdt(2)= a21*Tin+a22*Tinterior;
dxdt(4)=0;


  

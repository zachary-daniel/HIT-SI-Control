clear; close all; clc;

L_s=2.47e-6; %series inductance (with source)
R_s=5e-3; %series resistance (with source)
L_c=100e-9; %capacitor leg inductance
C=95e-6; %capacitance
R_c=2e-3; %capacitor leg resistance
L_cab=0; %cabling inductance
R_cab=0; %cabling resistanceL_p=8.21e-6; %parallel (to load) inductance/primary inductance
R_p=5e-3; %parallel (to load) resistance/primary resistance 
L=1.85e-6; %coil inductance/(transformed) plasma inductance (vacuum)
R=20e-3;
L_p=8.21e-6;
%Tank Circuit Dynx derived quantities
 a=1/(L_s+L_c);
 b=1/(L_cab+L_p+L_c);
 eh=(R_s+R_c)*a;
 f=R_c*a;
 g=L_c*a;
 h=R_c*b;
 k=(R_cab+R_p+R_c)*b;
l=R_p*b;
m=L_c*b;
n=L_p*b;
p=R_p/(L+L_p);
q=(R+R_p)/(L+L_p);
v=L_p/(L+L_p);


 %Tank Circuit Dynx SS (no flux)
 A1=[g*h-eh*(1-n*v),f*(1-n*v)-g*(k-n*p),g*(l-n*q),g*b-a*(1-n*v)];
 A2=[h-m*eh,m*f-k+n*p,l-n*q,b-m*a];
 A3=[v*(h-m*eh),p*(1-m*g)+v*(m*f-k),v*l-q*(1-m*g),v*(b-m*a)];
 A4=[(1/C)*(1-m*g-n*v),(-1/C)*(1-m*g-n*v),0,0];

kfAvac=(1/(1-m*g-n*v)).*[A1;A2;A3;A4];
kfBvac=(1/(1-m*g-n*v)).*[a*(1-n*v);m*a;v*m*a;0];
kfH=[0,1,0,0];
sysvac=ss(kfAvac,kfBvac,kfH,[0]);

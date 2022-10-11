clear; close all; clc;
% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
mdsopen('hitsiu', 220926010);
Amplitude = 600;
Amplitude1 = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7; 
L1 = (8.0141e-7); %H
L2  =2.0462e-6; %H
Lp = L2/2; % H Placeholder value until the impedance of the plasma can be calculated
M = .161*L2; % Coupling coefficient
Mw = .1346*L2;% Coupling coefficient
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
Mp = L2/5; % H This is a placeholder value until an actual measurment of the mutual inductance of the plasma can be measured


%Get inductance of the plasma
i_plasma = mdsvalue('\i_tor_spaavg');
%v_plasma = mdsvalue('\');


%%%State space model that includes the Plasma. A will be 13x13.

%Entries for x3
scaler_x3 = -1/(4*Mp*(L2-Mw)*(L2-2*M+Mw) );
x3a = -3*L2*Mp*R2  + 4 * M * Mp*R2 - Mp*Mw*R2;
x3b = -3 * L2 * Mp+ 4 * M * Mp - Mp * Mw;
x3c = 3 * L2 * Mp * R2 - 4 * M * Mp * R2 + Mp * Mw * R2 + 3 * L2 * Mp * R3 - 4 * M * Mp * R3 + Mp * Mw * R3;
x3d = L2 * Mp * R2 - Mp * Mw * R2;
x3e = L2 * Mp - Mp * Mw;
x3f = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x3g = L2 * Mp * R2 - Mp * Mw * R2;
x3h = L2 * Mp - Mp * Mw;
x3i = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x3j = L2 * Mp * R2- 4 * M * Mp * R2 + 3 * Mp * Mw * R2;
x3k = L2 * Mp  - 4 * M * Mp + 3 * Mp * Mw;
x3l = -L2 * Mp * R2 + 4 * M * Mp * R2 - 3 * Mp * Mw * R2 - L2 * Mp * R3 + 4 * M * Mp * R3 - 3 * Mp * Mw * R3;
x3m = (L2.^2) * Lp - 2 * L2 * Lp * M + 2 * Lp * M * Mw - Lp * (Mw.^2);


%entries for x6
scaler_x6 = scaler_x3;
x6a = L2 * Mp * R2 - Mp * Mw * R2;
x6b = L2 * Mp - Mp * Mw;
x6c = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x6d = -3 * L2 * Mp * R2 + 4 * M * Mp * R2 - Mp * Mw * R2;
x6e = -3 * L2 * Mp + 4 * M * Mp - Mp * Mw;
x6f = 3 * L2 * Mp * R2 - 4 * M * Mp * R2 + Mp * Mw * R2 + 3 * L2 * Mp * R3 - 4 * M * Mp * R3 + Mp * Mw * R3;
x6g = L2 * Mp * R2 - 4 * M * Mp * R2 + 3 * Mp * Mw * R2;
x6h = L2 * Mp - 4 * M * Mp + 3 * Mp * Mw;
x6i = -L2 * Mp * R2 + 4 * M * Mp * R2 - 3 * Mp * Mw * R2 - L2 * Mp * R3 + 4 * M * Mp * R3 - 3 * Mp * Mw * R3;
x6j = L2 * Mp * R2 - Mp * Mw * R2;
x6k = L2 * Mp  - Mp * Mw;
x6l = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x6m = (L2.^2) * Lp  - 2 * L2 * Lp * M + 2 * Lp * M * Mw - Lp * (Mw.^2);

%entries for x9
scaler_x9 = scaler_x3;
x9a = L2 * Mp * R2 - Mp * Mw * R2;
x9b = L2 * Mp - Mp * Mw;
x9c = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x9d = L2 * Mp * R2 - 4 * M * Mp * R2 + 3 * Mp * Mw * R2;
x9e = L2 * Mp  - 4 * M * Mp + 3 * Mp * Mw;
x9f = -L2 * Mp * R2 + 4 * M * Mp * R2 - 3 * Mp * Mw * R2 - L2 * Mp * R3 + 4 * M * Mp * R3 - 3 * Mp * Mw * R3;
x9g = -3 * L2 * Mp * R2 + 4 * M * Mp * R2 - Mp * Mw * R2;
x9h = -3 * L2 * Mp + 4 * M * Mp - Mp * Mw;
x9i = 3 * L2 * Mp * R2 - 4 * M * Mp * R2 + Mp * Mw * R2 + 3 * L2 * Mp * R3 - 4 * M * Mp * R3 + Mp * Mw * R3;
x9j = L2 * Mp * R2 - Mp * Mw * R2;
x9k = L2 * Mp - Mp * Mw;
x9l = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x9m = (L2.^2) * Lp - 2 * L2 * Lp * M + 2 * Lp * M * Mw - Lp * (Mw.^2);

%x12 entries
scaler_x12 = scaler_x3;
x12a = L2 * Mp * R2 - 4 * M * Mp * R2 + 3 * Mp * Mw * R2;
x12b = L2 * Mp - 4 * M * Mp + 3 * Mp * Mw;
x12c = -L2 * Mp * R2 + 4 * M * Mp * R2 - 3 * Mp * Mw * R2 - L2 * Mp * R3 + 4 * M * Mp * R3 - 3 * Mp * Mw * R3;
x12d = L2 * Mp * R2 - Mp * Mw * R2;
x12e = L2 * Mp - Mp * Mw;
x12f = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x12g = L2 * Mp * R2 - Mp * Mw * R2;
x12h = L2 * Mp - Mp * Mw;
x12i = -L2 * Mp * R2 + Mp * Mw * R2 - L2 * Mp * R3 + Mp * Mw * R3;
x12j = -3 * L2 * Mp * R2 + 4 * M * Mp * R2 - Mp * Mw * R2;
x12k = -3 * L2 * Mp + 4 * M * Mp - Mp * Mw;
x12l = 3 * L2 * Mp * R2 - 4 * M * Mp * R2 + Mp * Mw * R2 + 3 * L2 * Mp * R3 - 4 * M * Mp * R3 + Mp * Mw * R3;
x12m = (L2.^2) * Lp - 2 * L2 * Lp * M + 2 * Lp * M * Mw - Lp * (Mw.^2);

%x13 entries
scaler_x13 = -1/(4*(Mp.^2));
x13a = -Mp * R2;
x13b = -Mp;
x13c = Mp * R2 + Mp * R3;
x13d = -Mp * R2;
x13e = -Mp;
x13f = Mp * R2 + Mp * R3;
x13g = -Mp*R2;
x13h = -Mp;
x13i = Mp * R2 + Mp * R3;
x13j = -Mp * R2;
x13k = -Mp;
x13l = Mp * R2  + Mp * R3;
x13m = -L2 * Lp - 2 * Lp * M  - Lp * Mw;

%build A matrix
A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    scaler_x3*x3a, scaler_x3*x3b,scaler_x3*x3c, scaler_x3*x3d, scaler_x3*x3e, scaler_x3*x3f, scaler_x3*x3g, scaler_x3*x3h, scaler_x3*x3i, scaler_x3*x3j, scaler_x3*x3k, scaler_x3*x3l, scaler_x3*x3m;
    0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,  1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0;
    scaler_x6*x6a, scaler_x6*x6b, scaler_x6*x6c, scaler_x6*x6d, scaler_x6*x6e, scaler_x6*x6f, scaler_x6*x6g, scaler_x6*x6h, scaler_x6*x6i, scaler_x6*x6j, scaler_x6*x6k, scaler_x6*x6l, scaler_x6*x6m;
    0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0;
   scaler_x9*x9a, scaler_x9*x9b, scaler_x9*x9c, scaler_x9*x9d, scaler_x9*x9e, scaler_x9*x9f, scaler_x9*x9g, scaler_x9*x9h, scaler_x9*x9i, scaler_x9*x9j, scaler_x9*x9k, scaler_x9*x9l, scaler_x9*x9m;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap;
   scaler_x12*x12a, scaler_x12*x12b, scaler_x12*x12c, scaler_x12*x12d,scaler_x12*x12e, scaler_x12*x12f, scaler_x12*x12g, scaler_x12*x12h, scaler_x12*x12i, scaler_x12*x12j, scaler_x12*x12k, scaler_x12*x12l, scaler_x12*x12m;
    scaler_x13*x13a, scaler_x13*x13b, scaler_x13*x13c, scaler_x13*x13d,scaler_x13*x13e, scaler_x13*x13f, scaler_x13*x13g, scaler_x13*x13h, scaler_x13*x13i, scaler_x13*x13j, scaler_x13*x13k, scaler_x13*x13l, scaler_x13*x13m];

%Build B matrix
B = [1/L1, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 1/L1, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 1/L1, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 1/L1;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0];
%build C matrix where we observe all L2 currents and the plasma current
C = [0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1];

%D matrix is all zeros
D = zeros( size(C,1), size(B,2) );

%build continuous time model
sys_plasma = ss(A,B,C,D);
%discrete time model with dT 
sys_plasma_d = c2d(sys_plasma, dT, 'zoh');
%Create discrete time matrices
Ad = sys_plasma_d.A;
Bd = sys_plasma_d.B;
Cd = sys_plasma_d.C;
Dd = sys_plasma_d.D;

%check ctrb and obsv
control = rank(ctrb(Ad, Bd)) == size(A,1); %gucci
observe = rank(obsv(Ad, Cd)) == size(A,1); %gucci
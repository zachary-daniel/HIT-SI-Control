clear; close all; clc;
%Plasma model for the 4 injector circuits coupled to an inductor with a
%voltage source and a resistor. Will tweak as needed... :) Get this shit
%done

%constants
Amplitude = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7;
dT = SampleTime;
Start = 0;
FormationTime = .002;
L1 = (8.0141e-7); %Henry
L2 = 2.0462e-6; %Henry
M = .161*L2; % Coupling coefficient
Mw = .1346*L2;% Coupling coefficient
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
Lp = L2/5; %Henry
Mp = 1e-10; % Henry
Rp = 1e-13;
Vp = 0;
size_A = 12; %dimension of A matrix


state_variables = [(-R1-R2)/L1,-1/L1,R2/L1,0,0,0,0,0,0,0,0,0;
                   1/Cap,0,-1/Cap,0,0,0,0,0,0,0,0,0;
                   -R2,-1,R3+R2,0,0,0,0,0,0,0,0,0;
                    0,0,0,(-R1-R2)/L1,-1/L1,R2/L1,0,0,0,0,0,0;
                    0,0,0,1/Cap,0,-1/Cap,0,0,0,0,0,0;
                    0,0,0,-R2,-1,R3+R2,0,0,0,0,0,0;
                    0,0,0,0,0,0,(-R1-R2)/L1,-1/L1,R2/L1,0,0,0;
                    0,0,0,0,0,0,1/Cap,0,-1/Cap,0,0,0;
                    0,0,0,0,0,0, -R2,-1,R3+R2,0,0,0;
                    0,0,0,0,0,0,0,0,0,(-R1-R2)/L1,-1/L1,R2/L1;
                    0,0,0,0,0,0,0,0,0,1/Cap,0,-1/Cap;
                    0,0,0,0,0,0,0,0,0,-R2,-1,R3+R2;];

state_derivative_coeff = [1,0,0,0,0,0,0,0,0,0,0,0;
                          0,1,0,0,0,0,0,0,0,0,0,0;
                          0,0,-L2,0,0,-M,0,0,-M,0,0,-Mw;
                          0,0,0,1,0,0,0,0,0,0,0,0;
                          0,0,0,0,1,0,0,0,0,0,0,0;
                          0,0,-M,0,0,-L2,0,0,-Mw,0,0,-M;
                          0,0,0,0,0,0,1,0,0,0,0,0;
                          0,0,0,0,0,0,0,1,0,0,0,0;
                          0,0,-M,0,0,-Mw,0,0,-L2,0,0,-M;
                          0,0,0,0,0,0,0,0,0,1,0,0;
                          0,0,0,0,0,0,0,0,0,0,1,0;
                          0,0,-Mw,0,0,-M,0,0,-M,0,0,-L2];
A = state_derivative_coeff\state_variables;


scalar1 = 1/((L2-Mw)*( (L2.^2) - (4*M.^2) + 2*L2*Mw+ (Mw.^2) )); %Scale factor in front of the entries to the A matrix that are affected by mutual inductance

x3a =  (-L2.^2)*R2+(2*M.^2)*R2-L2*Mw*R2;
x3b = (-L2.^2)+(2*M.^2)-L2*Mw;
x3c = (L2.^2)*R2-(2*M.^2)*R2+L2*Mw*R2+(L2.^2)*R3-(2*M.^2)*R3+L2*Mw*R3;
x3d = L2*M*R2 - M*Mw*R2;
x3e = L2*M-M*Mw;
x3f = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3g = L2*M*R2-M*Mw*R2;
x3h = L2*M-M*Mw;
x3i = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3j = -2*(M.^2)*R2+L2*Mw*R2+(Mw.^2)*R2;
x3k = -2*(M.^2)+L2*Mw+Mw.^2;
x3l = 2*(M.^2)*R2-L2*Mw*R2-(Mw.^2)*R2+2*R3*(M.^2)-L2*Mw*R3-R3*Mw.^2;

%Entries for x6 in A matrix
x6a = -L2*M*R2+M*Mw*R2;
x6b = -L2*M+M*Mw;
x6c = L2*M*R2-M*Mw*R2+L2*M*R3-M*Mw*R3;
x6d = R2*(L2.^2)-2*R2*(M.^2)+L2*Mw*R2;
x6e = (L2.^2)-2*(M.^2)+L2*Mw;
x6f = -R2*(L2.^2)+2*R2*(M.^2)-L2*Mw*R2-R3*(L2.^2)+2*R3*(M.^2)-L2*Mw*R3;
x6g = 2*R2*(M.^2)-L2*Mw*R2-R2*(Mw.^2); 
x6h = 2*(M.^2)-L2*Mw-(Mw.^2);
x6i = -2*R2*(M.^2)+L2*Mw*R2+R2*(Mw.^2)-2*R3*(M.^2)+L2*Mw*R3+R3*(Mw.^2);
x6j = -L2*M*R2+M*Mw*R2;
x6k = -L2*M+M*Mw;
x6l = L2*M*R2-M*Mw*R2+L2*M*R3-M*Mw*R3;

%Entries for x9 in A matrix
x9a = -L2*M*R2 + M*Mw*R2;
x9b = -L2* M + M* Mw;
x9c = L2* M* R2 - M* Mw* R2 + L2* M *R3  - M* Mw* R3 ;
x9d = 2* (M.^2) *R2 - L2*Mw* R2 - (Mw.^2) *R2;
x9e = 2 *(M.^2) - L2* Mw - (Mw.^2);
x9f = -2* (M.^2)* R2+ L2 *Mw *R2 + (Mw.^2) *R2 - 2* (M.^2) *R3 + L2 *Mw* R3 + (Mw.^2) *R3;
x9g =(L2.^2) *R2 - 2* (M.^2) *R2 + L2 *Mw *R2;
x9h = (L2.^2) - 2* (M.^2) + L2 *Mw;
x9i = -(L2.^2) *R2 + 2* (M.^2)* R2 - L2 *Mw *R2 - (L2.^2) *R3 + 2* (M.^2)* R3- L2 *Mw *R3;
x9j = -L2 *M *R2 + M*Mw*R2;
x9k = -L2 *M + M *Mw ;
x9l = L2*M*R2 - M* Mw *R2 + L2 *M *R3 - M* Mw *R3;

%Entries for x12 in A matrix
x12a = -2* (M.^2) *R2+ L2* Mw* R2 +(Mw.^2)* R2;
x12b = -2 *(M.^2) + L2 *Mw + (Mw.^2);
x12c = 2 *(M.^2) *R2 - L2 *Mw *R2 - (Mw.^2) *R2 + 2 *(M.^2) *R3 - L2 *Mw *R3 - (Mw.^2) *R3;
x12d = L2 *M *R2 - M *Mw *R2;
x12e = L2 *M - M *Mw;
x12f = -L2 *M *R2 + M *Mw *R2 - L2 *M *R3 + M *Mw* R3;
x12g = L2 *M *R2 - M *Mw *R2;
x12h = L2 *M - M *Mw;
x12i = -L2 *M *R2 + M *Mw *R2 - L2 *M *R3 + M *Mw *R3;
x12j = (-L2.^2) *R2 + 2 *(M.^2) *R2 - L2 *Mw *R2;
x12k = (-L2.^2) + 2 *(M.^2) - L2* Mw;
x12l = (L2.^2) *R2 - 2 *(M.^2) *R2 + L2 *Mw *R2 + (L2.^2) *R3 - 2*(M.^2)* R3 + L2 *Mw *R3;

A2 = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   -scalar1*x3a, -scalar1*x3b,-scalar1*x3c, -scalar1*x3d, -scalar1*x3e, -scalar1*x3f, -scalar1*x3g, -scalar1*x3h, -scalar1*x3i, -scalar1*x3j, -scalar1*x3k, -scalar1*x3l;
   0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0;
   0, 0, 0,  1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0;
   scalar1*x6a, scalar1*x6b, scalar1*x6c, scalar1*x6d, scalar1*x6e, scalar1*x6f, scalar1*x6g, scalar1*x6h, scalar1*x6i, scalar1*x6j, scalar1*x6k, scalar1*x6l;
   0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap, 0, 0, 0;
   scalar1*x9a, scalar1*x9b, scalar1*x9c, scalar1*x9d, scalar1*x9e, scalar1*x9f, scalar1*x9g, scalar1*x9h, scalar1*x9i, scalar1*x9j, scalar1*x9k, scalar1*x9l;
   0, 0, 0, 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap;
   -scalar1*x12a, -scalar1*x12b, -scalar1*x12c, -scalar1*x12d,-scalar1*x12e, -scalar1*x12f, -scalar1*x12g, -scalar1*x12h, -scalar1*x12i, -scalar1*x12j, -scalar1*x12k, -scalar1*x12l];



clear; close all; clc;
% declare time and voltage for our data
Amplitude = 600;
Frequency = 19000;
RunTime = .004;
SampleTime = 1e-7;
Lp = 1.85e-6; 
L1 = 1.4e-6; %H
L2 = 1.5e-6; %H
C = 95e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 100e-9;
A = [(-1/L1*(R1+R2)), -1/L1, R2*1/L1;
     1/C, 0, -1/C;
     (-1/L2)*R2, 1/L2, (-1/L2)*(R3-R2)];
B = [1/L1;
    0;
    0;];
C = [0 0 1];
D = 0;
states = {'x1', 'x2', 'x3'};
inputs = {'v(t)'}; 
outputs = {'x3'};
%continuous time system
sysc = ss(A,B,C,D, "statename", states,"inputname", inputs, "outputname", outputs);

%Discrete Time System
sys_d = c2d(sysc, dT, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

G = eye(3);
H = zeros(3,3);

Q = diag(.1*ones(1,3));
R = .1;


% LQR Shit
Q_ctr = [1 0 0;
          0 1 0
          0 0 100];
R_ctr = .01;

K = lqr(Ad,Bd ,Q_ctr,R_ctr);
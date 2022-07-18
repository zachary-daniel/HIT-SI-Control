clear;close all;clc;
load("plasma_load_model.mat") % load ideal L2 current and sVduare wave from power supply. 
% declare time and voltage for our data
Amplitude = 600;
Frequency = 19000;
RunTime = .004;
SampleTime = 1e-7;
time = (0:SampleTime:RunTime)';
Lp = 1.85e-6; 
L1 = 1.4e-6; %H
L2 = 1.5e-6; %H
Cap = 95e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 100e-9;
t = 0:dT:RunTime;
A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
     1/Cap, 0, -1/Cap;
     (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2)];
B = [1/L1;
    0;
    0;];
C = [0,0,1];

D = zeros(size(C,1), size(B,2));
states = {'x1', 'x2', 'x3'};
inputs = {'v(t)'}; 
outputs = {'x3'};
%continuous time system
sysc = ss(A,B,C,D, "statename", states,"inputname", inputs, "outputname", outputs);

% Discrete Time System
sys_d = c2d(sysc, dT, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

% Kalman Filter system

Vd = .1*eye(3); % disturbance covariance
Vn = 1; % Noise covariance
BF = [B Vd 0*B];

sysc_noise = ss(A, BF, C, [0 0 0 0 Vn]);

sysFullOutput = ss(A, BF, eye(3), zeros(3, size (BF, 2)),"statename", states,"inputname", inputs, "outputname", outputs);
uDist = randn(3, size(t,2));
uNoise = randn(size(t));
vNoise = [newVoltages'; Vd*Vd*uDist; uNoise];
% sysFullOutput = ss(A,B,eye(3),D, "statename", states,"inputname", inputs, "outputname", outputs);


%Kalman Filter
[Kf, P, E] = lqe(A, Vd, C, Vd, Vn);
sys_kf = ss(A-Kf*C, [B Kf], eye(3), 0*[B Kf]);


% LVdR Shit
% Vd_ctr = [1 0 0;
%           0 1 0
%           0 0 100];
% R_ctr = .01; %Small value of R means that control is cheap



% Setting up set points to approach


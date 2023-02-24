clear; close all; clc;
%Plasma model for the 4 injector circuits coupled to an inductor with a
%voltage source and a resistor. Will tweak as needed... :) Get this shit
%done

%constants
Amplitude = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .002;
SampleTime = 1e-7;
dT = SampleTime;
Start = 0;
L1 = (8.0141e-7); %Henry
L2 = 2.0462e-6; %Henry
M = .161*L2; % Coupling coefficient
Mw = .1346*L2;% Coupling coefficient
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
Lp = L2/5; %Henry
Mp = Lp/2; % Henry
Rp = R2;
Vp = 10;
size_A = 13; %dimension of A matrix
num_inputs = 5;

NoisePower = .1;
%LQR cost matrices
accuracy_penalty = 10000;
Q_cost = diag(ones(size_A,1));
Q_cost(3,3) = accuracy_penalty;
Q_cost(6,6) = accuracy_penalty;
Q_cost(9,9) = accuracy_penalty;
Q_cost(12,12) = accuracy_penalty;
R_cost = .01;

%Kalman Filter Cost matrices
Q = diag(.001*ones(1,num_inputs)); % disturbance covariance
R = diag(1*ones(1,num_inputs)); % Noise covariance


sys_d_vacuum = load("sys_d_vacuum.mat"); %discrete time vacuum circuit model
syskf_vacuum = load('syskf_vacuum.mat'); %discrete timekalman filter for vacuum circuit
K_vacuum = load('lqr_gain_vacuum.mat'); %LQR gain matrix for vacuum model;

%coeff in front of the states x
state_coeff = [-(R1+R2)/L1,-1/L1,R2/L1,0,0,0,0,0,0,0,0,0,0;
               1/Cap,0,-1/Cap,0,0,0,0,0,0,0,0,0,0;
               -R2,-1,R2+R3,0,0,0,0,0,0,0,0,0,0;
               0,0,0,((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0,0;
                0, 0, 0,  1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0,0;
                0,0,0, -R2,-1,R2+R3,0,0,0,0,0,0,0;
                 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0,0;
                  0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap, 0, 0, 0,0;
                  0,0,0,0,0,0,-R2,-1,R2+R3,0,0,0,0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1,0;
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap,0;
                    0,0,0,0,0,0,0,0,0,-R2,-1,R2+R3,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,Rp];

%coeff in front of state derivatives x_dot
state_derivative_coeff = [1,0,0,0,0,0,0,0,0,0,0,0,0;
                          0,1,0,0,0,0,0,0,0,0,0,0,0;
                          0,0,-L2,0,0,-M,0,0,-M,0,0,-Mw,-Mp
                          0,0,0,1,0,0,0,0,0,0,0,0,0;
                          0,0,0,0,1,0,0,0,0,0,0,0,0;
                          0,0,-M,0,0,-L2,0,0,-Mw,0,0,-M,-Mp;
                          0,0,0,0,0,0,1,0,0,0,0,0,0;
                          0,0,0,0,0,0,0,1,0,0,0,0,0;
                          0,0,-M,0,0,-Mw,0,0,-L2,0,0,-M,-Mp;
                          0,0,0,0,0,0,0,0,0,1,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,1,0,0;
                          0,0,-Mw,0,0,-M,0,0,-M,0,0,-L2,-Mp;
                          0,0,-Mp,0,0,-Mp,0,0,-Mp,0,0,-Mp,-Lp];

%compute A matrix
A = (state_derivative_coeff) \ state_coeff;

%B matrix

%coeff in front of inputs for state vector
B =                 [1/L1,0,0,0,0;
                     0,0,0,0,0;
                     0,0,0,0,0;
                     0,1/L1,0,0,0;
                     0,0,0,0,0;
                     0,0,0,0,0;
                     0,0,1/L1,0,0;
                     0,0,0,0,0;
                     0,0,0,0,0;
                     0,0,0,1/L1,0;
                     0,0,0,0,0;
                     0,0,0,0,0;
                     0,0,0,0,-Vp/Lp;];

C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;
     0,0,0,0,0,0,0,0,0,0,0,0,1];

D = zeros(size(C,1),size(B,2));

%Create system
sysc= ss(A,B,C,D);

sys_d = c2d(sysc,dT,'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

%observable
observable = size(A,1) == rank(obsv(Ad,Cd));

%controllable

controllable = size(A,1) == rank(ctrb(Ad,Bd));

time = Start:dT:RunTime;
voltage = Amplitude*sin(2*pi*Frequency*time);

newVoltage = toSquare(voltage, Amplitude,SampleTime,time);

plasma_voltage = Vp*ones(size(time))';

inputs = [newVoltage, newVoltage, newVoltage, newVoltage, plasma_voltage];

out = lsim(sys_d,inputs,time);

%Build Kalman Filter
[kalmf, L, P] = kalman(sysc, Q, R, 0);

syskf_plasma = ss(Ad-L*Cd, [Bd L], eye(size_A), 0*[Bd L], dT);

%Build LQR controller
K_plasma = lqr(sys_d,Q_cost,R_cost);

% 
% figure()
% plot(time,out(:,1)) %L2 flux 1
% title("L2 Current Flux 1")
% xlabel('Time (s)')
% ylabel('Current (Amps)')
% 
% figure()
% plot(time,out(:,5)) %Plasma current
% title('Plasma Current')
% xlabel('Time (s)')
% ylabel('Current (Amps)')
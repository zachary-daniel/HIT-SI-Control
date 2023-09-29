clear; close all; clc;
%Plasma model for the 4 injector circuits coupled to an inductor with a
%voltage source and a resistor. Will tweak as needed... :) Get this shit
%done

%constants

% Start to switch values of circuit components to the voltage circuit
% values instead of the flux circuit values.
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
Mp = .5*sqrt(L2*Lp); % Henry order 150 nH and a resistive load order 30 mOhms
Rp = R2;
Vp = 0;
size_A = 13; %dimension of A matrix
num_inputs = 4;

NoisePower = 0;
%LQR cost matrices
accuracy_penalty = 10000;
Q_cost = diag(ones(size_A,1));
Q_cost(3,3) = accuracy_penalty;
Q_cost(6,6) = accuracy_penalty;
Q_cost(9,9) = accuracy_penalty;
Q_cost(12,12) = accuracy_penalty;
R_cost = .01;

%Kalman Filter Cost matrices
Q = diag(1*ones(1,num_inputs)); % disturbance covariance
R = diag(1*ones(1,num_inputs)); % Noise covariance

%Desired wave for LQR
s = load('desired_L2_wave.mat');
desired_L2_wave = s.L2_Current_Flux_1;


sys_d_vacuum = load("sys_d_vacuum.mat").sys_d; %discrete time vacuum circuit model
syskf_vacuum = load('syskf_vacuum.mat').syskf_vacuum; %discrete timekalman filter for vacuum circuit
K_vacuum = load('lqr_gain_vacuum.mat').K_vacuum; %LQR gain matrix for vacuum model;

% %augment A,B,C,D matrices for vacuum to include an extra state for plasma
% A_vacuum = sys_d_vacuum.A;
% B_vacuum = sys_d_vacuum.B;
% C_vacuum = sys_d_vacuum.C;
% D_vacuum = sys_d_vacuum.D;
% 
% A_vacuum = [A_vacuum zeros(12,1)];
% A_vacuum = [A_vacuum; zeros(1,13)];
% 
% B_vacuum = [B_vacuum zeros(12,1)];
% B_vacuum = [B_vacuum; zeros(1,5)];
% 
% C_vacuum = [C_vacuum zeros(4,1)];
% C_vacuum = [C_vacuum; zeros(1,13)];
% 
% D_vacuum = zeros(5,5);
% 
% %Augment sys_kf_vacuum matrices to include an extra row and col for plasma
% 
% A_kf_vacuum = syskf_vacuum.A;
% B_kf_vacuum = syskf_vacuum.B;
% C_kf_vacuum = syskf_vacuum.C;
% D_kf_vacuum = syskf_vacuum.D;
% 
% C_kf_vacuum = eye(13);
% D_kf_vacuum = zeros(13,10);
% A_kf_vacuum = [A_kf_vacuum zeros(12,1)];
% A_kf_vacuum = [A_kf_vacuum; zeros(1,13)];
% B_kf_vacuum = [B_kf_vacuum zeros(12,2)];
% B_kf_vacuum = [B_kf_vacuum; zeros(1,10)];
% 
% %Augment LQR gain K_vaccum to account for plasma
% 

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
                          0,0,-M,0,0,-Mw,0,0,-L2,0,0,-M, -Mp;
                          0,0,0,0,0,0,0,0,0,1,0,0,0;
                          0,0,0,0,0,0,0,0,0,0,1,0,0;
                          0,0,-Mw,0,0,-M,0,0,-M,0,0,-L2,-Mp;
                          0,0,-Mp,0,0,-Mp,0,0,-Mp,0,0,-Mp,-Lp];

%compute A matrix
A = (state_derivative_coeff) \ state_coeff;

%B matrix

%coeff in front of inputs for state vector
B =                 [1/L1,0,0,0;
                     0,0,0,0;
                     0,0,0,0;
                     0,1/L1,0,0;
                     0,0,0,0;
                     0,0,0,0;
                     0,0,1/L1,0;
                     0,0,0,0;
                     0,0,0,0;
                     0,0,0,1/L1;
                     0,0,0,0;
                     0,0,0,0;
                     0,0,0,0;];

C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;];

D = zeros(size(C,1),size(B,2));

%Create system
sysc= ss(A,B,C,D); 

sys_d_plasma = c2d(sysc,dT,'zoh');

Ad = sys_d_plasma.A;
Bd = sys_d_plasma.B;
Cd = sys_d_plasma.C;
Dd = sys_d_plasma.D;

%observable
observable = size(A,1) == rank(obsv(Ad,Cd));

%controllable

controllable = size(A,1) == rank(ctrb(Ad,Bd));

time = Start:dT:RunTime;
voltage = Amplitude*sin(2*pi*Frequency*time);

newVoltage = toSquare(voltage, Amplitude,SampleTime,time);

plasma_voltage = Vp*ones(size(time))';

inputs = [newVoltage, newVoltage, newVoltage, newVoltage];

[out_plasma,t,x] = lsim(sys_d_plasma,inputs,time);

out_vacuum = lsim(sys_d_vacuum,inputs,time);

%Build Kalman Filter
[kalmf, L, P] = kalman(sys_d_plasma, Q, R, 0);

syskf_plasma = ss(Ad-L*Cd, [Bd L], eye(size_A), 0*[Bd L], dT);

%Build LQR controller
K_plasma = lqr(sys_d_plasma,Q_cost,R_cost);



desired1.signals.values = desired_L2_wave;
desired1.time = time;

desired2.signals.values = phaseShift(desired_L2_wave,90);
desired2.time = time;


desired3.signals.values = phaseShift(desired_L2_wave,180);
desired3.time = time;

desired4.signals.values = phaseShift(desired_L2_wave,270);
desired4.time = time;

time_plasma = (0:dT:(RunTime-FormationTime));

simin.signals.values = desired_L2_wave(1:length(time_plasma),:);
simin.time = time_plasma;

dc_plasma_voltage.signals.values = Vp*(0:dT:(RunTime-FormationTime))';
dc_plasma_voltage.time = time_plasma;

%% Simulate Vacuum Portion
sim('All_injectors_LQG_just_ss_model.slx', 'StopTime', 'FormationTime')
init_vals_plasma = ans.KalmanFilter.signals.values(end,:);


%% Simulate Plasma
sim('All_injectors_LQG_plasma_ss_model.slx', 'StopTime', 'FormationTime')

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

%%
[~,~,x] = lsim(sysc,inputs,time);
figure()
for k = 1:12
    subplot(4,3,k)
    plot(time,x(:,k))
end
figure()
plot(x(:,13))
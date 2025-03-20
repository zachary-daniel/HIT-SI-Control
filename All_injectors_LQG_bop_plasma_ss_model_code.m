clear; close all; clc;
%Code for file All_Injectors_LQG_bop_plasma_ss_model. Fuck these names... 
%This file loads in and defines all of the necessary variables used in the
%simulink model. We're going to b evaluating how well the reduced order bop
%model can be used in an LQG loop

%constants
load("flux_trajectory_arr_220816009.mat")
% Going to use BOP-DMD A matrix here so we will load that in
A_bop_dmd = load("Python_Stuff\data\bop_dmd_A_matrix_220816009.mat").atilde_bop_no_conj;
% Start to switch values of circuit components to the voltage circuit
% values instead of the flux circuit values.
Amplitude = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7;
dT = SampleTime;
Start = 0;
FormationTime = .002;
L1 = (9.78e-7); %Henry
% L2 = 2.0462e-6; %Henry
% M = .161*L2; % Coupling coefficient
% Mw = .1346*L2;% Coupling coefficient
% Cap = 96e-6; % F
% R1 = .0025; %Ohm
% R2 = .005; % Ohm
% R3 = .005;% Ohm
% Lp = L2; %Henry
% Mp = .1*sqrt(L2*Lp); % Henry order 150 nH and a resistive load order 30 mOhms
% Rp = R2;
% Vp = 0;
size_A = 13; %dimension of A matrixR
num_inputs = 4;
r = 5; % 9 modes dmd modes are present, but since tehre are 4 eigenvectors that are 0, the rank of the model is 5
NoisePower = 0;

% Take SVD of the training data to get U matrix which will project to and
% from DMD basis
%DMD model is trained on the last 484 samples of the shot. Time when power
%supplies shut off through the end
% [U,sigma,~] = svd(flux_trajectory_arr(:,end-484+1:end),'econ');
[U,sigma,~] = svd(A_bop_dmd,'econ');
U_r = U(:,1:r); 
% U_r = U(:,1:r); 

%% Based on change of variables, we will change these penalty matrices

%LQR cost matrices
accuracy_penalty = 1000000;
Q_cost = diag(ones(size_A,1));
Q_cost(3,3) = accuracy_penalty;
Q_cost(6,6) = accuracy_penalty;
Q_cost(9,9) = accuracy_penalty;
Q_cost(12,12) = accuracy_penalty;

Q_cost_tilde = U_r'*Q_cost*U_r;
% Q_cost_tilde = eye(r);
% Q_cost_tilde(1,1) = accuracy_penalty;
% for k = 2:r-1
%     Q_cost_tilde(k,k) = accuracy_penalty/10;
% end


%R does not change with a reduced state variable
R_cost = 1;

%Kalman Filter Cost matrices
Q = 1;% disturbance covariance


R = 1; % Noise covariance

%Desired wave for LQR
s = load('desired_L2_wave.mat');
desired_L2_wave = s.L2_Current_Flux_1/3;



%%


%Project BOP matrix into ROM
Atilde_bop_dmd = U_r'*A_bop_dmd*U_r;


%B matrix. Need to resize for ROM. Has to be augmented to include a lot of
%zeros to account for augmented state and input vectors.
%From here onwards, will use 'tilde' to indicate augmented matrices

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
%B_bar = [B, zeros(13,9)];


C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;];
%C_bar = [C; zeros(9,13)];

D = zeros(size(C,1),size(B,2));

% We will now transform all of our matrices into the POD/DMD basis
% D_bar_tilde = U_r'*D_bar*U_r;
Btilde = U_r'*B;
Ctilde = C*U_r;



%Create system matrices
sysc= ss(Atilde_bop_dmd,Btilde,Ctilde,D); 

sys_d_plasma = c2d(sysc,dT,'zoh');

Ad = sys_d_plasma.A;
Bd = sys_d_plasma.B;
Cd = sys_d_plasma.C;
Dd = sys_d_plasma.D;

%observable
observable = size(Atilde_bop_dmd,1) == rank(obsv(Ad,Cd));

%controllable

controllable = size(Atilde_bop_dmd,1) == rank(ctrb(Ad,Bd));

time = Start:dT:RunTime;
voltage = Amplitude*sin(2*pi*Frequency*time);

newVoltage = toSquare(voltage, Amplitude,SampleTime,time);


inputs = [newVoltage, newVoltage, newVoltage, newVoltage]; 
     

sys_full = ss(real(A_bop_dmd), B,C,D);
[out_plasma,t,x] = lsim(sys_d_plasma,inputs,time);

[sys_bar,T] = minreal(sys_full, 1e-10);

% Abar = T*A*T';
% Bbar = T*B;
% Cbar = C*T';
% [U_c,S_c,V_c] = svd(ctrb(sys_bar.A,sys_bar.B),'econ');
% 
% Abar = T*real(A_bop_dmd)*T';
% Bbar = T*B;
% Cbar = C*T';
% 
% 
% [kalmf_bar, L, P] = kalman(sys_bar,Q,R,0);
%%
%Build Kalman Filter
[kalmf, L, P] = kalman(sys_d_plasma, Q, R, 0);

syskf_plasma = ss(Ad-L*Cd, [Bd L], eye(r), 0*[Bd L], dT);

%Build LQR controller
K_plasma = dlqr(sys_d_plasma.A, sys_d_plasma.B ,Q_cost_tilde,R_cost, 0);

%Need to also change our desired signals as we would for our LQR equations

desired1.signals.values = 0;
desired1.time = time;

%Declare simin object for desired signal
desired1tilde.signals.values=0;
desired1tilde.time=0;
% desired2.signals.values = phaseShift(desired_L2_wave,90);
desired2.signals.values = 0;
desired2.time = time;

%Declare simin object for desired signal
desired2tilde.signals.values=0;
desired2tilde.time=0;
% desired3.signals.values = phaseShift(desired_L2_wave,180);
desired3.signals.values = 0;
desired3.time = time;

%Declare simin object for desired signal
desired3tilde.signals.values=0;
desired3tilde.time=0;
% desired4.signals.values = phaseShift(desired_L2_wave,270);
desired4.signals.values = 0;
desired4.time = time;

%Declare simin object for desired signal
desired4tilde.signals.values=0;
desired4tilde.time=0;

%Declare simin object for desired signal
desired5tilde.signals.values=0;
desired5tilde.time=0;

% %  

%Stack desired signals into an array and project into reduced state
%dimension
desired_signals_arr = [zeros(size(desired_L2_wave)),zeros(size(desired_L2_wave)), desired_L2_wave, ...
zeros(size(desired_L2_wave)), zeros(size(desired_L2_wave)), desired_L2_wave, zeros(size(desired_L2_wave)), zeros(size(desired_L2_wave))...
    desired_L2_wave, zeros(size(desired_L2_wave)), zeros(size(desired_L2_wave)),desired_L2_wave,zeros(size(desired_L2_wave))];

desired_signals_arr_tilde = (U_r'*desired_signals_arr')'; %Desired signals in POD basis

time_plasma = (0:dT:(RunTime-FormationTime));

simin.signals.values = desired_L2_wave(1:length(time_plasma),:);
simin.time = time_plasma;

desiredtilde_struct = {desired1tilde, desired2tilde, desired3tilde, desired4tilde, desired5tilde};

%Larger desired waveform If I want to try LQG with a full state vector
% desiredtilde_struct = {desired1tilde, desired2tilde, desired3tilde, desired4tilde, desired5tilde, desired6tilde, desired7tilde, 
%     desired8tilde, desired9tilde,desired10tilde,desired11tilde, desired12tilde, desired13tilde};


for k = 1:length(desiredtilde_struct)
    desiredtilde_struct{k}.signals.values = real(desired_signals_arr_tilde(:,k)); %Take only the real part of the desired signal
    desiredtilde_struct{k}.time = time;
end


%% Plot desired wave
% figure()
% plot(time,desired1.signals.values)



%% Check if Kalman filter system and plasma system agree
% [xout_plasma,~,yout_plasma] = lsim(sys_d_plasma, inputs, time);
% 
% [~,~,xout_kf] = lsim(syskf_plasma, [inputs,xout_plasma], t);
% 
% for k = 1:r
%     figure()
%     plot(time,xout_kf(:,k)')
%     hold on
%     plot(time,yout_plasma(:,k)')
%     legend('kalman', 'true')
% end

%% Run plasma simulation with LQG
sim("All_injectors_LQG_bop_plasma_ss_model.slx")
%% Analysis of LQG loop
%need to convert the outputs of the simulation that are in the reduced
%basis to the full order basis
kf_outputs_tilde = (ans.KalmanFilter.signals.values); % Pull kalman filter outputs. Since the system is observed through a C matrix with
ss_outputs_tilde = (ans.StateSpaceModel.signals.values);
lqr_gains = ans.LQRoutputs.signals.values; %LQR outputs from the simulation

% off diagonal entries, we need to undo the affects of this matrix on the
% system. 
%% Uncomment later
% ss_outputs = (U_r*ss_outputs_tilde)'; %transpose for plotting purposes
%% Compare Kalman filter outputs with true model outputs from simulink
% for k = 1:r
%     figure()
%     plot(ans.KalmanFilter.signals.values(:,k))
%     hold on
%     plot(ss_outputs_tilde(k,:)')
%     legend('kalman', 'true output')
% end
%% Desired signals
figure()
for k = 1:5
    plot(desiredtilde_struct{1,k}.signals.values, 'DisplayName', sprintf('%f', (k)))
    hold on
end

%% Desired vs. actual signals
figure()
for k = 1:5
    subplot(5,1,k)
    plot(kf_outputs_tilde(:,k), 'LineWidth',1)
    hold on
    plot(desiredtilde_struct{1,k}.signals.values)
    legend('ss outputs', 'desired')

end

%% LQR outputs
figure('Name', 'LQR outputs')
for k = 1:4
   subplot(5,1,k)
   plot(time,ans.LQRoutputs.signals.values(:,k))
end

sysd_full = ss(A_bop_dmd, B, C, D);

%% Full state desired signals
desired_regala = real((U_r*desired_signals_arr_tilde'))'; %Garbage to get
% this 
% desired_regala = real((pinv(U_r')*kf_outputs_tilde')'); Similar garbage
modified_lqr_gains = [lqr_gains(:,3),lqr_gains(:,3),lqr_gains(:,3),lqr_gains(:,3)];

% [~,~,desired_regala] = lsim(sysd_full,modified_lqr_gains, time);
figure('Name','LQR outputs on full order system')
plot(time, desired_regala(:,3))
hold on
plot(time, desired_signals_arr(:,3))
legend('outputs of lqg loop projected into full state space','desired waveform in full state');

%%



clear; close all; clc;
%load("All_injectors_vaccum_full_workspace.mat")
%load("DMDc_on_old_shot.mat")
load("All_injectors_vaccum_full_workspace.mat")
% sys_dmd_on_shot_c = sys_dmd;
% Let's first check that DMDc can somewhat accurately recreate the A and B
% matrices that come out of the state space model I derived
dT = 1e-7;
time = 0:dT:.004;
B = sysc.B;
Q = .001; %diag(.001*ones(1,size(A, 1))); % disturbance covariance
R = diag(10*ones(1,size(B,2))); % Noise covariance

%% Stuff for single injector model
% Ad = sys_d.A;
% Bd = sys_d.B;
% Cd = sys_d.C;
% Dd = sys_d.D;
% [kalmf, L, P] = kalman(sys_d, Q, R, 0);
% syskf = ss(Ad-L*Cd, [Bd L],eye(3), 0*[Bd L], dT);
% [y,t] = lsim(sys_d, [newVoltages], time);
% [yk, tk] = lsim(syskf, [newVoltages, y], time);
% 
% L1_sim = yk(:,1);
% C_sim = yk(:,2);
% L2_sim = yk(:,3);


%%


% L1_Current_Flux_4 = ans.L1CurrentFlux4.signals.values; % catch L1 Current in Flux 4
% L2_Current_Flux_3 = ans.L2CurrentFlux3.signals.values; % catch L2 Current in Flux 3
% L2_Current_Flux_4 = ans.L2CurrentFlux4.signals.values; % catch L2 current in Flux 4
% C_Voltage_Flux_3 = ans.CVoltageFlux3.signals.values; %Catch C voltage from flux 3
% C_Voltage_Flux_4 = ans.CVoltageFlux4.signals.values; % Catch C voltage from flux 4
% First construct the X matrix that goes from time = 1 to time = end-1
% stack in all measurments

% B = [1/L1, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 1/L1, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 1/L1, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 1/L1;
%     0, 0, 0, 0;
%     0, 0, 0, 0;];



% X = [L1_Current_Flux_1(1:end-1,1)'; L1_Current_Flux_2(1:end-1,1)'; L1_Current_Flux_3(1:end-1,1)'; L1_Current_Flux_4(1:end-1,1)';
%     C_Voltage_Flux_1(1:end-1,1)'; C_Voltage_Flux_2(1:end-1,1)'; C_Voltage_Flux_3(1:end-1,1)'; C_Voltage_Flux_4(1:end-1,1)';
%    L2_Current_Flux_1(1:end-1,1)'; L2_Current_Flux_2(1:end-1,1)'; L2_Current_Flux_3(1:end-1,1)'; L2_Current_Flux_4(1:end-1,1)'];
% % data stacked from the second position to the end
% X2 = [L1_Current_Flux_1(2:end,1)'; L1_Current_Flux_2(2:end,1)'; L1_Current_Flux_3(2:end,1)'; L1_Current_Flux_4(2:end,1)';
%     C_Voltage_Flux_1(2:end,1)'; C_Voltage_Flux_2(2:end,1)'; C_Voltage_Flux_3(2:end,1)'; C_Voltage_Flux_4(2:end,1)';
%     L2_Current_Flux_1(2:end,1)'; L2_Current_Flux_2(2:end,1)'; L2_Current_Flux_3(2:end,1)'; L2_Current_Flux_4(2:end,1)';];

% Gonna try to do this as if there is only one injector circuit firing
%HOLY SHIT IT WORKED!!!!
% X = [L1_Current(1:end-1,1)'; C_Voltage(1:end-1,1)'; L2_Current(1:end-1,1)'];
% 
% X2 = [L1_Current(2:end,1)'; C_Voltage(2:end,1)'; L2_Current(2:end,1)'];
% Upsilon = [newVoltages(1:end-1,1)'];
% 
% Omega = [X;Upsilon];
% 
% [Utilde, Stilde Vtilde] = svd(Omega, 'econ');
% 
% U1tilde = Utilde(1:3,:);
% U2tilde = Utilde(4,:);
% 
% [Uhat, Shat, Vhat] = svd(X2, 'econ');
% 
% 
% Atilde = X2*Vtilde*pinv(Stilde)*U1tilde';
% 
% Btilde = X2*Vtilde*pinv(Stilde)*U2tilde';
% 
% sys_dmd_from_ss = ss(Atilde, Btilde, Cd, Dd, dT);

%%

% Omega = [X; Upsilon]; %omega matrix with the data and the inputs stacked on each other
% 
% %perform SVD on Omega
% 
% [Utilde, Stilde, Vtilde] = svd(Omega, 'econ'); %SVD omega matrix
% 
% [Uprime, Sprime, Vprime] = svd(X2, 'econ'); %Svd snapshots from the second half of the data
% 
% U1tilde = Utilde( 1:length(X(:,1)) , :); % U1 tilde which captures how the states progress and will be used for the A matrix
% U2tilde = Utilde( length( U1tilde(:,1))+1:end ,:); % U2 tilde which captures how inputs propegate and will be used for the B matrix. 
%  
% Atilde = Uprime'*X2*Vtilde*inv(Stilde)*U1tilde'*Uprime; %construct Atilde matrix which should be an approximation of the dynamics
%                                  % Eigs of A will show temporal behavior
% 
% 
% Btilde = Uprime'*X2*Vtilde*inv(Stilde)*U2tilde';
% 
% % sys_dmd_on_sim = ss((X2-Bd*Upsilon)*pinv(X), Btilde, C, D,dT); % this is
% % the case where B is known
% sys_dmd_on_sim = ss(Atilde, Btilde, Cd,Dd,dT);
% 
% 
% lsim(sys_dmd_on_sim, [newVoltages], time); 
% 

%% Let's try doing this on a small sample of the shot data to do a quick time prediction opposed to trying to capture all of the dynamics
start = 1;
stop = 3500; % As long as the stop sample is after the circuit goes to steady state this thing works like a charm

[kalmf, L, P] = kalman(sys_d, Q, R, 0);

syskf = ss(Ad-L*Cd, [Bd L],eye(12), 0*[Bd L], dT);
[y,t] = lsim(sys_d, [newVoltages newVoltages newVoltages newVoltages], time);
[yk, tk] = lsim(syskf, [newVoltages newVoltages newVoltages newVoltages, y], time);

[newVoltageShift1] = phaseShift(newVoltages, PhaseAngle1, loc_nada);
[newVoltageShift2] = phaseShift(newVoltages, PhaseAngle2, loc_nada);
[newVoltageShift3] = phaseShift(newVoltages, PhaseAngle3, loc_nada);

%Replace data from simulink with data from Kalman filter. If this is
%commented out, data from Simulink is going to be used instead

% L1_Current_Flux_1 = yk(:,1);
% C_Voltage_Flux_1 = yk(:,2);
% L2_Current_Flux_1 = yk(:,3);
% L1_Current_Flux_2 = yk(:,4);
% C_Voltage_Flux_2 = yk(:,5);
% L2_Current_Flux_2 = yk(:,6);
% L1_Current_Flux_3 = yk(:,7);
% C_Voltage_Flux_3 = yk(:,8);
% L2_Current_Flux_3 = yk(:,9);
% L1_Current_Flux_4 = yk(:,10);
% C_Voltage_Flux_4 = yk(:,11);
% L2_Current_Flux_4 = yk(:,12);


X = [L1_Current_Flux_1(start:stop-1,1)'; C_Voltage_Flux_1(start:stop-1,1)'; L2_Current_Flux_1(start:stop-1,1)'; L1_Current_Flux_2(start:stop-1,1)';
    C_Voltage_Flux_2(start:stop-1,1)'; L2_Current_Flux_2(start:stop-1,1)'; L1_Current_Flux_3(start:stop-1,1)'; C_Voltage_Flux_3(start:stop-1,1)';
   L2_Current_Flux_3(start:stop-1,1)'; L1_Current_Flux_4(start:stop-1,1)'; C_Voltage_Flux_4(start:stop-1,1)'; L2_Current_Flux_4(start:stop-1,1)'];
% data stacked from the second position to the end
X2 = [L1_Current_Flux_1(start+1:stop,1)'; C_Voltage_Flux_1(start+1:stop,1)'; L2_Current_Flux_1(start+1:stop,1)'; L1_Current_Flux_2(start+1:stop,1)'; 
    C_Voltage_Flux_2(start+1:stop,1)'; L2_Current_Flux_2(start+1:stop,1)'; L1_Current_Flux_3(start+1:stop,1)'; C_Voltage_Flux_3(start+1:stop,1)';
   L2_Current_Flux_3(start+1:stop,1)'; L1_Current_Flux_4(start+1:stop,1)'; C_Voltage_Flux_4(start+1:stop,1)'; L2_Current_Flux_4(start+1:stop,1)';];


Upsilon = [newVoltages(start:stop-1,1)'; newVoltageShift1(start:stop-1,1)'; newVoltageShift2(start:stop-1,1)'; newVoltageShift3(start:stop-1,1)';];

Omega = [X; Upsilon];

[Utilde, Stilde, Vtilde] = svd(Omega, 'econ');

[Uhat, Shat, Vhat] = svd(X2, 'econ');

U1tilde = Utilde(1:12,:);

U2tilde = Utilde(13:16, :);
%Abar and Bbar are the reconstructed matrices without a truncation in the
%SVD
Abar = X2*Vtilde*pinv(Stilde)*U1tilde';

Bbar = X2*Vtilde*pinv(Stilde)*U2tilde';

sys_dmd = ss(Abar, Bbar, Cd, Dd, dT);

% Now let's try with truncation

trunc1 =4; %size(Stilde,1);

trunc2 = 7;

% To truncate in SVD we take the first r cols of U, an rxr square of S, and
% the first r cols of V
Utilde_t = Utilde(:,1:trunc1);
Stilde_t = Stilde(1:trunc1, 1:trunc1);
Vtilde_t = Vtilde(:,1:trunc1); 
U1tilde_t = Utilde_t(1:12,:);
U2tilde_t = Utilde_t(13:16,:);

%%
Abar_t = X2*Vtilde_t*pinv(Stilde_t)*U1tilde_t';
Bbar_t = X2*Vtilde_t*pinv(Stilde_t)*U2tilde_t';
sys_dmd_t = ss(Abar_t, Bbar_t, Cd, Dd, dT);
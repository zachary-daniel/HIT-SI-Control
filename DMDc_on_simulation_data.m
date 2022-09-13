clear; close all; clc;
%load("All_injectors_vaccum_full_workspace.mat")
%load("DMDc_on_old_shot.mat")
load("Vaccum_load_model.mat")
% sys_dmd_on_shot_c = sys_dmd;
% Let's first check that DMDc can somewhat accurately recreate the A and B
% matrices that come out of the state space model I derived
dT = 1e-7;
time = 0:dT:.004;
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

X = [L1_Current(1:end-1,1)'; C_Voltage(1:end-1,1)'; L2_Current(1:end-1,1)'];

X2 = [L1_Current(2:end,1)'; C_Voltage(2:end,1)'; L2_Current(2:end,1)'];


Upsilon = [newVoltages(1:end-1,1)'];

Omega = [X; Upsilon]; %omega matrix with the data and the inputs stacked on each other

%perform SVD on Omega

[Utilde, Stilde, Vtilde] = svd(Omega, 'econ'); %SVD omega matrix

[Uprime, Sprime, Vprime] = svd(X2, 'econ'); %Svd snapshots from the second half of the data

U1tilde = Utilde( 1:length(X(:,1)) , :); % U1 tilde which captures how the states progress and will be used for the A matrix
U2tilde = Utilde( length( U1tilde(:,1))+1:end ,:); % U2 tilde which captures how inputs propegate and will be used for the B matrix. 
 
Atilde = Uprime'*X2*Vtilde*pinv(Stilde)*U1tilde'*Uprime; %construct Atilde matrix which should be an approximation of the dynamics
                                 % Eigs of A will show temporal behavior


Btilde = Uprime'*X2*Vtilde*pinv(Stilde)*U2tilde';

% sys_dmd_on_sim = ss((X2-Bd*Upsilon)*pinv(X), Btilde, C, D,dT); % this is
% the case where B is known
sys_dmd_on_sim = ss(Atilde, Btilde, Cd,Dd,dT);


lsim(sys_dmd_on_sim, [newVoltages], time); % this shit ain't accurate


%% Let's try doing this on a small sample of the shot data to do a quick time prediction opposed to trying to capture all of the dynamics
start = 1;
stop = start + 128;


X = [L1_Current_Flux_1(start:stop-1,1)'; L1_Current_Flux_2(start:stop-1,1)'; L1_Current_Flux_3(start:stop-1,1)'; L1_Current_Flux_4(start:stop-1,1)';
    C_Voltage_Flux_1(start:stop-1,1)'; C_Voltage_Flux_2(start:stop-1,1)'; C_Voltage_Flux_3(start:stop-1,1)'; C_Voltage_Flux_4(start:stop-1,1)';
   L2_Current_Flux_1(start:stop-1,1)'; L2_Current_Flux_2(start:stop-1,1)'; L2_Current_Flux_3(start:stop-1,1)'; L2_Current_Flux_4(start:stop-1,1)'];
% data stacked from the second position to the end
X2 = [L1_Current_Flux_1(start+1:stop,1)'; L1_Current_Flux_2(start+1:stop,1)'; L1_Current_Flux_3(start+1:stop,1)'; L1_Current_Flux_4(start+1:stop,1)';
    C_Voltage_Flux_1(start+1:stop,1)'; C_Voltage_Flux_2(start+1:stop,1)'; C_Voltage_Flux_3(start+1:stop,1)'; C_Voltage_Flux_4(start+1:stop,1)';
    L2_Current_Flux_1(start+1:stop,1)'; L2_Current_Flux_2(start+1:stop,1)'; L2_Current_Flux_3(start+1:stop,1)'; L2_Current_Flux_4(start+1:stop,1)';];


Upsilon = [newVoltages(start:stop-1,1)'; newVoltages(start:stop-1,1)'; newVoltages(start:stop-1,1)'; newVoltages(start:stop-1,1)';];

Omega = [X; Upsilon]; %omega matrix with the data and the inputs stacked on each other

%perform SVD on Omega

[Utilde, Stilde, Vtilde] = svd(Omega, 'econ'); %SVD omega matrix

[Uprime, Sprime, Vprime] = svd(X2, 'econ'); %Svd snapshots from the second half of the data

U1tilde = Utilde(1:12, :); % U1 tilde which captures how the states progress and will be used for the A matrix
U2tilde = Utilde(13:16,:); % U2 tilde which captures how inputs propegate and will be used for the B matrix. 

Atilde = Uprime'*X2*Vtilde*pinv(Stilde)*U1tilde'*Uprime; %construct Atilde matrix which should be an approximation of the dynamics
                                 % Eigs of A will show temporal behavior


Btilde = Uprime'*X2*Vtilde*pinv(Stilde)*U2tilde';

sys_dmd_on_sim = ss(Atilde, Btilde, C, D, dT);

lsim( sys_dmd_on_sim, [newVoltages(stop:stop+stop) newVoltages(stop:stop+stop) newVoltages(stop:stop+stop) newVoltages(stop:stop+stop)], time(stop:stop+stop), X(:,end));

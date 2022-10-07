clear; close all; clc;
mdsconnect('landau.hit')
mdsopen('hitsiu', 220802016) %Shot in vacuum
load('All_injectors_vaccum_full_workspace.mat')
%%
time = mdsvalue('dim_of(\v_div_1_fspa)');
time = time(4095:end); %get rid of negative time values
C = [0 ,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

D = zeros( size(C,1), 4 );
dT = time(2)- time(1);

% L2 currents
L2_Current_Flux_1 = mdsvalue('\i_fcoil_1');
L2_Current_Flux_2 = mdsvalue('\i_fcoil_2');
L2_Current_Flux_3 = mdsvalue('\i_fcoil_3');
L2_Current_Flux_4 = mdsvalue('\i_fcoil_4');

%L1 Currents
L1_Current_Flux_1 = mdsvalue('\i_spa_f1');
L1_Current_Flux_2 = mdsvalue('\i_spa_f2');
L1_Current_Flux_3 = mdsvalue('\i_spa_f3');
L1_Current_Flux_4 = mdsvalue('\i_spa_f4');

%Reduce number of samples of L1 Currents
L1_Current_Flux_1 = boxCarAvg(L1_Current_Flux_1, 2);
L1_Current_Flux_2 = boxCarAvg(L1_Current_Flux_2, 2);
L1_Current_Flux_3 = boxCarAvg(L1_Current_Flux_3, 2);
L1_Current_Flux_4 = boxCarAvg(L1_Current_Flux_4, 2);


%Reduce number of Samples of L2 currents

L2_Current_Flux_1 = boxCarAvg(L2_Current_Flux_1, 2);
L2_Current_Flux_2 = boxCarAvg(L2_Current_Flux_2, 2);
L2_Current_Flux_3 = boxCarAvg(L2_Current_Flux_3, 8);
L2_Current_Flux_4 = boxCarAvg(L2_Current_Flux_4, 2);


%Pseudo C Voltage. L2 is in parallel with C so the voltage should be the
%same
C_Voltage_Flux_1 = mdsvalue('\v_div_1_fc');
C_Voltage_Flux_2 = mdsvalue('\v_div_2_fc');
C_Voltage_Flux_3 = mdsvalue('\v_div_3_fc');
C_Voltage_Flux_4 = mdsvalue('\v_div_4_fc');

% SPA waveforms
SPA_Voltage_1 = mdsvalue('\v_div_1_fspa');
SPA_Voltage_2 = mdsvalue('\v_div_2_fspa');
SPA_Voltage_3 = mdsvalue('\v_div_3_fspa');
SPA_Voltage_4 = mdsvalue('\v_div_4_fspa');



%% DMD
start = 4095; %start at t = 0
stop = length(L2_Current_Flux_1);

X = [L1_Current_Flux_1(start:stop-1,1)'; C_Voltage_Flux_1(start:stop-1,1)'; L2_Current_Flux_1(start:stop-1,1)'; L1_Current_Flux_2(start:stop-1,1)';
    C_Voltage_Flux_2(start:stop-1,1)'; L2_Current_Flux_2(start:stop-1,1)'; L1_Current_Flux_3(start:stop-1,1)'; C_Voltage_Flux_3(start:stop-1,1)';
   L2_Current_Flux_3(start:stop-1,1)'; L1_Current_Flux_4(start:stop-1,1)'; C_Voltage_Flux_4(start:stop-1,1)'; L2_Current_Flux_4(start:stop-1,1)'];

X2 = [L1_Current_Flux_1(start+1:stop,1)'; C_Voltage_Flux_1(start+1:stop,1)'; L2_Current_Flux_1(start+1:stop,1)'; L1_Current_Flux_2(start+1:stop,1)'; 
    C_Voltage_Flux_2(start+1:stop,1)'; L2_Current_Flux_2(start+1:stop,1)'; L1_Current_Flux_3(start+1:stop,1)'; C_Voltage_Flux_3(start+1:stop,1)';
   L2_Current_Flux_3(start+1:stop,1)'; L1_Current_Flux_4(start+1:stop,1)'; C_Voltage_Flux_4(start+1:stop,1)'; L2_Current_Flux_4(start+1:stop,1)';];

Upsilon = [SPA_Voltage_1(start:stop-1,1)'; SPA_Voltage_2(start:stop-1,1)'; SPA_Voltage_3(start:stop-1,1)'; SPA_Voltage_4(start:stop-1,1)';];

Omega = [X;Upsilon];
[Utilde, Stilde, Vtilde] = svd(Omega, 'econ');

U1tilde = Utilde(1:size(X,1),:);
U2tilde = Utilde(size(X,1)+1:size(Utilde,1),:);

A = X2*Vtilde*pinv(Stilde)*U1tilde';
B = X2*Vtilde*pinv(Stilde)*U2tilde';

sys = ss(A,B,C,D, dT);


A_sub_matrices = zeros(100*size(A,1), size(A,1));
B_sub_matrices = zeros(100*size(B,1), size(B,2));

A_avg = zeros(size(A));
B_avg = zeros(size(B));

Avg_eig_A = zeros(size(12,1));

trunc1 = 4; %Truncation points to make matrix algebra faster

trunc2 = 12;
num_iters = 0;

for i = 1:1000
    X_sub_matrix = X(:, i:5:length(time)-1); %Build sub matrices of X, X2, and Upsilon
    X2_sub_matrix = X2(:, i:5:length(time)-1);
    Upsilon_sub_matrix = Upsilon(:, i:5:length(time)-1);
  
    Omega_sub = [X_sub_matrix;Upsilon_sub_matrix;]; %Stack X and Upsilon
      %perform SVD
    [Utilde_sub, Stilde_sub, Vtilde_sub] = svd(Omega_sub, 'econ');
   % [Uhat_sub, Shat_sub, Vhat_sub] = svd(X_sub_matrix, 'econ');
    %Break matrices up to truncation points
     Utilde_sub_t= Utilde_sub(:,1:trunc1);
     Stilde_sub_t = Stilde_sub(1:trunc1, 1:trunc1);
     Vtilde_sub_t = Vtilde_sub(:,1:trunc1); 
    %Break up Utilde_sub into U1 and U2tilde. U1 is the dyanamics of the
    %input and U2 is the dynamics of the input
    U1tilde_sub = Utilde_sub(1:12,:);
    U2tilde_sub = Utilde_sub(13:16,:);
    
    A_sub = X2_sub_matrix*Vtilde_sub*pinv(Stilde_sub)*U1tilde_sub'; %construct A and B matrices
    B_sub = X2_sub_matrix*Vtilde_sub*pinv(Stilde_sub)*U2tilde_sub';

    Avg_eig_A = eig(A_sub) + Avg_eig_A;

    start_point = (num_iters*12)+1;
    len_a = length(A_sub);
    len_b = length(B_sub);

    A_sub_matrices(start_point:start_point+(len_a-1) ,:) = A_sub;
    B_sub_matrices(start_point:start_point+(len_b-1) ,:) = B_sub;
    A_avg = A_avg + A_sub;
    B_avg = B_avg + B_sub;
    num_iters = num_iters+1;
   % sys_temp = ss(A_sub, B_sub, Cd, Dd, dT);

end

Avg_eig_A = Avg_eig_A/i;
A_avg = A_avg/i;
B_avg = B_avg/i;
sys_temp = ss(A_avg, B_avg, C, D, dT);


plot(eig(A), '*', "LineStyle", "none")
title("Eigenvalues of A matrix from DMDc on Vacuum Shot")
xlabel('Real part')
ylabel('Imaginary part')
hold on
plot(eig(Ad), 'r*', 'LineStyle','none')
legend('Eigenvalues from DMDc', "Eigenvalues from statespace model")
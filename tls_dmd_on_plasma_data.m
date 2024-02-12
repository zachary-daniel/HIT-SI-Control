clear;close all;clc;
load('shots\high_itor_plasma_shots\matlab data\221129065_plasma_flux.mat')
%Two good shots are 220816009,221129011
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')

%%

flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,...
                      i_L1_2,v_cap_2,i_fcoil_2,...
                      i_L1_3,v_cap_3,i_fcoil_3,...
                      i_L1_4,v_cap_4,i_fcoil_4];
num_rows = size(flux_trajectory_arr,2);
 %Alright so my knowledge of the system is that it is still obsv and ctrb with only
%measurments of the flux coil.
%I'm gonna neglect the cap
%voltages because I don't really trust them tbh
%   flux_trajectory_arr = [i_fcoil_1,...                                          
%                       i_fcoil_2,...
%                       i_fcoil_3,...
%                       i_fcoil_4];


time = time';


[~,t_equal_0] = min(abs(time));
time = time(:,t_equal_0:end);
dT = time(2)-time(1); %time step for discrete time model
flux_trajectory_arr = flux_trajectory_arr';
flux_trajectory_arr = flux_trajectory_arr(:,t_equal_0:end);

inputs = [v_spa_1,v_spa_2,v_spa_3,v_spa_4]';
inputs = inputs(:,t_equal_0:end);
%% Build B, C, and D matrices
L1 = (8.0141e-7); %Henry. Could this matrix change when we enter plasma?
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
B = B(1:num_rows,:);


C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;];
C = C(:,1:num_rows);
C_12 = C(:,1:12);
D = zeros(size(C,1),size(B,2));
%Some data smoothing to improve performance
flux_trajectory_arr_smooth = smoothdata(flux_trajectory_arr);
%% From Rowley De-Biasing DMD paper
%{ Step 1) Stack input and output as Z. Z = [X;Y]. Where X is data from 0
%to end -1, and Y is data from 1 to end

start = 700; %starting time sample for training
stop = 1300; %ending time sample for training
X = flux_trajectory_arr_smooth(:,start:stop-1);
Y = flux_trajectory_arr_smooth(:,start+1:stop);
% We will now subtract off the actuation from the X data set. 
X_sub = X; % This is our data with the actuation subtracted off

%We are now going to stack our X_sub and Y matrices
Z = [X_sub;Y];

%Now we take the SVD of Z and construct our projection operator
[U,S,V] = svd(Z,'econ');

%now we are going to retain the first r cols of V, where r is our
%truncation rank.
r = 9;
V_r = V(:,1:r);
%now our projection operator
P = V_r*V_r';

%Now we project our data matrices onto this operator
X_sub_bar = X*P;
Y_bar = Y*P;

%Now we calculate another SVD, this time on our projected input data
[U_bar,S_bar,V_bar] = svd(X_sub_bar,'econ');

A_tilde = U_bar'*Y_bar*V_bar*pinv(S_bar);
%% Let's do some simulating now

eigs = eig(A_tilde);

scatter(real(eigs),imag(eigs))

sys_d_dmd = ss(A_tilde,U_bar'*B,C,D,dT);


[~,~,xout_dmd] = lsim(sys_d_dmd,inputs,time);
xout_dmd = xout_dmd';
%%
figure()
for k = 1:length(A_tilde)-1
    subplot(4,3,k);
    plot(time,xout_dmd(k,:))
    hold on
    plot(time,flux_trajectory_arr(k,:))
    legend('OPT-DMD','Test')
    
end
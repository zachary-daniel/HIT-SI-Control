clear;close all;clc;

%File for training BOP-DMD model on one shot, and then seeing how well that
%model extrapolates to the shot immediately before it. In this case we are
%going to train on 221201009, and then test on 221201010. Or vice versa.
%Doesn't really matter...
%221201009 will be referred to as 'shot 11'
shot_221201009 = load('shots\high_itor_plasma_shots\matlab data\221201009_plasma_flux.mat');
%shot 221201010 will be referred to as 'shot 12'
shot_221201010 = load('shots\high_itor_plasma_shots\matlab data\221201010_plasma_flux.mat');

%array of shot 11 circuit measurments
flux_trajectory_arr_9 = [shot_221201009.i_L1_1 , shot_221201009.v_cap_1, shot_221201009.i_fcoil_1,...
                      shot_221201009.i_L1_2 , shot_221201009.v_cap_2 , shot_221201009.i_fcoil_2,...
                      shot_221201009.i_L1_3 , shot_221201009.v_cap_3 , shot_221201009.i_fcoil_3,...
                      shot_221201009.i_L1_4 , shot_221201009.v_cap_4 , shot_221201009.i_fcoil_4,...
                      shot_221201009.i_tor];



%array of shot 12 circuit measurments
flux_trajectory_arr_10 = [shot_221201010.i_L1_1 , shot_221201010.v_cap_1, shot_221201010.i_fcoil_1,...
                      shot_221201010.i_L1_2 , shot_221201010.v_cap_2 , shot_221201010.i_fcoil_2,...
                      shot_221201010.i_L1_3 , shot_221201010.v_cap_3 , shot_221201010.i_fcoil_3,...
                      shot_221201010.i_L1_4 , shot_221201010.v_cap_4 , shot_221201010.i_fcoil_4,...
                      shot_221201010.i_tor];


time = shot_221201009.time'; %manipulate array shape to conform to standard
[~,t_equal_0] = min(abs(time)); %find when time = 0
time = time(:,t_equal_0:end); %eliminate all samples before this time
flux_trajectory_arr_9 = flux_trajectory_arr_9'; %flip arrays to conform to standard
flux_trajectory_arr_9 = flux_trajectory_arr_9(:,t_equal_0:end); %start array at t = 0

flux_trajectory_arr_10 = flux_trajectory_arr_10'; %flip array to conform to standard
flux_trajectory_arr_10 = flux_trajectory_arr_10(:,t_equal_0:end); %start array at t = 0

%Power supply waveforms for shot 11
inputs_9 = [shot_221201009.v_spa_1 , shot_221201009.v_spa_2 , ...
    shot_221201009.v_spa_3,shot_221201009.v_spa_4]';

%power supply waveforms for shot 12
inputs_10 = [shot_221201010.v_spa_1 , shot_221201010.v_spa_2 , ...
    shot_221201010.v_spa_3,shot_221201010.v_spa_4]';

%Set time series to start at 0
inputs_9 = inputs_9(:,t_equal_0:end);
inputs_10 = inputs_10(:,t_equal_0:end);

%Some of these shots have a dc offset as an artifact so we're going to
%remove that before doing any sort of training or testing.

% mean_9 = mean(flux_trajectory_arr_9');
% flux_trajectory_arr_9(1:12,:) = flux_trajectory_arr_9(1:12,:) - mean_9(1:12)';
% mean_10 = mean(flux_trajectory_arr_10');
% flux_trajectory_arr_10(1:12,:) = flux_trajectory_arr_10(1:12,:) - mean_10(1:12)';



%% This section is for training OPT-DMD on shot 11
train_len = length(time)-462; %This is rougly where the power supplies turn off 
r = 5; %r = 5, and r = 7 seem to give the best data. 
imode = 2;
train= flux_trajectory_arr_10(1:13,train_len:end);




shape = size(train(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];
train_time = time(train_len:end);
copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w_init,e_init,b,converged,atilde,u,afull] = optdmd(train,train_time,r,imode,[],[],[],copts);
figure()
scatter(real(e_init),imag(e_init))
xlabel('real')
ylabel('imag')
title('opt-dmd eigs')




%%


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

B_12 = B(1:12,:);

C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;];

C_12 = C(:,1:12);
D = zeros(size(C,1),size(B,2));

B_tilde = B(1:length(afull),:);

if shape(1) == 12
    sys_opt = ss(real(afull),B_tilde,C_10,D);
else
    sys_opt = ss(real(afull),B_tilde,C(:,1:length(afull)),D);
end
[y,t,xout_opt] = lsim(sys_opt,inputs_9',time);
xout_opt = xout_opt';
%% BOP-DMD training
bop_train = train; %matrix of training data
bop_time = train_time; %training time vector
num_trials = 20; %number of trials to run
percent = .8;
[w_avg_no_conj,e_avg_no_conj,b_avg_no_conj,atilde_bop_no_conj] = bop_dmd_func(bop_train,bop_time,r,percent,num_trials,[],[],[e_init],[]);

%BOP with complex conjugate enforced
[w_avg,e_avg,b_avg,atilde_bop] = bop_dmd_func(bop_train,bop_time,r,percent,num_trials,[],[],[e_init],true);

sys_bop_no_conj = ss(atilde_bop_no_conj,B_tilde(1:shape,:),C(:,1:shape),D); %bop with no complex conjugate eforcement

[~,~,xout_bop_no_conj] = lsim(sys_bop_no_conj,inputs_9',time); %simulate bop with no conjugte

sys_bop = ss(atilde_bop,B_tilde(1:shape,:),C(:,1:shape),D); %bop with no complex conjugate eforcement

[~,~,xout_bop] = lsim(sys_bop,inputs_9',time); %simulate bop with no conjugte

xout_bop_no_conj = xout_bop_no_conj'; %fix dimensions to match convention

xout_bop = xout_bop';


%% Plot results compared to test data
figure()
scatter(real(e_init),imag(e_init))
hold on
scatter(real(e_avg_no_conj),imag(e_avg_no_conj))
xlabel('real')
ylabel('imag')
legend('opt eigs','bop eigs')
title('BOP-dmd eigs')

%%

figure() %OPT vs. test data
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time,xout_opt(k,:))
    hold on
    plot(time,flux_trajectory_arr_9(k,:))
    legend('OPT-DMD','Test')
    
end
%%
figure() %BOP no conj vs. test
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time,xout_bop_no_conj(k,:))
    hold on
    plot(time,flux_trajectory_arr_9(k,:))
    legend('BOP-DMD no conj','Test')
end

figure() %BOP with conjugate vs. test
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time,xout_bop(k,:))
    hold on
    plot(time,flux_trajectory_arr_9(k,:))
    legend('BOP-DMD','Test')
    
end
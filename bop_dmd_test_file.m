clear;close all;clc;
load("vacuum_model_A_matrix.mat")

full_eigs = eigs(A);


%%

load('Python_Stuff\data\plasma_shots\flux_data\220816009\220816009_flux.mat')
%Two good shots are 220816009,221129011
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')

%%

flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,...
                      i_L1_2,v_cap_2,i_fcoil_2,...
                      i_L1_3,v_cap_3,i_fcoil_3,...
                      i_L1_4,v_cap_4,i_fcoil_4,i_tor];



time = time';
[~,t_equal_0] = min(abs(time));
time = time(:,t_equal_0:end);
flux_trajectory_arr = flux_trajectory_arr';
flux_trajectory_arr = flux_trajectory_arr(:,t_equal_0:end);

inputs = [v_spa_1,v_spa_2,v_spa_3,v_spa_4]';
inputs = inputs(:,t_equal_0:end);
train_len = length(time)-483; %This is rougly where the power supplies turn off 

%%
r = 9;
[U,S,V] = svd(A,'econ');
S(r+1:end,r+1:end) = 0;
A_rom =  U*S*V';
r_eigs = eigs(A_rom);

%%

imode = 1;
train= flux_trajectory_arr(1:end,train_len:end);


e_init_vacuum_aug = [r_eigs;r_eigs(1)*ones(r-length(r_eigs),1)];

shape = size(train(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];
train_time = time(train_len:end);
copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w,e_init,b,converged,atilde,u,afull] = optdmd(train,train_time,r,imode,[],[],[],copts);
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
    sys_opt = ss(real(afull),B_tilde,C_12,D);
else
    sys_opt = ss(real(afull),B_tilde,C(:,1:length(afull)),D);
end
[y,t,xout_opt] = lsim(sys_opt,inputs',time);
xout_opt = xout_opt';


%% Plot result of opt-DMD fit compared to test data
figure()
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time(),xout_opt(k,:))
    hold on
    plot(time(),flux_trajectory_arr(k,:))
    legend('OPT-DMD','Test')
    
end
%% See how BOP-DMD fairs
bop_train = train(:,1:length(train)-0);
bop_time = train_time(1:length(train_time)-0);
num_trials = 20;
[w_avg,e_avg,b_avg,atilde_bop] = bop_dmd_func(bop_train,bop_time,r,.9,num_trials,[],[],[e_init],true);
[w_avg_no_conj,e_avg_no_conj,b_avg_no_conj,atilde_bop_no_conj] = bop_dmd_func(bop_train,bop_time,r,.75,num_trials,[],[],[e_init],[]);

%%

figure()

scatter(real(e_avg),imag(e_avg))
xlabel('real')
ylabel('imag')
title('bop-dmd eigs')
hold on
scatter(real(e_init),imag(e_init))
legend('BOP','OPT')





sys_bop = ss(atilde_bop,B,C,D);
sys_bop_no_conj = ss(atilde_bop_no_conj,B,C,D);

[~,~,xout_bop] = lsim(sys_bop,inputs',time);
xout_bop = xout_bop';

[~,~,xout_bop_no_conj] = lsim(sys_bop_no_conj,inputs',time);
xout_bop_no_conj = xout_bop_no_conj';

%Plot Bop results against ground truth
figure()
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time,xout_bop(k,:))
    hold on
    plot(time,flux_trajectory_arr(k,:))
    legend('BOP-DMD','Test')
    
end

figure()
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time(),xout_bop_no_conj(k,:))
    hold on
    plot(time(),flux_trajectory_arr(k,:))
    legend('BOP-DMD no conj','Test')
    
end
%This will save the current A matrix to the given file. Uncomment to update
%this
% save('Python_Stuff/data/bop_dmd_plasma_shot_220816009.mat',"xout_bop_no_conj");
% save("Python_Stuff/data/bop_dmd_A_matrix_220816009.mat","atilde_bop_no_conj");
%Error calculations
BOP_error = norm(sqrt( (xout_bop(1:12,:)).^2 - flux_trajectory_arr(1:12,:).^2 ));
OPT_error = norm(sqrt( (xout_opt(1:12,:)).^2 - flux_trajectory_arr(1:12,:).^2 ));
BOP_no_conj_error = norm(sqrt( (xout_bop_no_conj(1:12,:)).^2 - flux_trajectory_arr(1:12,:).^2 ));
disp('BOP-DMD LSTSQ error:')
disp(BOP_error)
disp('OPT-DMD LSTSQ error: ')
disp(OPT_error)
disp('BOP-DMD no conj LSTSQ error:')
disp(BOP_no_conj_error)


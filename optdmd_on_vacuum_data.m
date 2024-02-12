clear;close all;clc;
load('Python_Stuff/data/vacuum_shots_220816/good_shots/220816005.mat')
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')



flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,i_L1_2,v_cap_2,i_fcoil_2,i_L1_3,v_cap_3,i_fcoil_3,i_L1_4,v_cap_4,i_fcoil_4];

time = time';
dT = time(2) - time(1);
flux_trajectory_arr = flux_trajectory_arr';


%%
r = 9;

imode = 2;
train = flux_trajectory_arr(1:12,1773:end);
shape = size(train(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
train_time = time(1773:end);
[w,e_init,b,converged,atilde,u,afull] = optdmd(train,train_time,r,imode,[],[],[],copts);
num_trials = 10;
[w_avg_no_conj,e_avg_no_conj,b_avg_no_conj,atilde_bop_no_conj] = bop_dmd_func(train,train_time,r,.8,num_trials,[],[],[e_init],[]);
%%
inputs = [v_spa_1,v_spa_2,v_spa_3,v_spa_4];

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



if shape(1) == 12
    sys_opt = ss(real(afull),B_12,C_12,D);
else
    sys_opt = ss(real(afull),B_12,C,D);
end



% [y,t,xout_plasma] = lsim(sys_plasma,inputs(252:end,:),time(252:end));
% r= 12
% [w,e,b,atilde,u,afull] = optdmd(train,time(1773:end),r,imode);

[~,~,xout_opt] = lsim(sys_opt,inputs(252:end,:),time(252:end));
xout_opt = xout_opt';
sys_bop = ss(atilde_bop_no_conj,B_12,C_12,D);
[~,~,xout_bop] = lsim(sys_bop,inputs(252:end,:),time(252:end));
xout_bop = xout_bop';

%% Plot result of opt-DMD fit compared to test data
figure()
for k = 1:12
    subplot(4,3,k);
    plot(time(252:end),xout_opt(k,:),'r','LineWidth',2)
    hold on
    plot(time(252:end),flux_trajectory_arr(k,252:end),'-k')
    legend('OPT-DMD','Test')
    
end

figure()
for k = 1:12
    subplot(4,3,k);
    plot(time(252:end),xout_bop(k,:),'r','LineWidth',2)
    hold on
    plot(time(252:end),flux_trajectory_arr(k,252:end),'-k')
    legend('BOP-DMD','Test')
    
end

figure()

scatter(real(e_avg_no_conj),imag(e_avg_no_conj))
xlabel('real')
ylabel('imag')
title('bop-dmd eigs')
hold on
scatter(real(e_init),imag(e_init))
legend('BOP','OPT')






BOP_error = norm(sqrt( (xout_bop(1:12,:)).^2 - flux_trajectory_arr(1:12,252:end).^2 ));
OPT_error = norm(sqrt( (xout_opt(1:12,:)).^2 - flux_trajectory_arr(1:12,252:end).^2 ));

disp('BOP-DMD LSTSQ error:')
disp(BOP_error)
disp('OPT-DMD LSTSQ error: ')
disp(OPT_error)
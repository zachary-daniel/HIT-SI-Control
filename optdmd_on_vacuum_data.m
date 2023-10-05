clear;close all;clc;
load('Python_Stuff/data/vacuum_shots_220816/good_shots/220816005.mat')
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')



flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,i_L1_2,v_cap_2,i_fcoil_2,i_L1_3,v_cap_3,i_fcoil_3,i_L1_4,v_cap_4,i_fcoil_4];

time = time';
flux_trajectory_arr = flux_trajectory_arr';

%%
r = 5;
imode = 2;
train = flux_trajectory_arr(1:12,1773:end);
shape = size(train(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w,e,b,converged,atilde,u,afull] = optdmd(train,time(1773:end),r,imode,[],[],[],copts);
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

% [y,t,xout_plasma] = lsim(sys_plasma,inputs(252:end,:),time(252:end));
% r= 12
% [w,e,b,atilde,u,afull] = optdmd(train,time(1773:end),r,imode);
if shape(1) == 12
    sys_opt = ss(real(afull),B_12,C_12,D);
else
    sys_opt = ss(real(afull),B,C,D);
end
[y,t,xout_test] = lsim(sys_opt,inputs(252:end,:),time(252:end));
%% Plot result of opt-DMD fit compared to test data
figure()
for k = 1:12
    subplot(4,3,k);
    plot(time(252:end),xout_test(:,k),'r','LineWidth',2)
    hold on
    plot(time(252:end),flux_trajectory_arr(k,252:end),'-k')
    legend('OPT-DMD','Test')
    
end
%%

[w_avg,e_avg,b_avg,atilde_avg,flag] = bop_dmd_func_avg(10,train(:,1:270),.6,5,imode,time(1773:end-114),30,1e-10);

figure()
scatter(real(e_avg),imag(e_avg));
hold on
scatter(real(e),imag(e));
legend('BOP','OPT')
%%
sys_bop_avg = ss(atilde_avg,B_12,C_12,D);
[~,~,xout_bop_avg] = lsim(sys_bop_avg,inputs(252:end,:),time(252:end));


figure()
for k = 1:12
    subplot(4,3,k);
    plot(time(252:end),xout_bop_avg(:,k))
    hold on
    plot(time(252:end),flux_trajectory_arr(k,252:end))
    legend('BOP-DMD','Test')
    
end


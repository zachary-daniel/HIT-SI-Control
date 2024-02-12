clear;close all;clc;
load('Python_Stuff/data/plasma_shots/flux_data/220816009/220816009_flux.mat')
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')



flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,...
    i_L1_2,v_cap_2,i_fcoil_2,i_L1_3,v_cap_3,...
    i_fcoil_3,i_L1_4,v_cap_4,i_fcoil_4,i_tor];

time = time';
flux_trajectory_arr = flux_trajectory_arr';

%%
r = 9;
imode = 2;
train = flux_trajectory_arr(1:13,1773:end);
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
    plot(time(252:end),xout_test(:,k))
    hold on
    plot(time(252:end),flux_trajectory_arr(k,252:end))
    legend('OPT-DMD','Test')
    
end


%% Let's try and see if this bagging wrapper I wrote for the opt-DMD function helps at all

% [w_cluster,e_cluster,b_cluster,atilde_bop,flag] = bop_dmd_func_cluster(1,train,1,4,imode,time(1773:end),40);
% 
% afull_bop = w_cluster*diag(e_cluster)*pinv(w_cluster);
% 
% % w_avg = 0;
% % e_avg = 0;
% % b_avg = 0;
% % for k = 1:4
% %     [w,e,b,atilde,u,afull_bop] = optdmd(train,time(1773:end),r,imode);
% %     w_avg = w_avg + w;
% %     b_avg = b_avg + b;
% %     e_avg = e_avg + e;
% % end
% % w_avg = w_avg/k;
% % e_avg = e_avg/k;
% % b_avg = b_avg/k;
% % 
% % afull_bop = w_avg*diag(e_avg)*pinv(w_avg);
% if shape(1) == 12
%     sys_bop = ss(atilde_bop,B_12,C_12,D);
% else
%     sys_bop = ss(atilde_bop,B,C,D);
% end
%%
% [t,y,xout_bop] = lsim(sys_bop,inputs(252:end,:),time(252:end));
% 
% % plot(time(252:end),recon')
% figure()
% for k = 1:12
%     subplot(4,3,k);
%     plot(time(252:end),xout_bop(:,k))
%     hold on
%     plot(time(252:end),flux_trajectory_arr(k,252:end))
%     legend('OPT-DMD','Test')
%     
% end

% figure()
% for k = 1:12
%     subplot(4,3,k);
%     plot(time(252:end),xout_bop(:,k))
%     hold on
%     plot(time(252:end),xout_test(:,k))
%     plot(time(252:end),flux_trajectory_arr(k,252:end))
%     legend('BOP-DMD','Opt-DMD','Test')
%     
% end

%% Let's try BOP-DMD with some regular averaging now instead of the clustering
% train_norm = train./(max(train'))';
% 
% 
% %%
% [w_avg,e_avg,b_avg,atilde_avg,flag] = bop_dmd_func_avg(20,train,.63,4,imode,time(1773:end),40);
% 
% figure()
% scatter(real(e_avg),imag(e_avg));
% hold on
% scatter(real(e),imag(e));
% legend('BOP','OPT')
% %%
% sys_bop_avg = ss(atilde_avg,B,C,D);
% [~,~,xout_bop_avg] = lsim(sys_bop_avg,inputs(252:end,:),time(252:end));
% 
% 
% figure()
% for k = 1:12
%     subplot(4,3,k);
%     plot(time(252:end),xout_bop_avg(:,k))
%     hold on
%     plot(time(252:end),flux_trajectory_arr(k,252:end))
%     legend('BOP-DMD','Test')
%     
% end
% 
% %% Fuck around with normalized data
% % [w,e,b,converged,atilde,u,afull] = optdmd(train_norm,time(1773:end),r,imode,[],[],[],copts);
% % 
% % sys_bop_norm = ss(afull,B,C,D);
% % [~,~,xout_bop_norm] = lsim(sys_bop_avg,inputs(1773:end,:),time(1773:end));
% % 
% % figure()
% % for k = 1:12
% %     subplot(4,3,k);
% %     plot(time(1773:end),xout_bop_norm(:,k))
% %     hold on
% %     plot(time(1773:end),flux_trajectory_arr(k,1773:end))
% %     legend('OPT-Norm-Data','Test')
% %     
% % end
% %% I'm going to take this model with the slightly different B matrix and try it on a new plasma shot
% load('shots\220816008_flux.mat');
% 
% %% Load in second shot
% flux_trajectory_arr_2= [i_L1_1,v_cap_1,i_fcoil_1,i_L1_2,v_cap_2,i_fcoil_2,i_L1_3,v_cap_3,i_fcoil_3,i_L1_4,v_cap_4,i_fcoil_4,i_tor];
% inputs_2 = [v_spa_1,v_spa_2,v_spa_3,v_spa_4];
% time = time';
% flux_trajectory_arr_2 = flux_trajectory_arr_2';
% 
% train = flux_trajectory_arr_2(1:13,1773:end);
% %% I'm going to take the BOP model from one shot and apply it to another shot. I'm also going to see if this changed value 
% % of L1 which goes into the B matrix improves performance across multiple
% % shots
% 
% 
% figure()
% for k = 1:12
%     subplot(4,3,k);
%     plot(time(252:end),flux_trajectory_arr(k,252:end))
%     hold on
%     plot(time(252:end),flux_trajectory_arr_2(k,252:end))
%     legend('220816009','220816008')
%     
% end
% 
% figure()
% for k = 1:4
%     subplot(4,1,k)
%     plot(time(252:end),inputs(252:end,k))
%     hold on
%     plot(time(252:end),inputs_2(252:end,k))
%     legend('220816009','220816008')
% end
% %% OK so this B-matrix changing around is valid. Makes performance better across different shots from same day
% [~,~,xout_bop_new_shot] = lsim(sys_bop_avg,inputs_2(252:end,:),time(252:end));
% figure()
% for k = 1:12
%     subplot(4,3,k);
%     plot(time(252:end),xout_bop_new_shot(:,k))
%     hold on
%     plot(time(252:end),flux_trajectory_arr_2(k,252:end))
%     legend('BOP-DMD','New Shot')
%     
% end
% %% Compare BOP-DMD to training data
% % [~,~,xout_train_fit] = lsim(sys_bop_avg,inputs(1773:end,:),time(1773:end));
% % 
% % figure()
% % for k = 1:12
% %     subplot(4,3,k);
% %     plot(time(1773:end),xout_train_fit(:,k))
% %     hold on
% %     plot(time(1773:end),flux_trajectory_arr(k,1773:end))
% %     legend('BOP-DMD','Training Data')
% %     
% % end
% 
% %% Maybe my B matrix is wrong for when a plasma is present. Let's compute one analytically from the data
% % t_vec = time.*ones(13,length(time));
% beh = ( train + 1 - expm(sys_bop_avg.A.*time(:,1773:end)) - expm(sys_bop_avg.A.*t(:,1773:end))*train(:,1) )/inputs(1773,1);
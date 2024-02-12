clear;close all;clc;
load('Python_Stuff\data\plasma_shots\flux_data\220816009\220816009_flux.mat')
%Two good shots are 220816009,221129011
% load('Python_Stuff/data/plasma_ss_A_matrix.mat')

%%

flux_trajectory_arr = [i_L1_1,v_cap_1,i_fcoil_1,...
                      i_L1_2,v_cap_2,i_fcoil_2,...
                      i_L1_3,v_cap_3,i_fcoil_3,...
                      i_L1_4,v_cap_4,i_fcoil_4,i_tor];
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
flux_trajectory_arr = flux_trajectory_arr';
flux_trajectory_arr = flux_trajectory_arr(:,t_equal_0:end);

inputs = [v_spa_1,v_spa_2,v_spa_3,v_spa_4]';
inputs = inputs(:,t_equal_0:end);
train_len = length(time)-483; %This is rougly where the power supplies turn off 

%%
r = 9;
imode = 2;
train_rough = flux_trajectory_arr(1:end,train_len:end);
%%


train = (train_rough);
shape = size(train(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];
train_time = time(train_len:end);
copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w,e_init,b,converged,atilde,u,afull] = optdmd(train,train_time,r,imode,[],[],[],copts);
% e_init(6) = real(e_init(6));
% e_init(7) = real(e_init(7));
% w(:,7) = real(w(:,7));
% afull = w*(diag(e_init))*pinv(w);
figure()
scatter(real(e_init),imag(e_init))
xlabel('real')
ylabel('imag')
title('opt-dmd eigs')
% [w,e,b,converged,atilde,u,afull] = optdmd(train,train_time,r,imode,[],e_init,[],copts);

%Take the SVD of the training data. We are going to use the first r rows of
%the U matrix to modify our B matrix and get a good reduced order actuation
%matrix
[U,S,V] = svd(flux_trajectory_arr);

U(r+1:end,:) = 0;

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
[y,t,xout_test] = lsim(sys_opt,inputs',time);
xout_test = xout_test';




%% Plot result of opt-DMD fit compared to test data
figure()
for k = 1:length(afull)-1
    subplot(4,3,k);
    plot(time,xout_test(k,:))
    hold on
    plot(time,flux_trajectory_arr(k,:))
    legend('OPT-DMD','Test')
    
end



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
train_norm = train./(max(train'))';


%%
% [w_avg,e_avg,b_avg,atilde_avg,flag] = bop_dmd_func_avg(10,train(:,1:end),.6,r,imode,train_time(1:end),40,1e-10,true,e_init);
% % for k = 3:2:7
% %     imag_avg = .5*(abs(imag(e_avg(k))) + abs(imag(e_avg(k+1))));
% %     real_avg = .5*(real(e_avg(k)) + real(e_avg(k+1)));
% %     e_avg(k) = real_avg + 1i*imag_avg;
% %     e_avg(k+1) = real_avg - 1i*imag_avg;
% % %     real_avg_w = .5*(real(w_avg(:,k)) + real(w_avg(:,k+1)));
% % %     imag_avg_w = .5*(abs(imag(w_avg(:,k))) + abs(imag(w_avg(:,k+1))));
% % %     w_avg(:,k) = real_avg_w -1i*imag_avg_w;
% % %     w_avg(:,k+1) = real_avg_w +1i*imag_avg_w;
% % end
% atilde_avg = w_avg*diag(e_avg)*pinv(w_avg);
% 
% sys_bop = ss(atilde_avg,B_tilde,C(:,1:length(afull)),D);
% %%
% figure()
% scatter(real(e_avg),imag(e_avg));
% hold on
% scatter(real(e_init),imag(e_init));
% legend('BOP','OPT')
% 
% afull_vecnorm = vecnorm(afull);
% atilde_avg_vecnorm = vecnorm(atilde_avg);
% ratio = atilde_avg_vecnorm./afull_vecnorm;
% 
% w_opt = vecnorm(w*pinv(w));
% w_bop = vecnorm(w_avg*pinv(w_avg));
% ratio_w = w_bop./w_opt;
% %%
% 
% [~,~,xout_bop] = lsim(sys_bop,inputs,time);
% xout_bop = xout_bop';
% 
% figure()
% for k = 1:length(afull)-1
%     subplot(4,3,k);
%     plot(time,xout_bop(k,:))
%     hold on
%     plot(time,flux_trajectory_arr(k,:))
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

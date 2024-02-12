clear;close all;clc;
% load('Python_Stuff/data/plasma_shots/flux_data/220816009/220816009_flux.mat')

% Here we load in both the vacuum and plasma dynamics matrices
A = load('Python_Stuff/data/plasma_ss_A_matrix.mat').A;
A_vacuum = load('vacuum_model_A_matrix.mat').A';
load('Python_Stuff\data\voltages.mat');
e = eigs(A);
e_vacuum = eigs(A_vacuum);

%%

% time = linspace(0,.004,length(voltage_from_python));


%Make the B, C, and D matrix for both the vacuum and plasma models

L1 = (8.0141e-7); %Henry
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

B_12 = B(1:12,:); %for vacuum

C = [0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0;];

C_12 = C(:,1:12); %for vacuum
D = zeros(size(C,1),size(B,2));

sys_vacuum = ss(A_vacuum,B_12,C_12,D); %vacuum model continuous time system
sys_plasma = ss(A,B,C,D); %plasma model continuous time system
voltages = [newVoltage,newVoltage,newVoltage,newVoltage]; %Stack the voltages up. This is the order that's used in the corresponding python notebook
%%                                   %though it shouldn't matter
voltages(30001:end,:) = 0; %Set the last millisecond of the voltages  = 0 which is the same time on the experiment when the SPA's turn off

[t,y,xout_plasma] = lsim(sys_plasma,voltages,time); %simulate our plasma system

[t,y,xout_vacuum] = lsim(sys_vacuum,voltages,time); %simulate acuum system

figure()
title('Plasma Simulated Data')
%plot plasma system
for k = 1:12
    
    subplot(4,3,k)
    plot(xout_plasma(:,k))
end

figure()
title('Vacuum Simulated Data')
%plot vacuum system
for k = 1:12
    
    subplot(4,3,k)
    plot(xout_vacuum(:,k))
end
%% OPT-DMD
r = 12; %number of modes
imode = 2;
train = xout_vacuum(30001:end,:)'; %get training data. Has to be states X snapshots
shape = size(train(:,1));

train_no_control = (xout_vacuum - (B_12*voltages')')';

lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w,e,b,converged,atilde,u,afull] = optdmd(train,time(30001:end),r,imode,[],[],[],copts); %Fit to vacuum sim data
%%
[~,~,~,converged,~,~,afull_no_control] = optdmd(train_no_control(:,20000:30000),time(20000:30000),r,imode,[],[],[],copts); %Fit to vacuum sim data with control subtracted off

sys_no_control = ss(afull_no_control,B_12,C_12,D);
%%


% [y,t,xout_plasma] = lsim(sys_plasma,inputs(252:end,:),time(252:end));
% r= 12
% [w,e,b,atilde,u,afull] = optdmd(train,time(1773:end),r,imode);
if shape(1) == 12
    sys_opt = ss(afull,B_12,C_12,D);
else
    sys_opt = ss(afull,B,C,D);
end
[y,t,xout_test] = lsim(sys_opt,voltages,time); %get test data. Just simulate the whole shot that it was trained on. 
[~,~,xout_no_control] = lsim(sys_no_control,voltages,time);

%%
figure() %plot DMD model vs. test data
title('OPT-DMD fit of clean simulated vacuum data')
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_test(:,k),"LineWidth", 2)
    hold on
    plot(time,xout_vacuum(:,k))
    legend('DMD','Vacuum Test')
end
%%
figure() %plot DMD model vs. test data
title('OPT-DMD trained on data with control subtracted off')
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_no_control(:,k),"LineWidth", 2)
    hold on
    plot(time,xout_vacuum(:,k))
    legend('DMD no control','Vacuum Test')
end
%%
opt_vacuum_eigs = eigs(sys_opt.A);
opt_no_control_eigs = eigs(sys_no_control.A);
figure()

scatter(real(e_vacuum),imag(e_vacuum),100)
hold on
scatter(real(opt_vacuum_eigs),imag(opt_vacuum_eigs),'d')
hold on
scatter(real(opt_no_control_eigs),imag(opt_no_control_eigs));
legend('Model','DMD','Control subtracted off')
%% Now we do the same for the plasma data
r = 13; %number of modes. Full rank in this case
imode = 2;
train_plasma = xout_plasma(30001:end,:)'; %get training data. Has to be states X snapshots
shape = size(train_plasma(:,1));
lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
ubc = [zeros(r,1); Inf*ones(r,1)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
[w,e,b,converged,atilde,u,afull_plasma] = optdmd(train_plasma,time(30001:end),r,imode,[],[],[],copts); %Fit to plasma sim data

%%

if shape(1) == 12
    sys_opt_plasma = ss(afull_plasma,B_12,C_12,D);
else
    sys_opt_plasma = ss(afull_plasma,B,C,D);
end
[y,t,xout_test] = lsim(sys_opt_plasma,voltages,time); %get test data. Just simulate the whole shot that it was trained on. 

%%
figure() %plot DMD model vs. test data for plasma 
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_test(:,k),"LineWidth", 2)
    hold on
    plot(time,xout_plasma(:,k))
    legend('DMD','Plasma Test')
end

%% Comparrison of eigenvalues from OPT-DMD and true model
opt_plasma_eigs = eigs(sys_opt_plasma.A);

figure()
scatter(real(e),imag(e),100)
hold on
scatter(real(opt_plasma_eigs),imag(opt_plasma_eigs),'d')

legend('Model eigs','DMD')

%% Simualte opt plasma model with new voltage waveforms

new_input = [newVoltage,newVoltageShift2,newVoltageShift1,newVoltageShift3];
new_input(40000:end) = 0;
[t,y,xout_test] = lsim(sys_opt_plasma,voltages,time);

% [t,y,xout_plasma] = lsim(sys_plasma,voltages,time);

%%
figure() %plot DMD model vs. test data for plasma 
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_test(:,k),"LineWidth", 2)
    hold on
    plot(time,xout_plasma(:,k))
    legend('DMD','Plasma Test')
end
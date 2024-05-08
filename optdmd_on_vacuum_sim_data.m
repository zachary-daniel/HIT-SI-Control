clear;close all;clc;
% load('Python_Stuff/data/plasma_shots/flux_data/220816009/220816009_flux.mat')

% Here we load in both the vacuum and plasma dynamics matrices
A = load('Python_Stuff/data/plasma_ss_A_matrix.mat').A;
A_vacuum = load('vacuum_model_A_matrix.mat').A';
load('Python_Stuff\data\voltages.mat');
e_plasma = eigs(A);
e_vacuum = eigs(A_vacuum);
SampleTime = 1e-7;
time = 0:SampleTime:.004;

voltage_inputs = [newVoltage,newVoltageShift1,newVoltageShift2,newVoltageShift3];
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

%%                                   %though it shouldn't matter
voltage_inputs(:,:) = 0; %Set the last millisecond of the voltage_inputs  = 0 which is the same time on the experiment when the SPA's turn off
voltage_inputs(1,1) = 100; %Set first value of input array to one to simulate an impulse response
[t,y,xout_plasma] = lsim(sys_plasma,voltage_inputs,time); %simulate our plasma system

[t,y,xout_vacuum] = lsim(sys_vacuum,voltage_inputs,time); %simulate acuum system

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
train = xout_vacuum(1:10000,:)'; %get training data. Has to be states X snapshots
shape = size(train(:,1));


lbc = [-Inf*ones(r,1); -Inf*ones(r,1)]; %Stability constraints
ubc = [zeros(r,1); Inf*ones(r,1)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

[w,e,b,converged,atilde,u,afull_vacuum] = optdmd(train,time(1,1:10000),r,imode,[],[],[],copts); %Fit to vacuum sim data
afull_vacuum = w*diag(e)*pinv(w);
%%

sys_opt = ss(afull_vacuum,B_12,C_12,D);

[y,t,xout_dmd] = lsim(sys_opt,voltage_inputs,time); %get test data. Just simulate the whole shot that it was trained on. 


%%
figure() %plot DMD model vs. test data
title('OPT-DMD fit of clean simulated vacuum data')
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_dmd(:,k),"LineWidth", 2)
    hold on
    plot(time,xout_vacuum(:,k))
    legend('DMD','Vacuum Test')
end

figure() %plot DMD model vs. test data
title('OPT-DMD fit of clean simulated vacuum data')
for k = 1:12
    subplot(4,3,k);
    plot(time,xout_dmd(:,k),"LineWidth", 2)

end


%% Plot the residual of the model as a function of time

% figure() %plot DMD model vs. test data
% title('OPT-DMD trained on data with control subtracted off')
% for k = 1:12
%     subplot(4,3,k);
%     plot(time,xout_no_control(:,k),"LineWidth", 2)
%     hold on
%     plot(time,xout_vacuum(:,k))
%     legend('DMD no control','Vacuum Test')
% end

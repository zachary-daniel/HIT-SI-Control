clear;close all;clc;
load("files_for_dmd.mat")
%load("DMDc_on_old_shot.mat")

% sys_dmd_on_shot_c = sys_dmd;
% Let's first check that DMDc can somewhat accurately recreate the A and B
% matrices that come out of the state space model I derived
dT = 1e-7;
time = 0:dT:.004;
B = sysc.B;
Q = .001; %diag(.001*ones(1,size(A, 1))); % disturbance covariance
R = diag(10*ones(1,size(B,2))); % Noise covariance
noisePower = 0;
theta1 = 0;
theta2 = 0;
theta3 = 0;
newVoltageShift1 = phaseShift(newVoltage,theta1);
newVoltageShift2 = phaseShift(newVoltage,theta2);
newVoltageShift3 = phaseShift(newVoltage,theta3);
input_voltages = [newVoltage,newVoltageShift1,newVoltageShift2,newVoltageShift3];
[y,tout,xout] = lsim(sysc,input_voltages,time);

%xout is full state measurment

%get the shapes of everything all set up
xout = xout';
input_voltages = input_voltages';
%stack inputs and states
omega = [xout(:,1:end-1);input_voltages(:,1:end-1)];
%get inputs from 2nd position to the end
x2 = xout(:,2:end);
x = xout(:,1:end-1);
%how many singular values we are going to use
trunc = 12;

%take an svd of our stacked inputs and states
[utilde,stilde,vtilde] = svd(omega,'econ');

%truncate all of our resulting matrices accordingly
utilde_trunc = utilde(:,1:trunc);

stilde_trunc = stilde(1:trunc,1:trunc);

vtilde_trunc = vtilde(:, 1:trunc);

utilde_trunc_1 = utilde_trunc(1:12,:);
utilde_trunc_2 = utilde_trunc(13:16,:);

Atilde = x2*vtilde_trunc*pinv(stilde_trunc)*utilde_trunc_1';
Btilde = x2*vtilde_trunc*pinv(stilde_trunc)*utilde_trunc_2';

sysd = c2d(sysc,dT);
%This does a poor job, but maybe if we do foreward backwards dmd we will
%get a better result?

sys_dmdc_d = ss(Atilde,Btilde,sysd.C,sysd.D,dT);

[y_dmdc,tout_dmdc,xout_dmdc] = lsim(sys_dmdc_d,input_voltages,time,dT);


%stack voltages in reverse

% input_voltages_backwards  = input_voltages(:,2:end);
% omega_backwards = [x2;input_voltages_backwards];
% 
% [utilde_backwards,stilde_backwards,vtilde_backwards] = svd(omega_backwards, 'econ');
% 
% 
% utilde_backwards_trunc = utilde_backwards(:,1:trunc);
% 
% stilde_backwards_trunc = stilde_backwards(1:trunc,1:trunc);
% 
% vtilde_backwards_trunc = vtilde_backwards(:, 1:trunc);
% 
% utilde_backwards_trunc_1 = utilde_backwards_trunc(1:12,:);
% utilde_backwards_trunc_2 = utilde_backwards_trunc(13:16,:);
% 
% 
% Atilde_backwards = x*vtilde_backwards_trunc*pinv(stilde_backwards_trunc)*utilde_backwards_trunc_1';
% 
% Btilde_backwards = x*vtilde_backwards_trunc*pinv(stilde_backwards_trunc)*utilde_backwards_trunc_2';
% 
% A_avg = sqrt(Atilde_backwards*pinv(Atilde));
% 
% B_avg = sqrt(Btilde_backwards*pinv(Btilde));
% 
% sys_dmdc_fb = ss(A_avg,B_avg, sysc.C,sysc.D);

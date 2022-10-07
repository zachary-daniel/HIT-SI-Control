clear; close all; clc;
mdsconnect('landau.hit')
mdsopen('hitsiu', 220926010) %Shot with a noice sphereomak

time = mdsvalue('dim_of(\v_div_1_fspa)');


time = time(4095:length(time),:);
% L2 currents
L2_Current_Flux_1 = mdsvalue('\i_fcoil_1');
L2_Current_Flux_2 = mdsvalue('\i_fcoil_2');
L2_Current_Flux_3 = mdsvalue('\i_fcoil_3');
L2_Current_Flux_4 = mdsvalue('\i_fcoil_4');

%L1 Currents
L1_Current_Flux_1 = mdsvalue('\i_spa_f1');
L1_Current_Flux_2 = mdsvalue('\i_spa_f2');
L1_Current_Flux_3 = mdsvalue('\i_spa_f3');
L1_Current_Flux_4 = mdsvalue('\i_spa_f4');

%Reduce number of samples of L1 Currents
L1_Current_Flux_1 = boxCarAvg(L1_Current_Flux_1, 2);
L1_Current_Flux_2 = boxCarAvg(L1_Current_Flux_2, 2);
L1_Current_Flux_3 = boxCarAvg(L1_Current_Flux_3, 2);
L1_Current_Flux_4 = boxCarAvg(L1_Current_Flux_4, 2);


%Reduce number of Samples of L2 currents

L2_Current_Flux_1 = boxCarAvg(L2_Current_Flux_1, 2);
L2_Current_Flux_2 = boxCarAvg(L2_Current_Flux_2, 2);
L2_Current_Flux_3 = boxCarAvg(L2_Current_Flux_3, 8);
L2_Current_Flux_4 = boxCarAvg(L2_Current_Flux_4, 2);


%Pseudo C Voltage. L2 is in parallel with C so the voltage should be the
%same
C_Voltage_Flux_1 = mdsvalue('\v_div_1_fc');
C_Voltage_Flux_2 = mdsvalue('\v_div_2_fc');
C_Voltage_Flux_3 = mdsvalue('\v_div_3_fc');
C_Voltage_Flux_4 = mdsvalue('\v_div_4_fc');

% SPA waveforms
SPA_Voltage_1 = mdsvalue('\v_div_1_fspa');
SPA_Voltage_2 = mdsvalue('\v_div_2_fspa');
SPA_Voltage_3 = mdsvalue('\v_div_3_fspa');
SPA_Voltage_4 = mdsvalue('\v_div_4_fspa');

%Toroidal currents. This becomes non-zero when breakdown occurs
i_tor_1 = mdsvalue('\I_inj_1');
i_tor_2 = mdsvalue('\I_inj_2');
i_tor_3 = mdsvalue('\I_inj_3');
i_tor_4 = mdsvalue('\I_inj_4');

%Get rid of extra i_tor samples
i_tor_1 = i_tor_1(1:65536, 1);
i_tor_2 = i_tor_2(1:65536, 1);
i_tor_3 = i_tor_3(1:65536, 1);
i_tor_4 = i_tor_4(1:65536, 1);

%Reduce number of i_tor samples
i_tor_1 = boxCarAvg(i_tor_1,2);
i_tor_2 = boxCarAvg(i_tor_2,2);
i_tor_3 = boxCarAvg(i_tor_3,2);
i_tor_4 = boxCarAvg(i_tor_4,2);

%Poloidal angle i_tor measurments
i_tor_000 = mdsvalue('\I_inj_000');
i_tor_090 = mdsvalue('\I_inj_090');
i_tor_180 = mdsvalue('\I_inj_180');
i_tor_270 = mdsvalue('\I_inj_270');

%Reduce number of i_tor_angle samples
i_tor_000 = boxCarAvg(i_tor_000,4);
i_tor_090 = boxCarAvg(i_tor_090,4);
i_tor_180 = boxCarAvg(i_tor_180,4);
i_tor_270 = boxCarAvg(i_tor_270,4);

start = 1;
stop = length(L2_Current_Flux_1);

%%
%Form snapshot matrices
% This is from 1:end-1
X = [L1_Current_Flux_1(start:stop-1,1)'; C_Voltage_Flux_1(start:stop-1,1)'; L2_Current_Flux_1(start:stop-1,1)'; L1_Current_Flux_2(start:stop-1,1)';
    C_Voltage_Flux_2(start:stop-1,1)'; L2_Current_Flux_2(start:stop-1,1)'; L1_Current_Flux_3(start:stop-1,1)'; C_Voltage_Flux_3(start:stop-1,1)';
   L2_Current_Flux_3(start:stop-1,1)'; L1_Current_Flux_4(start:stop-1,1)'; C_Voltage_Flux_4(start:stop-1,1)'; L2_Current_Flux_4(start:stop-1,1)'; 
   i_tor_1(start:stop-1)'; i_tor_2(start:stop-1)'; i_tor_3(start:stop-1)'; i_tor_4(start:stop-1)'; i_tor_000(start:stop-1)'; i_tor_090(start:stop-1)';
   i_tor_180(start:stop-1)'; i_tor_270(start:stop-1)';];

%from 2:end
X2 = [L1_Current_Flux_1(start+1:stop,1)'; C_Voltage_Flux_1(start+1:stop,1)'; L2_Current_Flux_1(start+1:stop,1)'; L1_Current_Flux_2(start+1:stop,1)'; 
    C_Voltage_Flux_2(start+1:stop,1)'; L2_Current_Flux_2(start+1:stop,1)'; L1_Current_Flux_3(start+1:stop,1)'; C_Voltage_Flux_3(start+1:stop,1)';
   L2_Current_Flux_3(start+1:stop,1)'; L1_Current_Flux_4(start+1:stop,1)'; C_Voltage_Flux_4(start+1:stop,1)'; L2_Current_Flux_4(start+1:stop,1)';
    i_tor_1(start+1:stop)'; i_tor_2(start+1:stop)'; i_tor_3(start+1:stop)'; i_tor_4(start+1:stop)'; i_tor_000(start+1:stop)'; i_tor_090(start+1:stop)';
   i_tor_180(start+1:stop)'; i_tor_270(start+1:stop)';];
%%
% Matrix of inputs from 1:end-1
Upsilon = [SPA_Voltage_1(start:stop-1)'; SPA_Voltage_1(start:stop-1)'; SPA_Voltage_1(start:stop-1)'; SPA_Voltage_1(start:stop-1)';];

Omega = [X; Upsilon];

[Utilde, Stilde, Vtilde] = svd(Omega, 'econ');

[Uhat, Shat, Vhat] = svd(X2, 'econ');

U1tilde = Utilde(1:20,:);
U2tilde = Utilde(21:24,:);

A_tilde = X2*Vtilde*pinv(Stilde)*U1tilde';
B_tilde = X2*Vtilde*pinv(Stilde)*U2tilde';

C = [0 ,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];

D = zeros( size(C,1), 4 );

% get rid of spa vals before time = 0
SPA_Voltage_1 = SPA_Voltage_1(4095:length(SPA_Voltage_1));
SPA_Voltage_2 = SPA_Voltage_2(4095:length(SPA_Voltage_2));
SPA_Voltage_3 = SPA_Voltage_3(4095:length(SPA_Voltage_3));
SPA_Voltage_4 = SPA_Voltage_4(4095:length(SPA_Voltage_4));

dT = time(2) - time(1);
sys = ss(A_tilde, B_tilde,C,D, dT);

clear; close all; clc;
mdsconnect('landau.hit')
mdsopen('hitsiu', 220816009)

%I think these are the spa waveforms and not the l1 voltages despite what
%it says in mds
v_spa_f1 = mdsvalue('\v_div_1_fspa');
v_spa_f1 = double(v_spa_f1);
v_spa_f2 = mdsvalue('\v_div_2_fspa');
v_spa_f2 = double(v_spa_f2);
v_spa_f3 = mdsvalue('\v_div_3_fspa');
v_spa_f3 = double(v_spa_f3);
v_spa_f4 = mdsvalue('\v_div_4_fspa');
v_spa_f4 = double(v_spa_f4);
% Current through L2
i_l2_f1 = double(mdsvalue('\i_fcoil_1'));
i_l2_f2 = double(mdsvalue('\i_fcoil_2')); % This current has more channels
i_l2_f3 = double(mdsvalue('\i_fcoil_3'));
i_l2_f4 = double(mdsvalue('\i_fcoil_4'));
%%
% Voltage across L2
v_l2_f1 = double(mdsvalue('\v_div_1_fc'));
v_l2_f2 = double(mdsvalue('\v_div_2_fc'));
v_l2_f3 = double(mdsvalue('\v_div_3_fc'));
v_l2_f4 = double(mdsvalue('\v_div_4_fc'));

%Time variable for whole data set
time = double(mdsvalue('dim_of(\v_div_1_fc)'));

% I think the L1 current isn't right but im gonna check anyways?
i_l1_f1 = double(mdsvalue('\i_spa_f1'));
i_l1_f2 = double(mdsvalue('\i_spa_f2'));
i_l1_f3 = double(mdsvalue('\i_spa_f3'));
i_l1_f4 = double(mdsvalue('\i_spa_f4'));


% boxcar average the l1 currents so they have the same number of samples
% as the voltage dividers
i_l1_f1 = boxCarAvg(i_l1_f1,2);
i_l1_f2 = boxCarAvg(i_l1_f2,2);
i_l1_f3 = boxCarAvg(i_l1_f3,2);
i_l1_f4 = boxCarAvg(i_l1_f4,2);

%boxcar average l2 currents as above
i_l2_f1 = boxCarAvg(i_l2_f1,2);
i_l2_f2 = boxCarAvg(i_l2_f2,2);
i_l2_f3 = boxCarAvg(i_l2_f3,8); % this thing is kinda fucked up so maybe leave it out?
i_l2_f4 = boxCarAvg(i_l2_f4,2);

%Reconstruct cap voltage. Oh wait! You can't because they don't measure the
%L1 voltage correctly :))))))))

%%
% Here goes DMDc

%This ends up giving me an a and a B matrix, but they're prediction for thev_spa_f4
%shot isn't very good. Next step is to try to do shorter intervals and
%mulitple recalculations of teh A and B matrices
%stack data in rows from index 1 to index end-1
start = 10000;
stop_point = start+500;
X = [i_l1_f1(start:stop_point-1,1)'; i_l1_f2(start:stop_point-1,1)'; i_l1_f3(start:stop_point-1,1)'; i_l1_f4(start:stop_point-1,1)'; i_l2_f1(start:stop_point-1,1)'; i_l2_f2(start:stop_point-1,1)'; i_l2_f3(start:stop_point-1,1)'; i_l2_f4(start:stop_point-1,1)'; v_l2_f1(start:stop_point-1,1)'; v_l2_f2(start:stop_point-1,1)'; v_l2_f3(start:stop_point-1,1)'; v_l2_f4(start:stop_point-1,1)';];
x0 = [i_l1_f1(start,1)'; i_l1_f2(start,1)'; i_l1_f3(start,1)'; i_l1_f4(start,1)'; i_l2_f1(start,1)'; i_l2_f2(start,1)'; i_l2_f3(start,1)'; i_l2_f4(start,1)'; v_l2_f1(start,1)'; v_l2_f2(start,1)'; v_l2_f3(start-1,1)'; v_l2_f4(start,1)';];
%stack data from index 2 to end
X2 = [i_l1_f1(start+1:stop_point,1)'; i_l1_f2(start+1:stop_point,1)'; i_l1_f3(start+1:stop_point,1)'; i_l1_f4(start+1:stop_point,1)'; i_l2_f1(start+1:stop_point,1)'; i_l2_f2(start+1:stop_point,1)'; i_l2_f3(start+1:stop_point,1)'; i_l2_f4(start+1:stop_point,1)'; v_l2_f1(start+1:stop_point,1)'; v_l2_f2(start+1:stop_point,1)'; v_l2_f3(start+1:stop_point,1)'; v_l2_f4(start+1:stop_point,1)';];

%stack inputs from 1:end-1
Upsilon = [v_spa_f1(start:stop_point-1)'; v_spa_f2(start:stop_point-1)'; v_spa_f3(start:stop_point-1)'; v_spa_f4(start:stop_point-1)';];

Omega = [X; Upsilon]; %omega matrix with the data and the inputs stacked on each other

%perform SVD on Omega

[Utilde, Stilde, Vtilde] = svd(Omega, 'econ'); %SVD omega matrix

[Uprime, Sprime, Vprime] = svd(X2, 'econ'); %Svd snapshots from the second half of the data

U1tilde = Utilde(1:12,:); % data has 12 rows so take first 12 rows of svd from omega matrix
U2tilde = Utilde(13:16,:); %inputs is four rows so the remaining rows are assigned here
Atilde = Uprime'*X2*Vtilde*pinv(Stilde)*U1tilde'*Uprime; %construct Atilde matrix which should be an approximation of the dynamics
                                 % Eigs of A will show temporal behavior


Btilde = Uprime'*X2*Vtilde*Stilde'*U2tilde';
                         
C = [0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

D = zeros( size(C,1), size(Btilde,2) );

dT = time(2)-time(1);
sys_dmd = ss(Atilde,Btilde,C,D);                                      

sys_dmd_d = c2d(sys_dmd, dT, 'zoh');

%functions below this line

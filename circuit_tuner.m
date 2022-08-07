clear; close all; clc;
mdsconnect('landau.hit');
shot = 220802016;
mdsopen('hitsiu', shot)
Amplitude = 600;
Frequency = 18500;
RunTime = .004;
SampleTime = 1e-7;
Lp = 1.85e-6; 
L1 = 1.4e-6; %H
L1_V = 4.832e-7; % H
L2 = 1.5e-6; %H
%M = L2/5;
%Mw = M/5;
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
% This is for a shot in which the second flux circuit is open and the first
% flux circuit is run @ 600 V.
index = 9946;
index2 = 17943;
t = mdsvalue('dim_of(\v_div_1_fc)');
v_fc_1 = mdsvalue('\v_div_1_fc');
v_fc_2 = mdsvalue('\v_div_2_fc');
v_l1 = mdsvalue('\v_div_1_fspa');
i_fc_1 = mdsvalue('\i_fcoil_1');
i_fc_2 = mdsvalue('\i_fcoil_2');
i_l1 = mdsvalue('\i_spa_f1');
t_i = mdsvalue('dim_of(\i_spa_f1)');
Cap_and_R2_Voltage = v_fc_1(index)+R3*(i_fc_1(index2));
% [nada_v_fc_1, loc_nada_v_fc_1] = findpeaks(-abs(v_fc_1));
di1 = (i_fc_2(index2+1) - i_fc_2(index2))/( t_i(index2+1) - t_i(index2) );
M = (-v_fc_2(index) - i_fc_2(index2)*(-R2-R3))/di1; % This is the mutual inductance of this dumb fucking circuit

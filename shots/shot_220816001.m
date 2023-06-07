clear; close all; clc;

shot = "220816009.hdf5";

time = h5read(shot,"/time");

i_tor = h5read(shot,"/itor");

v_spa_1 = h5read(shot,"/inj/circuit/v_1_fspa");
v_spa_2 = h5read(shot,"/inj/circuit/v_1_fspa");
v_spa_3 = h5read(shot,"/inj/circuit/v_1_fspa");
v_spa_4 = h5read(shot,"/inj/circuit/v_1_fspa");



i_fcoil_1 = h5read(shot,"/inj/circuit/i_fcoil_1");
i_fcoil_2 = h5read(shot,"/inj/circuit/i_fcoil_2");
i_fcoil_3 = h5read(shot,"/inj/circuit/i_fcoil_3");
i_fcoil_4 = h5read(shot,"/inj/circuit/i_fcoil_4");

v_fcoil_1 = h5read(shot,"/inj/circuit/v_fcoil_1");
v_fcoil_2 = h5read(shot,"/inj/circuit/v_fcoil_2");
v_fcoil_3 = h5read(shot,"/inj/circuit/v_fcoil_3");
v_fcoil_4 = h5read(shot,"/inj/circuit/v_fcoil_4");

%Capacitor is in parallel with the flux coil so the voltages are equal
v_cap_1 = v_fcoil_1;
v_cap_2 = v_fcoil_2;
v_cap_3 = v_fcoil_3;
v_cap_4 = v_fcoil_4;

%L1 (series inductor) voltages
v_L1_1 = v_spa_1 - v_cap_1;
v_L1_2 = v_spa_2 - v_cap_2;
v_L1_3 = v_spa_3 - v_cap_3;
v_L1_4 = v_spa_4 - v_cap_4;

% All NaN's in this array?
plot(time,v_fcoil_1)



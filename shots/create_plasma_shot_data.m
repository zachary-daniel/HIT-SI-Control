function [data,shot_num] = create_plasma_shot_data(shot) %pass in shot as form "shot#.hdf5", 
% saves a .mat file where the name of the file is the shot number
% For analyzing plasma shots, we're going to need the voltage circuit data,
% and not the flux circuit data.
str_length = length(shot);
shot_num = extractBefore(shot,str_length-4);

time = h5read(shot,"/time");

i_tor = h5read(shot,"/itor");

v_spa_1 = h5read(shot,"/inj/circuit/v_1_vspa");
v_spa_2 = h5read(shot,"/inj/circuit/v_1_vspa");
v_spa_3 = h5read(shot,"/inj/circuit/v_1_vspa");
v_spa_4 = h5read(shot,"/inj/circuit/v_1_vspa");



i_vcoil_1 = h5read(shot,"/inj/circuit/i_vcoil_1");
i_vcoil_2 = h5read(shot,"/inj/circuit/i_vcoil_2");
i_vcoil_3 = h5read(shot,"/inj/circuit/i_vcoil_3");
i_vcoil_4 = h5read(shot,"/inj/circuit/i_vcoil_4");

v_vcoil_1 = h5read(shot,"/inj/circuit/v_vcoil_1");
v_vcoil_2 = h5read(shot,"/inj/circuit/v_vcoil_2");
v_vcoil_3 = h5read(shot,"/inj/circuit/v_vcoil_3");
v_vcoil_4 = h5read(shot,"/inj/circuit/v_vcoil_4");

%Capacitor is in parallel with the flux coil so the voltages are equal
v_cap_1 = v_vcoil_1;
v_cap_2 = v_vcoil_2;
v_cap_3 = v_vcoil_3;
v_cap_4 = v_vcoil_4;

%L1 currents
i_L1_1 = h5read(shot,"/inj/circuit/i_spa_v1");
i_L1_2 = h5read(shot,"/inj/circuit/i_spa_v2");
i_L1_3 = h5read(shot,"/inj/circuit/i_spa_v3");
i_L1_4 = h5read(shot,"/inj/circuit/i_spa_v4");


%capacitor current
i_cap_1 = i_L1_1 - i_vcoil_1;
i_cap_2 = i_L1_2 - i_vcoil_2;
i_cap_3 = i_L1_3 - i_vcoil_3;
i_cap_4 = i_L1_4 - i_vcoil_4;

%L1 (series inductor) voltages
v_L1_1 = v_spa_1 - v_cap_1;
v_L1_2 = v_spa_2 - v_cap_2;
v_L1_3 = v_spa_3 - v_cap_3;
v_L1_4 = v_spa_4 - v_cap_4;

data = [time,v_spa_1,v_spa_2,v_spa_3,v_spa_4,...
    i_L1_1,v_cap_1,i_vcoil_1,v_L1_1,i_cap_1,v_vcoil_1,i_L1_2,v_cap_2,i_vcoil_2,v_L1_2,i_cap_2,v_vcoil_2,...
    i_L1_3,v_cap_3,i_vcoil_3,v_L1_3,i_cap_3,v_vcoil_3,i_L1_4,v_cap_4,i_vcoil_4,v_L1_4,i_cap_4,v_vcoil_4,i_tor];
save(strcat(shot_num,'_voltage'))
end

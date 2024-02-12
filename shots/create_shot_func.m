function [data,shot_num] = create_shot_func(shot, plasma, flux,varargin) %pass in shot as form "shot#.hdf5", 
% saves a .mat file where the name of the file is the shot number
% shot: string of shot number followed by .hdf5. Eg. 123456789.hdf5
% plasma: boolean indicating whether or not the shot was a vacuum or plasma
% shot. If plasma, will include i_tor measurment
% flux: boolean. If true, will include flux circuit measurments, if false,
% will return voltage circuit measurments
%varargin: can include a file path as a string to pull shots from a
%different folder or directory


str_length = length(shot);
shot_num = extractBefore(shot,str_length-4); %extract the shot number before the .hdf5
plasma_or_vac = ''; % declaration. This will be tacked onto the end of the file name so the
                    % saved file will be denoted as vacuum or plasma
if nargin > 3 %concatenate string onto file name if provided

    shot = strcat(varargin,shot);
    shot = shot{1};
end
time = h5read(shot,"/time");
if plasma == true
    i_tor = h5read(shot,"/itor"); %grab toroidal plasma current
    plasma_or_vac = '_plasma'; %assign for plasma
    data = zeros(30,length(time)); %pre-allocate space for data array
    data(30,:) = i_tor'; %opt-dmd likes feature, then data as rows and cols
else 
    plasma_or_vac = '_vacuum';
    data = zeros(29,length(time)); %one smaller for no i_tor
end
file_name = strcat(shot_num,plasma_or_vac);

if flux == true
    
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
    
    %L1 currents
    i_L1_1 = h5read(shot,"/inj/circuit/i_spa_f1");
    i_L1_2 = h5read(shot,"/inj/circuit/i_spa_f2");
    i_L1_3 = h5read(shot,"/inj/circuit/i_spa_f3");
    i_L1_4 = h5read(shot,"/inj/circuit/i_spa_f4");
    
    
    %capacitor current KCL to get this relation
    i_cap_1 = i_L1_1 - i_fcoil_1;
    i_cap_2 = i_L1_2 - i_fcoil_2;
    i_cap_3 = i_L1_3 - i_fcoil_3;
    i_cap_4 = i_L1_4 - i_fcoil_4;
    
    %L1 (series inductor) voltages. KCL to get these relation
    v_L1_1 = v_spa_1 - v_cap_1;
    v_L1_2 = v_spa_2 - v_cap_2;
    v_L1_3 = v_spa_3 - v_cap_3;
    v_L1_4 = v_spa_4 - v_cap_4;
    
    %assign data
    data(1:29,:) = [time, v_spa_1, v_spa_2, v_spa_3, v_spa_4, ...
    i_L1_1, v_cap_1, i_fcoil_1, v_L1_1, i_cap_1, v_fcoil_1, i_L1_2, v_cap_2, i_fcoil_2, v_L1_2, i_cap_2, v_fcoil_2,...
    i_L1_3, v_cap_3, i_fcoil_3, v_L1_3, i_cap_3, v_fcoil_3, i_L1_4, v_cap_4, i_fcoil_4, v_L1_4, i_cap_4, v_fcoil_4]';
    save(strcat(file_name,'_flux'))
else 
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
    
    
    %capacitor current. KCL to get this relation
    i_cap_1 = i_L1_1 - i_vcoil_1;
    i_cap_2 = i_L1_2 - i_vcoil_2;
    i_cap_3 = i_L1_3 - i_vcoil_3;
    i_cap_4 = i_L1_4 - i_vcoil_4;
    
    %L1 (series inductor) voltages. KCL to get this relation
    v_L1_1 = v_spa_1 - v_cap_1;
    v_L1_2 = v_spa_2 - v_cap_2;
    v_L1_3 = v_spa_3 - v_cap_3;
    v_L1_4 = v_spa_4 - v_cap_4;
    %assign data 
    data(1:29,:) = [time,v_spa_1,v_spa_2,v_spa_3,v_spa_4,...
        i_L1_1,v_cap_1,i_vcoil_1,v_L1_1,i_cap_1,v_vcoil_1,i_L1_2,v_cap_2,i_vcoil_2,v_L1_2,i_cap_2,v_vcoil_2,...
        i_L1_3,v_cap_3,i_vcoil_3,v_L1_3,i_cap_3,v_vcoil_3,i_L1_4,v_cap_4,i_vcoil_4,v_L1_4,i_cap_4,v_vcoil_4]';
    save(strcat(file_name,'_voltage'))
end

end

clear;close all;clc;
files = dir("vacuum_shots_220810");
for i = size(files,1)
    file_name = files(i).name;
    disp(file_name)
    [data,shot] = create_shot_func(file_name,false,true,'vacuum_shots_220810/');
end

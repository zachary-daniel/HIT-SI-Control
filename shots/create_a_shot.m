clear;close all;clc;
files = dir("high_itor_plasma_shots");
for i = 3:size(files,1)
    file_name = files(i).name;
    disp(file_name)
    [data,shot] = create_shot_func(file_name,true,true,'high_itor_plasma_shots/');
end

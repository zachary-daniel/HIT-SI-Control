clear;close all;clc;
files = dir("*.*");
for i = 3:size(files,1)
    str = files(i).name;
    [data,shot] = create_shot_data(str);
end

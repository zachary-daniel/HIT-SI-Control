clear;close all;clc;
[data,shot] = create_shot_func('221201010.hdf5',true,true,'high_itor_plasma_shots/');
load('high_itor_plasma_shots\matlab data\221129011_plasma_flux.mat');
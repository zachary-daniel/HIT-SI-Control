clear;close all;clc;
load("vacuum_shots_220810/220810008_vacuum_flux.mat")
i_L2_1 = data(8,:);
v_L2_1 = data(11,:);
z_1 = hilbert(v_L2_1);
z_2 = hilbert(i_L2_1);
zrat_1 = (imag(z_1) ./ imag(z_2));

plot((zrat_1))
%Just wanted to see the step response of the simulink model to a step
%function. I guess this has reassured me? Whatever
clear; close all; clc; Amplitude1 = 600;
load('All_injectors_vaccum_full_workspace.mat')
sim("Circuit_with_step_inputs.slx");

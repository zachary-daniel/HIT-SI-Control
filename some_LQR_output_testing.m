clear; close all; clc;
SampleTime = 1e-7;
Amplitude = 300000;
RunTime = .004;
time = 0:SampleTime:RunTime;


load("All_injectors_with_LQG_code.mat", "LQR_Outputs")

vals = toSquare(LQR_Outputs(:,1), Amplitude, SampleTime, time);

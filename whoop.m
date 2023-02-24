clear;close all;clc;
load('voltages.mat');
SampleTime = 1e-7;
RunTime = .004;
time = 0:SampleTime:RunTime;
newVoltage = toSquare(voltage, 600, SampleTime, time);
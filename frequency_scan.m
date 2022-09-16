clear; close all; clc;
%load("All_injectors_vaccum_full_workspace.mat")
Amplitude = 600;
Amplitude1 = 600;
%Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7; 
L1 = (8.0141e-7); %H
L2 = 2.0462e-6; %H
M = 3.2951e-7;
Mw = 2.7548e-7;
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
NoisePower = 0;
PhaseAngle1 = 0;
PhaseAngle2 = 0;
PhaseAngle3 = 0;


[MaxCapVoltages, MaxLCurrent, Frequencies] = FrequencyScan(900, 600, 1900,2);

% A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
%      1/Cap, 0, -1/Cap;
%      (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2)];
% B = [1/L1;
%     0;
%     0;];
% C = [0,0,1];



%functions below this line
function [MaxCapVoltageArray, MaxLCurrentArray, frequencies] = FrequencyScan(initFreq, Amplitude, endFreq, steps)
    stepSize = round((endFreq-initFreq)/steps);
    Amp = Amplitude;
    Frequency = initFreq;
    assignin("base", "Frequency", Frequency);
    assignin("base", "Amp", Amp);
    MaxCapVoltageArray = zeros(steps,1);
    MaxLCurrentArray = zeros(steps,1);
    frequencies = zeros(steps,1);
    index = 1;
    for i = initFreq:stepSize:endFreq
        frequencies(index,1) = Frequency;
        Frequency = i;
        assignin("base", "Frequency", Frequency);
        assignin("base", "Amp", Amp);
        open("Use_for_freq_scan.m")
        MaxLCurrentArray(index,1) = ans.L2CurrentMax.signals.values(end);
        index = index+1;
    end
end
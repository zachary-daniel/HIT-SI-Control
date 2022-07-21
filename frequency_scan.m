clear; close all; clc;
Amplitude = 600;
Frequency = 26762;
RunTime = .004;
SampleTime = 1e-7;
Lp = 1.85e-6; 
L1 = 1.4e-6; %H
L1_V = 4.832e-7; % H
L2 = 1.5e-6; %H
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
NoisePower = 0;
M = L2/2;
[MaxCapVoltages, MaxLCurrent, Frequencies] = FrequencyScan(19000, 600, 19010,2);

% A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
%      1/Cap, 0, -1/Cap;
%      (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2)];
% B = [1/L1;
%     0;
%     0;];
% C = [0,0,1];

A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0;
    1/Cap, 0, -1/Cap, 0, 0, 0;
    (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2), 0, 0, -M
    0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
    0, 0, 0, 1/Cap, 0, -1/Cap;
    0, 0, -M, (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2);];
B = [1/L1;
    0;
    0;
    0;
    0;
    0;];
C = [0, 0, 1, 0, 0, 0];

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
        run("All_injectors_vaccum.m")
        MaxLCurrentArray(index,1) = ans.L2CurrentMax.signals.values(end);
        index = index+1;
    end
end
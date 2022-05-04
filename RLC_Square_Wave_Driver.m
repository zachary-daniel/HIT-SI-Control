%%
clear;close all;clc;

RunTime = .004;
% [width, height] = squareInput(19100,600);
open("RLC_model_1.slx")
% Set Variables
t = 0:1e-6:RunTime;
% start = 18000;
% stop = 22000;
i = 21724;
period = 1/i;
steps = 30;
amplitude = 600;
duty = 65;
Model = "RLC_model_1.slx";
input = 600*square(2*pi*i*t,duty);
simin.time = t';
simin.signals.values = input';
sim("RLC_model_1.slx", "StopTime", "RunTime")



%%
clear; close all; clc

%Sim init
RunTime = .004;
% [width, height] = squareInput(19100,600);
open("RLC_model_1.slx")
% Set Variables
t = 0:1e-6:RunTime;
start = 18000;
stop = 22000;
steps = 30;
amplitude = 600;
duty = 65;
Model = "RLC_model_1.slx";

stepSize = round((stop-start)/steps);
    open(Model)
    max_current = zeros(steps,1);
    frequencies = zeros(steps, 1);
    index = 1;
    for i = start:stepSize:stop
        Input = amplitude*square(2*pi*i*t,duty);
        period = 1/i;
        simin.time = t';
        simin.signals.values = Input';
        sim(Model, "StopTime", "RunTime");
        max_current(index, 1) = ans.L2CurrentMax.signals.values(end);
        frequencies(index,1) = i;
        index = index+1;
    end

plot(frequencies, max_current);

%  Input = 600*square(2*pi*18950*t,40);
%  plot(t,Input);
% % 
%  simin.time = t';
%  simin.signals.values = Input';
% 
%  sim("RLC_model_1.slx", "StopTime", "RunTime")

% functions below this line

function [width, height] = squareInput(frequency, amplitude)
    syms x;
    Period = 1/frequency;
    Area = int(amplitude*sin((2*pi*frequency)*x),0,Period/2);
    height = Area/(Period/2);
    width = Period/2;
end

function [max_current, frequencies] = frequency_scan(start, stop, steps, amplitude, duty, time, Modelstr)
    stepSize = round((stop-start)/steps);
    assignin("base", "time", time)
    open(Modelstr)
    max_current = zeros(steps,1);
    frequencies = zeros(steps, 1);
    index = 1;
    for i = start:stepSize:stop
        Input = amplitude*square(2*pi*i*time,duty);
        assignin("base", "Input", Input);
        simin.time = time';
        simin.signals.values = Input';
        sim(Modelstr, "StopTime", "RunTime");
        max_current(index, 1) = ans.L2CurrentMax.signals.values(end);
        frequencies(index,1) = i;
        index = index+1;
    end
end
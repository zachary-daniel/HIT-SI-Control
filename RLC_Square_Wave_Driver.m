clear;close all;clc;

[width, height] = squareInput(19100,600);
open("RLC_model_1.slx")
 % t = ans.tout;
t = 0:1e-6:.004;

Input = 600*square(2*pi*19100*t,50);
plot(t,Input);

simin.time = t';
simin.signals.values = Input';

sim("RLC_model_1.slx")

% functions below this line

function [width, height] = squareInput(frequency, amplitude)
    syms x;
    Period = 1/frequency;
    Area = int(amplitude*sin((2*pi*frequency)*x),0,Period/2);
    height = Area/(Period/2);
    width = Period/2;
end
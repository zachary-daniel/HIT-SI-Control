clear; close all; clc;
load("All_injectors_vaccum.mat")
dT = 1e-7;
RunTime = .004;
time = 0:dT:RunTime';
plot(time/10, voltage)
xlim([0, .0001])
[nada, nada_locs] = findpeaks(-abs(voltage));
[shifted_voltage] = phaseShift(voltage, 180, nada_locs);
plot(time, voltage)
hold on
plot(time, shifted_voltage)
xlim([.003 .004])
%functions below this line


function [newVals] = phaseShift(vals, shift, loc_nada) % takes shift in degrees
    shift = mod(shift, 360);
    newVals = zeros(size(vals));
    period = loc_nada(3) - loc_nada(1);
    cut_off = round(period*(shift/360));
    newVals(1:end-cut_off) = vals(cut_off+1:end);
    newVals(end-cut_off:end) = -vals(end:-1:end-cut_off);
    
end
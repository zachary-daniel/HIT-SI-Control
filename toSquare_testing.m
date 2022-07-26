clear; close all; clc;
load("All_injectors_vaccum.mat")
SampleTime = 1e-7;
RunTime = .004;
time = 0:SampleTime:RunTime;
Frequency = 19000;
Amplitude = 600;
voltage = 600*sin(2*pi*Frequency*time);
plot(time, voltage)
[peaks, loc_peaks] = findpeaks(voltage); % Peaks
[troughs, loc_troughs] = findpeaks(-voltage); %Troughs
[peak_times] = locsToTimes(loc_peaks, time); %Peak times
[trough_times] = locsToTimes(loc_troughs, time); %Trough times
[nada, loc_nada] = findpeaks(-(abs(voltage))); % nada
[nada_times] = locsToTimes(loc_nada, time); % nada times
troughs = -troughs;

[newVoltages] = toSquare(voltage, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude, SampleTime)

% Functions below this line


function [ts] = locsToTimes(locs, time)
    for i = 1:length(locs)
        t = locs(i);
        locs(i) = time(t);
    end
    ts = locs;
end


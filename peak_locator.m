close all; clc;
time = out.tout;
current = out.L2Current.signals.values;
[peaks, loc_peaks] = findpeaks(current);
[troughs, loc_troughs] = findpeaks(-current);
[peak_times] = locsToTimes(loc_peaks, time);
[trough_times] = locsToTimes(loc_troughs, time);
[zeros, loc_zeros] = findpeaks(-(abs(current)));
[zero_times] = locsToTimes(loc_zeros, time);
figure(1)
plot(time, current, "Color", "g")
hold on
plot(trough_times, -troughs, "ok")
hold on
plot(peak_times, peaks, "or")
hold on
plot(zero_times, zeros, "ob")
figure(2)
plot(time, -abs(current))

%Functions below this line

function [ts] = locsToTimes(locs, time)
    for i = 1:length(locs)
        t = locs(i);
        locs(i) = time(t);
    end
    ts = locs;
end
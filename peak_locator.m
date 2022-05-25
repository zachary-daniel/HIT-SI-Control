close all; clc;
%% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
Amplitude = 600;
Frequency = 19000;
RunTime = .004;
sim("RLC_Sin_To_Square_w_H_Bridge.slx", "StopTime", "RunTime")
time = ans.tout;
voltage = ans.simout.signals.values;

%Locate peaks, troughs, and nada, and change the location data into time
%data

[peaks, loc_peaks] = findpeaks(voltage); % Peaks
[troughs, loc_troughs] = findpeaks(-voltage); %Troughs
[peak_times] = locsToTimes(loc_peaks, time); %Peak times
[trough_times] = locsToTimes(loc_troughs, time); %Trough times
[nada, loc_nada] = findpeaks(-(abs(voltage))); % nada
[nada_times] = locsToTimes(loc_nada, time); % nada times
troughs = -troughs;
%Plot the wave with peaks, nada, troughs
% figure(1)
% title("Peaks, Troughs, nada of Sine Wave")
% xlabel("Time")
% ylabel("Amplitude")
% plot(time, voltage, "Color", "g")
% hold on
% plot(trough_times, troughs, "ok")
% hold on
% plot(peak_times, peaks, "or")
% hold on
% plot(nada_times, nada, "ob")
% legend("voltage", "Troughs", "Peaks", "nada")
% figure(2)
% plot(time, -abs(voltage))

% Equivalent Sine wave area in dis
[newVoltages] = toSquare(voltage, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude);
open("RLC_with_HBrdige_forward_circuit.slx")
simin.time = time;
simin.signals.values = newVoltages;
sim("RLC_with_HBrdige_forward_circuit.slx", "StopTime", "RunTime");
%  plot(time, newVoltages)
%  hold on
%  plot(time, voltage, "r")






%Functions below this line

function [ts] = locsToTimes(locs, time)
    for i = 1:length(locs)
        t = locs(i);
        locs(i) = time(t);
    end
    ts = locs;
end

function [newValues] = toSquare(values, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude)
    Magic_number = 4/pi;
    num_waves = length(nada);
    Areas = zeros(num_waves,1);
    first_maxima = 0;
    second_maxima = 0;
    first_maxima_times = 0;
    second_maxima_times = 0;
    if trough_times(1) < peak_times(1)
        first_maxima = troughs;
        first_maxima_times = trough_times;
        second_maxima = peaks;
        second_maxima_times = peak_times;
    else
        first_maxima = peaks;
        first_maxima_times = peak_times;
        second_maxima = troughs;
        second_maxima_times = trough_times;
    end
    
    extrema = zeros(length(troughs)+length(peaks), 1);
    extrema_times = zeros(length(trough_times)+length(peak_times),1);
    index = 1;
    for i = 1:2:length(extrema)-1
        extrema(i,1) = first_maxima(index);
        extrema(i+1,1) = second_maxima(index);
        extrema_times(i,1) = first_maxima_times(index);
        extrema_times(i+1,1) = second_maxima_times(index);
        index=index+1;
    end
    extrema(end) = first_maxima(end);
    extrema_times(end) = first_maxima_times(end);
    for i = 1:length(extrema)-1
        base = nada_times(i+1) - nada_times(i);
        Areas(i,1) = .5*base*extrema(i)*Magic_number;
    end
    Areas(end) = .5*(nada_times(end) - nada_times(length(nada_times)-1))*extrema(end)*Magic_number;
    %Make Square Waves here
    % Make new voltage 
    newValues = zeros(length(values),1);
    Amplitude = 600;
    for i = 1:length(Areas)
        width = abs(Areas(i)/Amplitude);
        b1 = round((extrema_times(i)-width/2)/1e-7);
        b2 = round((extrema_times(i) + width/2)/1e-7);
        if Areas(i) < 0
            newValues(b1:b2,1) = -Amplitude;
        else 
            newValues(b1:b2,1) = Amplitude;
        end
    end
end
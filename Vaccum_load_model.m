clear;close all; clc;
%% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
Amplitude = 600;
Frequency = 19000;
RunTime = .004;
SampleTime = 1e-7;
Lp = 1.85e-6; 
L1 = 1.4e-6; %H
L2 = 1.5e-6; %H
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
NoisePower = 0;
A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
     1/Cap, 0, -1/Cap;
     (1/L2)*R2, 1/L2, (-1/L2)*(R3+R2)];
B = [1/L1;
    0;
    0;];
C = [0,0,1];

D = zeros( size(C,1), size(B,2) );
states = {'x1', 'x2', 'x3'};
inputs = {'v(t)'}; 
outputs = {'x3'};
%continuous time system
sysc = ss(A,B,C,D, "statename", states,"inputname", inputs, "outputname", outputs);
A = sysc.A;
B = sysc.B;
C = sysc.C;
D = sysc.D;

%Discrete Time System
sys_d = c2d(sysc, dT, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

G = eye(3);
H = zeros(3,3);

Q = diag(.001*ones(1,3)); % disturbance covariance
R = 1; % Noise covariance
time = (0:SampleTime:RunTime);
backwards_vals = (Amplitude*sin(time*Frequency*2*pi));

scaling_amplitude = 1; %(sin(time*2*pi*Frequency/20)); This is a test that the method works for sin waves that scale up and down in amplitude

%%
backwards_circuit_simin.time = (time)';
time = 0:SampleTime:RunTime;
backwards_circuit_simin.signals.values = (scaling_amplitude.*backwards_vals)';
sim("RLC_Sin_To_Square_Backwards.slx", "StopTime", "RunTime")
% sim("RLC_Sin_To_Square_Backwards.slx")
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
open("Vaccum_circuit_w_load_forward_7_11.slx")
simin.time = time;
simin.signals.values = newVoltages;
sim("Vaccum_circuit_w_load_forward_7_11.slx", "StopTime", "RunTime");
figure()
plot(time, newVoltages)
hold on
plot(time, voltage, "r")
xlabel("Time")
ylabel("Input Voltage")
title("Conversion of Sin wave to Square Wave for H-Bridge")
legend("Output of Reverse Circuit", "Input to H-Bridge")

% Grab Klaman Filter outputs and get corresponding values
L2_Current = ans.L2Current.signals.values;
C_Voltage = ans.CVoltage.signals.values;
L1_Current = ans.L1Current.signals.values;
L2_Current_Approx = ans.KalmanFilterOutputs.signals.values(:,3);
C_Voltage_Approx = ans.KalmanFilterOutputs.signals.values(:,2);
L1_Current_Approx = ans.KalmanFilterOutputs.signals.values(:,1);

% Plot L2 vs. L2 Approx
figure()
plot(time, L2_Current, "LineWidth",2)
hold on
plot(time, L2_Current_Approx, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L2 Current")
legend("True L2 Current with Noise", "Denoised L2 Current from KF", "Location", "northwest")




% Plot C vs. C approx
figure()
plot(time, C_Voltage, "LineWidth", 5)
hold on
plot(time, C_Voltage_Approx, "Linewidth",.25)
xlabel("Time")
ylabel("Voltage")
title("Noisey Output vs. Kalman Filter Output for C Voltage")
legend("True C Voltage with Noise", "Denoised C Voltage from KF", "Location", "northwest")



% Plot L1 Current vs. L1 Approx
figure()
plot(time, L1_Current, "LineWidth", 2)
hold on
plot(time, L1_Current_Approx, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L1 Current")
legend("True L1 Current with Noise", "Denoised L1 Current from KF", "Location", "northwest")

% %%
% close all;
% figure()
% % plot(time, L2_Current)
% % hold on
% % plot(time, L1_Current)
% figure()
% sysFullOutput = ss(A,B,eye(3), D);
% [y,t] = lsim(sysFullOutput, newVoltages, time);
% plot(time, y(:,1));
% legend("calculated", "measured")
% hold on
% plot(time, L1_Current);
% title("L1")
% figure()
% plot(time, y(:,2))
% hold on
% plot(time, C_Voltage)
% title("C")
% legend("calculated", "measured")
% figure()
% plot(time, y(:,3))
% hold on
% plot(time, L2_Current)
% title("L2")
% legend("calculated", "measured")
% % legend("Lsim L1", "Lsim C", "Lsim L2", "Sim L1", "Sim C", "Sim L2");
% %Functions below this line

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
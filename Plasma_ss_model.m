clear; close all; clc;
% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
run('init_vals_for_plasma_model_code.m')
% mdsopen('hitsiu', 220926010);
Amplitude = 600;
Amplitude1 = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .002;
SampleTime = 1e-7; 
L1 = (8.0141e-7); %H
L2  =2.0462e-6; %H
Lp = L2/2; % H Placeholder value until the impedance of the plasma can be calculated
M = .161*L2; % Coupling coefficient
Mw = .1346*L2;% Coupling coefficient
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
Mp = M; % H This is a placeholder value until an actual measurment of the mutual inductance of the plasma can be measured
NoisePower = .0;
PhaseAngle1 = 0;
PhaseAngle2 = 0;
PhaseAngle3 = 0;

%Get inductance of the plasma
i_plasma = mdsvalue('\i_tor_spaavg');
%v_plasma = mdsvalue('\');



%%%State space model that includes the Plasma. A will be 13x13.

%Entries for x3
% scaler_x3 = 
% x3a = 


%entries for x6
scaler_x6 = -1/((L2 - Mw) *(L2 - 2* M + Mw) *(L2 *Lp + 2* Lp* M - 4* (Mp.^2) + Lp *Mw) );
x6a = L2 * Lp * M * R2 - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2;
x6b = L2 * Lp * M - L2 * (Mp.^2) - Lp * M * Mw + (Mp.^2) * Mw;
x6c = -L2 * Lp * M * R2 + L2 * Mp.^2 * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2 - L2 * Lp * M * R3 + L2 * (Mp.^2) * R3 + Lp * M * Mw * R3 - (Mp.^2) * Mw * R3;
x6d = -(L2.^2) * Lp * R2 + 2 * Lp * (M.^2) * R2 + 3 * L2 * (Mp.^2) * R2 - 4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + (Mp.^2) * Mw * R2;
x6e = -(L2.^2) * Lp + 2 * Lp * (M.^2) + 3 * L2 * (Mp.^2) - 4 * M * (Mp.^2) - L2 * Lp * Mw + (Mp.^2) * Mw;
x6f = (L2.^2) * Lp * R2 - 2 * Lp * (M.^2) * R2 - 3 * L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - (Mp.^2) * Mw * R2 + (L2.^2) * Lp * R3 - 2 * Lp * (M.^2) * R3 - 3 * L2 * (Mp.^2) * R3 + 4 * M * (Mp.^2) * R3 + L2 * Lp * Mw * R3 - (Mp.^2) * Mw * R3;
x6g = -2 * Lp * (M.^2) * R2 - L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - 3 * (Mp.^2) * Mw * R2 + Lp * (Mw.^2) * R2;
x6h = -2 * Lp * (M.^2) - L2 * (Mp.^2) + 4 * M * (Mp.^2) + L2 * Lp * Mw - 3 * (Mp.^2) * Mw + Lp * (Mw.^2);
x6i = 2 * Lp * (M.^2) * R2 + L2 * (Mp.^2) * R2 - 4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + 3 * (Mp.^2) * Mw * R2 - Lp * (Mw.^2) * R2 + 2 * Lp * M.^2 * R3 + L2 * (Mp.^2) * R3 - 4 * M * (Mp.^2) * R3 - L2 * Lp * Mw * R3 + 3 * (Mp.^2) * Mw * R3 - Lp * (Mw.^2) * R3;
x6j = L2 * Lp * M * R2 - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2;
x6k = L2 * Lp * M  - L2 * (Mp.^2) - Lp * M * Mw + (Mp.^2) * Mw;
x6l = -L2 * Lp * M * R2 + L2 * (Mp.^2) * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2 - L2 * Lp * M * R3 + L2 * (Mp.^2) * R3  + Lp * M * Mw * R3 - (Mp.^2) * Mw * R3;
x6m = 0;
%entries for x9
scaler_x9 = -1/ ( (L2 - Mw) * (L2 - 2 * M + Mw) * (L2 * Lp + 2 * Lp * M - 4 * (Mp.^2) + Lp * Mw) );
x9a = L2 * Lp * M * R2 - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2;
x9b = L2 * Lp * M - L2 * (Mp.^2) - Lp * M * Mw + (Mp.^2) * Mw;
x9c = -L2 * Lp * M * R2 + L2 * (Mp.^2) * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2 - L2 * Lp * M * R3 + L2 * (Mp.^2) * R3 + Lp * M * Mw * R3 - (Mp.^2) * Mw * R3;
x9d = -2 * Lp * (M.^2) * R2 - L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - 3 * (Mp.^2) * Mw * R2 + Lp * (Mw.^2) * R2;
x9e = -2 * Lp * (M.^2) - L2 * (Mp.^2) + 4 * M * (Mp.^2) + L2 * Lp * Mw - 3 * (Mp.^2) * Mw + Lp * (Mw.^2);
x9f = 2 * Lp * (M.^2) * R2 + L2 * (Mp.^2) * R2 - 4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + 3 * (Mp.^2) * Mw * R2 - Lp * (Mw.^2) * R2 + 2 * Lp * (M.^2) * R3 + L2 * (Mp.^2) * R3 - 4 * M * (Mp.^2) * R3 - L2 * Lp * Mw * R3 + 3 * (Mp.^2) * Mw * R3 - Lp * (Mw.^2) * R3;
x9g = -(L2.^2) * Lp * R2 + 2 * Lp * (M.^2) * R2 + 3 * L2 * (Mp.^2) * R2 - 4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + (Mp.^2) * Mw * R2;
x9h = -(L2.^2) * Lp + 2 * Lp * (M.^2) + 3 * L2 * (Mp.^2) - 4 * M * (Mp.^2) - L2 * Lp * Mw + (Mp.^2) * Mw;
x9i = (L2.^2) * Lp * R2 - 2 * Lp * (M.^2) * R2 - 3 * L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - (Mp.^2) * Mw * R2 + (L2.^2) * Lp * R3 - 2 * Lp * (M.^2) * R3 - 3 * L2 * (Mp.^2) * R3 + 4 * M * (Mp.^2) * R3 + L2 * Lp * Mw * R3 - (Mp.^2) * Mw * R3;
x9j = L2 * Lp * M * R2  - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2;
x9k = L2 * Lp * M - L2 * (Mp.^2) - Lp * M * Mw + (Mp.^2) * Mw;
x9l = -L2 * Lp * M * R2 + L2 * (Mp.^2) * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2 - L2 * Lp * M * R3 + L2 * (Mp.^2) * R3 + Lp * M * Mw * R3 - (Mp.^2) * Mw * R3;
x9m = 0;

%x12 entries
scaler_x12 = 1/( (L2 - Mw) * (L2 - 2 * M + Mw) * (L2 * Lp + 2 * Lp * M - 4 * (Mp.^2) + Lp * Mw) );

x12a = 2 * Lp * (M.^2) * R2 + L2 * (Mp.^2) * R2 - 4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + 3 * (Mp.^2) * Mw * R2 - Lp * (Mw.^2) * R2;
x12b = 2 * Lp * (M.^2) + L2 * (Mp.^2) - 4 * M * (Mp.^2) - L2 * Lp * Mw + 3 * (Mp.^2) * Mw - Lp * (Mw.^2);
x12c = -2 * Lp * (M.^2) * R2 - L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - 3 * (Mp.^2) * Mw * R2 + Lp * (Mw.^2) * R2 - 2 * Lp * (M.^2) * R3 - L2 * (Mp.^2) * R3 + 4 * M * (Mp.^2) * R3 + L2 * Lp * Mw * R3 - 3 * (Mp.^2) * Mw * R3 + Lp * (Mw.^2) * R3;
x12d = -L2 * Lp * M * R2 + L2 * (Mp.^2) * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2;
x12e = -L2 * Lp * M + L2 * (Mp.^2) + Lp * M * Mw - (Mp.^2) * Mw;
x12f = L2 * Lp * M * R2 - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2 +L2 * Lp * M * R3- L2 * (Mp.^2) * R3- Lp * M * Mw * R3 + (Mp.^2) * Mw * R3;
x12g = -L2 * Lp * M * R2 + L2 * (Mp.^2) * R2 + Lp * M * Mw * R2 - (Mp.^2) * Mw * R2;
x12h = -L2 * Lp * M + L2 * (Mp.^2) + Lp * M * Mw - (Mp.^2) * Mw;
x12i = L2 * Lp * M * R2 - L2 * (Mp.^2) * R2 - Lp * M * Mw * R2 + (Mp.^2) * Mw * R2 + L2 * Lp * M * R3 - L2 * (Mp.^2) * R3 - Lp * M * Mw * R3 + (Mp.^2) * Mw * R3;
x12j = (L2.^2) * Lp * R2 - 2 * Lp * (M.^2) * R2 - 3 * L2 * (Mp.^2) * R2 + 4 * M * (Mp.^2) * R2 + L2 * Lp * Mw * R2 - (Mp.^2) * Mw * R2;
x12k = (L2.^2) * Lp - 2 * Lp * (M.^2) - 3 * L2 * (Mp.^2) + 4 * M * (Mp.^2) + L2 * Lp * Mw - (Mp.^2) * Mw;
x12l =  -(L2.^2) * Lp * R2 + 2 * Lp * (M.^2) * R2 + 3 * L2 * (Mp.^2) * R2 -4 * M * (Mp.^2) * R2 - L2 * Lp * Mw * R2 + (Mp.^2) * Mw * R2 -(L2.^2) * Lp * R3 + 2 * Lp * (M.^2) * R3 + 3 * L2 * (Mp.^2) * R3 - 4 * M * (Mp.^2) * R3 - L2 * Lp * Mw * R3 + (Mp.^2) * Mw * R3;
x12m = 0;

%%
%x13 entries


%build A matrix
% A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%    1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%     scaler_x3*x3a, scaler_x3*x3b,scaler_x3*x3c, scaler_x3*x3d, scaler_x3*x3e, scaler_x3*x3f, scaler_x3*x3g, scaler_x3*x3h, scaler_x3*x3i, scaler_x3*x3j, scaler_x3*x3k, scaler_x3*x3l, scaler_x3*x3m;
%     0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0;
%     0, 0, 0,  1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0;
%     scaler_x6*x6a, scaler_x6*x6b, scaler_x6*x6c, scaler_x6*x6d, scaler_x6*x6e, scaler_x6*x6f, scaler_x6*x6g, scaler_x6*x6h, scaler_x6*x6i, scaler_x6*x6j, scaler_x6*x6k, scaler_x6*x6l, scaler_x6*x6m;
%     0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0;
%     0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0;
%    scaler_x9*x9a, scaler_x9*x9b, scaler_x9*x9c, scaler_x9*x9d, scaler_x9*x9e, scaler_x9*x9f, scaler_x9*x9g, scaler_x9*x9h, scaler_x9*x9i, scaler_x9*x9j, scaler_x9*x9k, scaler_x9*x9l, scaler_x9*x9m;
%     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
%     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap;
%    scaler_x12*x12a, scaler_x12*x12b, scaler_x12*x12c, scaler_x12*x12d,scaler_x12*x12e, scaler_x12*x12f, scaler_x12*x12g, scaler_x12*x12h, scaler_x12*x12i, scaler_x12*x12j, scaler_x12*x12k, scaler_x12*x12l, scaler_x12*x12m;
%     scaler_x13*x13a, scaler_x13*x13b, scaler_x13*x13c, scaler_x13*x13d,scaler_x13*x13e, scaler_x13*x13f, scaler_x13*x13g, scaler_x13*x13h, scaler_x13*x13i, scaler_x13*x13j, scaler_x13*x13k, scaler_x13*x13l, scaler_x13*x13m];
% 
% %Build B matrix
% B = [1/L1, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 1/L1, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 1/L1, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 1/L1;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0];
% %build C matrix where we observe all L2 currents and the plasma current
% C = [0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
%     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1;
%     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1;
%     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1];
% 
% %D matrix is all zeros
% D = zeros( size(C,1), size(B,2) );
% 
% %build continuous time model
% sys_plasma = ss(A,B,C,D);
% %discrete time model with dT 
% sys_plasma_d = c2d(sys_plasma, dT, 'zoh');
% %Create discrete time matrices
% Ad = sys_plasma_d.A;
% Bd = sys_plasma_d.B;
% Cd = sys_plasma_d.C;
% Dd = sys_plasma_d.D;
% 
% %check ctrb and obsv
% control = rank(ctrb(Ad, Bd)) == size(A,1); %gucci
% observe = rank(obsv(Ad, Cd)) == size(A,1); %gucci
% 
% Q = .001; %diag(.001*ones(1,size(A, 1))); % disturbance covariance
% R = diag(10*ones(1,size(B,2))); % Noise covariance
 time = (0:SampleTime:RunTime);
 backwards_vals = (Amplitude*sin(time*Frequency*2*pi));
% 
 scaling_amplitude = 1; %(sin(time*2*pi*Frequency/20)); This is a test that the method works for sin waves that scale up and down in amplitude
% 


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

[newVoltages] = toSquare(voltage, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude, SampleTime);
[newVoltageShift1] = phaseShift(newVoltages, PhaseAngle1, loc_nada);
[newVoltageShift2] = phaseShift(newVoltages, PhaseAngle2, loc_nada);
[newVoltageShift3] = phaseShift(newVoltages, PhaseAngle3, loc_nada);
injector1 = zeros(size(newVoltages));
%Kalman filter bitch
% [kalmf, L, P] = kalman(sys_plasma_d, Q, R, 0);
% %%
% syskf = ss(Ad-L*Cd, [Bd L],eye(size(Ad,1)), 0*[Bd L], dT);
% [y,t] = lsim(sys_plasma_d, [voltage voltage voltage voltage], time);
% [yk, tk] = lsim(syskf, [voltage voltage voltage voltage, y], time);


%%

open("All_injectors_with_plasma.slx")
shiftedSignal1.time = time;
shiftedSignal1.signals.values = newVoltageShift1;
shiftedSignal2.time = time;
shiftedSignal2.signals.values = newVoltageShift2;
shiftedSignal3.time = time;
shiftedSignal3.signals.values = newVoltageShift3;
simin.time = time;
simin.signals.values = newVoltages;
sim("All_injectors_with_plasma.slx", "StopTime", "RunTime");
figure()
plot(time, newVoltages)
hold on
plot(time, voltage, "r")
xlabel("Time")
ylabel("Input Voltage")
title("Conversion of Sin wave to Square Wave for H-Bridge")
legend("Input to H-Bridge", "Output of Reverse Circuit")

% % Grab Klaman Filter outputs and get corresponding values
% L2_Current_Flux_1 = ans.L2CurrentFlux1.signals.values;
% C_Voltage_Flux_1 = ans.CVoltageFlux1.signals.values;
% L1_Current_Flux_1 = ans.L1CurrentFlux1.signals.values;
% L2_Current_Approx_Flux_1= ans.KalmanFilterOutputsFlux1.signals.values(:,3);
% C_Voltage_Approx_Flux_1 = ans.KalmanFilterOutputsFlux1.signals.values(:,2);
% L1_Current_Approx_Flux_1 = ans.KalmanFilterOutputsFlux1.signals.values(:,1);
% 
% L2_Current_Flux_2 = ans.L2CurrentFlux2.signals.values;
% C_Voltage_Flux_2 = ans.CVoltageFlux2.signals.values;
% L1_Current_Flux_2 = ans.L1CurrentFlux2.signals.values;
% 
% L2_Current_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,6);
% C_Voltage_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,5);
% L1_Current_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,4);
% 
% L1_Current_Approx_Flux_3 = ans.KalmanFilterOutputsFlux1.signals.values(:,7);
% L1_Current_Flux_3 = ans.L1CurrentFlux3.signals.values;
% C_Voltage_Flux_3 = ans.CVoltageFlux3.signals.values;
% L2_Current_Flux_3 = ans.L2CurrentFlux3.signals.values;
% 
% L1_Current_Flux_4 = ans.L1CurrentFlux4.signals.values;
% C_Voltage_Flux_4 = ans.CVoltageFlux4.signals.values;
% L2_Current_Flux_4 = ans.L2CurrentFlux4.signals.values;

% Plot L2 vs. L2 Approx
% figure()
% plot(time, L2_Current_Flux_1, "LineWidth",2)
% hold on
% plot(time, L2_Current_Approx_Flux_1, "Linewidth", .25)
% xlabel("Time")
% ylabel("Current")
% title("Noisey Output vs. Kalman Filter Output for L2 Current (Flux 1)")
% legend("True L2 Current with Noise", "Denoised L2 Current from KF", "Location", "northwest")
% 
% 
% 
% %%  
% 
% % Plot C vs. C approx
% figure()
% plot(time, C_Voltage_Flux_1,  "LineWidth", 1)
% hold on
% plot(time, C_Voltage_Approx_Flux_1,'r--', "Linewidth",2)
% 
% xlabel("Time (s)", 'fontsize', 20)
% ylabel("Voltage (Volts)", 'fontsize', 20)
% title("Noisey Output vs. Kalman Filter Output for C Voltage (Flux 1)", 'fontsize', 30)
% legend("Cap Voltage with Noise", "Denoised Cap Voltage from KF", "Location", "northeast", 'fontsize', 15)
% 
% 
% %%
% % Plot L1 Current vs. L1 Approx
% figure()
% plot(time, L1_Current_Flux_1, "LineWidth", 2)
% hold on
% plot(time, L1_Current_Approx_Flux_1, "Linewidth", .25)
% xlabel("Time")
% ylabel("Current")
% title("Noisey Output vs. Kalman Filter Output for L1 Current (Flux 1)")
% legend("True L1 Current with Noise", "Denoised L1 Current from KF", "Location", "northwest")
% 
% % L2 Current Flux 2
% figure()
% plot(time, L2_Current_Flux_2, "LineWidth",2)
% hold on
% plot(time, L2_Current_Approx_Flux_2, "Linewidth", .25)
% xlabel("Time")
% ylabel("Current")
% title("Noisey Output vs. Kalman Filter Output for L2 Current (Flux 2)")
% legend("True L2 Current with Noise", "Denoised L2 Current from KF", "Location", "northwest")
% 
% %C voltage Flux 2
% figure()
% plot(time, C_Voltage_Flux_2, "LineWidth", 5)
% hold on
% plot(time, C_Voltage_Approx_Flux_2, "Linewidth",.25)
% xlabel("Time")
% ylabel("Voltage")
% title("Noisey Output vs. Kalman Filter Output for C Voltage (Flux 2)")
% legend("True C Voltage with Noise", "Denoised C Voltage from KF", "Location", "northwest")
% 
% % L1 Current Flux 2
% figure()
% plot(time, L1_Current_Flux_2, "LineWidth", 2)
% hold on
% plot(time, L1_Current_Approx_Flux_2, "Linewidth", .25)
% xlabel("Time")
% ylabel("Current")
% title("Noisey Output vs. Kalman Filter Output for L1 Current (Flux 2)")
% legend("True L1 Current with Noise", "Denoised L1 Current from KF", "Location", "northwest")
% 
% 
% 

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
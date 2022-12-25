clear; close all; clc;
% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
mdsopen('hitsiu', 220802016);
Amplitude = 600;
Amplitude1 = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7; 
L1 = (8.0141e-7); %H
L2 = 2.0462e-6; %H
M = .161*L2; % Coupling coefficient
Mw = .1346*L2;% Coupling coefficient
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
NoisePower = .1;
PhaseAngle1 = 90;
PhaseAngle2 = 180;
PhaseAngle3 = 270;

scalar1 = 1/((L2-Mw)*( (L2.^2) - (4*M.^2) + 2*L2*Mw+ (Mw.^2) )); %Scale factor in front of the entries to the A matrix that are affected by mutual inductance

x3a =  (-L2.^2)*R2+(2*M.^2)*R2-L2*Mw*R2;
x3b = (-L2.^2)+(2*M.^2)-L2*Mw;
x3c = (L2.^2)*R2-(2*M.^2)*R2+L2*Mw*R2+(L2.^2)*R3-(2*M.^2)*R3+L2*Mw*R3;
x3d = L2*M*R2 - M*Mw*R2;
x3e = L2*M-M*Mw;
x3f = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3g = L2*M*R2-M*Mw*R2;
x3h = L2*M-M*Mw;
x3i = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3j = -2*(M.^2)*R2+L2*Mw*R2+(Mw.^2)*R2;
x3k = -2*(M.^2)+L2*Mw+Mw.^2;
x3l = 2*(M.^2)*R2-L2*Mw*R2-(Mw.^2)*R2+2*R3*(M.^2)-L2*Mw*R3-R3*Mw.^2;

%Entries for x6 in A matrix
x6a = -L2*M*R2+M*Mw*R2;
x6b = -L2*M+M*Mw;
x6c = L2*M*R2-M*Mw*R2+L2*M*R3-M*Mw*R3;
x6d = R2*(L2.^2)-2*R2*(M.^2)+L2*Mw*R2;
x6e = (L2.^2)-2*(M.^2)+L2*Mw;
x6f = -R2*(L2.^2)+2*R2*(M.^2)-L2*Mw*R2-R3*(L2.^2)+2*R3*(M.^2)-L2*Mw*R3;
x6g = 2*R2*(M.^2)-L2*Mw*R2-R2*(Mw.^2); 
x6h = 2*(M.^2)-L2*Mw-(Mw.^2);
x6i = -2*R2*(M.^2)+L2*Mw*R2+R2*(Mw.^2)-2*R3*(M.^2)+L2*Mw*R3+R3*(Mw.^2);
x6j = -L2*M*R2+M*Mw*R2;
x6k = -L2*M+M*Mw;
x6l = L2*M*R2-M*Mw*R2+L2*M*R3-M*Mw*R3;

%Entries for x9 in A matrix
x9a = -L2*M*R2 + M*Mw*R2;
x9b = -L2* M + M* Mw;
x9c = L2* M* R2 - M* Mw* R2 * + L2* M *R3  - M* Mw* R3 ;
x9d = 2* (M.^2) *R2 - L2*Mw* R2 - (Mw.^2) *R2;
x9e = 2 *(M.^2) - L2* Mw - (Mw.^2);
x9f = -2* (M.^2)* R2+ L2 *Mw *R2 + (Mw.^2) *R2 - 2* (M.^2) *R3 + L2 *Mw* R3 + (Mw.^2) *R3;
x9g =(L2.^2) *R2 - 2* (M.^2) *R2 + L2 *Mw *R2;
x9h = (L2.^2) - 2* (M.^2) + L2 *Mw;
x9i = -(L2.^2) *R2 + 2* (M.^2)* R2 - L2 *Mw *R2 - (L2.^2) *R3 + 2* (M.^2)* R3- L2 *Mw *R3;
x9j = -L2 *M *R2 + M*Mw*R2;
x9k = -L2 *M + M *Mw ;
x9l = L2*M*R2 - M* Mw *R2 + L2 *M *R3 - M* Mw *R3;

%Entries for x12 in A matrix
x12a = -2* (M.^2) *R2+ L2* Mw* R2 +(Mw.^2)* R2;
x12b = -2 *(M.^2) + L2 *Mw + (Mw.^2);
x12c = 2 *(M.^2) *R2 - L2 *Mw *R2 - (Mw.^2) *R2 + 2 *(M.^2) *R3 - L2 *Mw *R3 - (Mw.^2) *R3;
x12d = L2 *M *R2 - M *Mw *R2;
x12e = L2 *M - M *Mw;
x12f = -L2 *M *R2 + M *Mw *R2 - L2 *M *R3 + M *Mw* R3;
x12g = L2 *M *R2 - M *Mw *R2;
x12h = L2 *M - M *Mw;
x12i = -L2 *M *R2 + M *Mw *R2 - L2 *M *R3 + M *Mw *R3;
x12j = (-L2.^2) *R2 + 2 *(M.^2) *R2 - L2 *Mw *R2;
x12k = (-L2.^2) + 2 *(M.^2) - L2* Mw;
x12l = (L2.^2) *R2 - 2 *(M.^2) *R2 + L2 *Mw *R2 + (L2.^2) *R3 - 2*(M.^2)* R3 + L2 *Mw *R3;

A = [((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   -scalar1*x3a, -scalar1*x3b,-scalar1*x3c, -scalar1*x3d, -scalar1*x3e, -scalar1*x3f, -scalar1*x3g, -scalar1*x3h, -scalar1*x3i, -scalar1*x3j, -scalar1*x3k, -scalar1*x3l;
   0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0, 0, 0, 0;
   0, 0, 0,  1/Cap, 0, -1/Cap, 0, 0, 0, 0, 0, 0;
   scalar1*x6a, scalar1*x6b, scalar1*x6c, scalar1*x6d, scalar1*x6e, scalar1*x6f, scalar1*x6g, scalar1*x6h, scalar1*x6i, scalar1*x6j, scalar1*x6k, scalar1*x6l;
   0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap, 0, 0, 0;
   scalar1*x9a, scalar1*x9b, scalar1*x9c, scalar1*x9d, scalar1*x9e, scalar1*x9f, scalar1*x9g, scalar1*x9h, scalar1*x9i, scalar1*x9j, scalar1*x9k, scalar1*x9l;
   0, 0, 0, 0, 0, 0, 0, 0, 0, ((-1/L1)*(R1+R2)), -1/L1, R2*1/L1;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap, 0, -1/Cap;
   -scalar1*x12a, -scalar1*x12b, -scalar1*x12c, -scalar1*x12d,-scalar1*x12e, -scalar1*x12f, -scalar1*x12g, -scalar1*x12h, -scalar1*x12i, -scalar1*x12j, -scalar1*x12k, -scalar1*x12l];







B = [1/L1, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 1/L1, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 1/L1, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 1/L1;
    0, 0, 0, 0;
    0, 0, 0, 0;];

% B = [1/L1, 0;
%     0, 0;
%     0, 0;
%     0, 1/L1;
%     0, 0;
%     0, 0;];

% C = [0 0 1 0 0 0;
%     0 0 0 0 0 1;];

C = [0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

D = zeros( size(C,1), size(B,2) );
states = {'x1', 'x2', 'x3', 'x1b', 'x2b', 'x3b', 'x1c', 'x2c', 'x3c', 'x1d', 'x2d', 'x3d'};
inputs = {'v(t)'}; 
outputs = {'x3'};
%%
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

Q = .001; %diag(.001*ones(1,size(A, 1))); % disturbance covariance
R = diag(10*ones(1,size(B,2))); % Noise covariance
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
[nada, loc_nada] = findpeaks(-(abs(voltage))); % nada

[newVoltages] = toSquare(voltage, Amplitude, SampleTime, time);
[newVoltageShift1] = phaseShift(newVoltages, PhaseAngle1, loc_nada);
[newVoltageShift2] = phaseShift(newVoltages, PhaseAngle2, loc_nada);
[newVoltageShift3] = phaseShift(newVoltages, PhaseAngle3, loc_nada);
injector1 = zeros(size(newVoltages));
%Kalman filter bitch
[kalmf, L, P] = kalman(sys_d, Q, R, 0);

syskf = ss(Ad-L*Cd, [Bd L],eye(12), 0*[Bd L], dT);
[y,t] = lsim(sys_d, [voltage voltage voltage voltage], time);
[yk, tk] = lsim(syskf, [voltage voltage voltage voltage, y], time);


%%

open("Vaccum_circuits_all_injectors.slx")
shiftedSignal1.time = time;
shiftedSignal1.signals.values = newVoltageShift1;
shiftedSignal2.time = time;
shiftedSignal2.signals.values = newVoltageShift2;
shiftedSignal3.time = time;
shiftedSignal3.signals.values = newVoltageShift3;
simin.time = time;
simin.signals.values = newVoltages;
sim("Vaccum_circuits_all_injectors.slx", "StopTime", "RunTime");
figure()
plot(time, newVoltages)
hold on
plot(time, voltage, "r")
xlabel("Time")
ylabel("Input Voltage")
title("Conversion of Sin wave to Square Wave for H-Bridge")
legend("Input to H-Bridge", "Output of Reverse Circuit")

% Grab Klaman Filter outputs and get corresponding values
L2_Current_Flux_1 = ans.L2CurrentFlux1.signals.values;
C_Voltage_Flux_1 = ans.CVoltageFlux1.signals.values;
L1_Current_Flux_1 = ans.L1CurrentFlux1.signals.values;
L2_Current_Approx_Flux_1= ans.KalmanFilterOutputsFlux1.signals.values(:,3);
C_Voltage_Approx_Flux_1 = ans.KalmanFilterOutputsFlux1.signals.values(:,2);
L1_Current_Approx_Flux_1 = ans.KalmanFilterOutputsFlux1.signals.values(:,1);

L2_Current_Flux_2 = ans.L2CurrentFlux2.signals.values;
C_Voltage_Flux_2 = ans.CVoltageFlux2.signals.values;
L1_Current_Flux_2 = ans.L1CurrentFlux2.signals.values;

L2_Current_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,6);
C_Voltage_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,5);
L1_Current_Approx_Flux_2 = ans.KalmanFilterOutputsFlux1.signals.values(:,4);

L1_Current_Approx_Flux_3 = ans.KalmanFilterOutputsFlux1.signals.values(:,7);
L1_Current_Flux_3 = ans.L1CurrentFlux3.signals.values;
C_Voltage_Flux_3 = ans.CVoltageFlux3.signals.values;
L2_Current_Flux_3 = ans.L2CurrentFlux3.signals.values;

L1_Current_Flux_4 = ans.L1CurrentFlux4.signals.values;
C_Voltage_Flux_4 = ans.CVoltageFlux4.signals.values;
L2_Current_Flux_4 = ans.L2CurrentFlux4.signals.values;

% Plot L2 vs. L2 Approx
figure()
plot(time, L2_Current_Flux_1, "LineWidth",2)
hold on
plot(time, L2_Current_Approx_Flux_1, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L2 Current (Flux 1)")
legend("True L2 Current with Noise", "Denoised L2 Current from KF", "Location", "northwest")




% Plot C vs. C approx
figure()
plot(time, C_Voltage_Flux_1,  "LineWidth", .5)
hold on
plot(time, C_Voltage_Approx_Flux_1,'r--', "Linewidth",2)
C_Voltage_No_Noise = ans.TrueCVoltage.signals.values; % Cap voltage with no noise
plot(time, C_Voltage_No_Noise, 'g:', 'Linewidth', 2)
xlabel("Time")
ylabel("Voltage")
title("Noisey Output vs. Kalman Filter Output for C Voltage (Flux 1)")
legend("Cap Voltage with Noise", "Denoised Cap Voltage from KF", 'Cap Voltage no Noise', "Location", "northwest")



% Plot L1 Current vs. L1 Approx
figure()
plot(time, L1_Current_Flux_1, "LineWidth", 2)
hold on
plot(time, L1_Current_Approx_Flux_1, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L1 Current (Flux 1)")
legend("True L1 Current with Noise", "Denoised L1 Current from KF", "Location", "northwest")

% L2 Current Flux 2
figure()
plot(time, L2_Current_Flux_2, "LineWidth",2)
hold on
plot(time, L2_Current_Approx_Flux_2, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L2 Current (Flux 2)")
legend("True L2 Current with Noise", "Denoised L2 Current from KF", "Location", "northwest")

%C voltage Flux 2
figure()
plot(time, C_Voltage_Flux_2, "LineWidth", 5)
hold on
plot(time, C_Voltage_Approx_Flux_2, "Linewidth",.25)
xlabel("Time")
ylabel("Voltage")
title("Noisey Output vs. Kalman Filter Output for C Voltage (Flux 2)")
legend("True C Voltage with Noise", "Denoised C Voltage from KF", "Location", "northwest")

% L1 Current Flux 2
figure()
plot(time, L1_Current_Flux_2, "LineWidth", 2)
hold on
plot(time, L1_Current_Approx_Flux_2, "Linewidth", .25)
xlabel("Time")
ylabel("Current")
title("Noisey Output vs. Kalman Filter Output for L1 Current (Flux 2)")
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






% This script will be used to preinitialize values for the plasma model
clear; close all; clc;
% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
mdsopen('hitsiu', 220802016);
Amplitude = 600;
Amplitude1 = 600;
Frequency = 19000;%double(mdsvalue('\sihi_freq'));
RunTime = .002;
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
PhaseAngle1 = 0;
PhaseAngle2 = 0;
PhaseAngle3 = 0;



time = (0:SampleTime:RunTime);
backwards_vals = (Amplitude*sin(time*Frequency*2*pi));

scaling_amplitude = 1; %(sin(time*2*pi*Frequency/20)); This is a test that the method works for sin waves that scale up and down in amplitude




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


[newVoltages] = toSquare(voltage, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude, SampleTime);
[newVoltageShift1] = phaseShift(newVoltages, PhaseAngle1, loc_nada);
[newVoltageShift2] = phaseShift(newVoltages, PhaseAngle2, loc_nada);
[newVoltageShift3] = phaseShift(newVoltages, PhaseAngle3, loc_nada);
injector1 = zeros(size(newVoltages));





shiftedSignal1.time = time;
shiftedSignal1.signals.values = newVoltageShift1;
shiftedSignal2.time = time;
shiftedSignal2.signals.values = newVoltageShift2;
shiftedSignal3.time = time;
shiftedSignal3.signals.values = newVoltageShift3;
simin.time = time;
simin.signals.values = newVoltages;
sim("init_vals_for_plasma_model.slx", "StopTime", "RunTime");

%get final vals
%flux 1
%%


%Get final value of each of the states

L2_Current_Flux_1_end_vac = ans.L2CurrentFlux1.signals.values(end);
C_Voltage_Flux_1_end_vac = ans.CVoltageFlux1.signals.values(end);
L1_Current_Flux_1_end_vac = ans.L1CurrentFlux1.signals.values(end);

%flux 2
L2_Current_Flux_2_end_vac = ans.L2CurrentFlux2.signals.values(end);
C_Voltage_Flux_2_end_vac = ans.CVoltageFlux2.signals.values(end);
L1_Current_Flux_2_end_vac = ans.L1CurrentFlux2.signals.values(end);

%flux 3
L2_Current_Flux_3_end_vac = ans.L2CurrentFlux3.signals.values(end);
C_Voltage_Flux_3_end_vac = ans.CVoltageFlux3.signals.values(end);
L1_Current_Flux_3_end_vac = ans.L1CurrentFlux3.signals.values(end);

%flux 4
L2_Current_Flux_4_end_vac = ans.L2CurrentFlux4.signals.values(end);
C_Voltage_Flux_4_end_vac = ans.CVoltageFlux4.signals.values(end);
L1_Current_Flux_4_end_vac = ans.L1CurrentFlux4.signals.values(end);




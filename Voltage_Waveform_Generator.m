clear; close all;clc
Amplitude = 400;
SampleTime = 1e-6;
time = (0:SampleTime:.004);
frequency = 19000;
voltage = sin(2*pi*frequency*time);
PhaseAngle1 = 90;
PhaseAngle2 = 180;
PhaseAngle3 = 270;
[newVoltages] = toSquare(voltage, Amplitude, SampleTime, time);
[newVoltageShift1] = phaseShift(newVoltages, PhaseAngle1);
[newVoltageShift2] = phaseShift(newVoltages, PhaseAngle2);
[newVoltageShift3] = phaseShift(newVoltages, PhaseAngle3);

voltage_inputs = [newVoltages,newVoltageShift1,newVoltageShift2,newVoltageShift3];
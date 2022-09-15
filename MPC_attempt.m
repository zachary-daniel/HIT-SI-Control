clear;close all; clc;
% Initialize data and plot to make sure everything is gucci
% declare time and voltage for our data
load("All_injectors_vaccum_full_workspace.mat")
% mdsopen('hitsiu', 220802016);
penalty = 0;
Amplitude = 600;
Amplitude1 = 600;
Frequency = 19100;%double(mdsvalue('\sihi_freq'));
RunTime = .004;
SampleTime = 1e-7; 
L1 = (8.0141e-7); %H
L2 = 2.0462e-6; %H
M = 3.2951e-7;
Mw = 2.7548e-7;
Cap = 96e-6; % F
R1 = .0025; %Ohm
R2 = .005; % Ohm
R3 = .005;% Ohm
dT = 1e-7;
NoisePower = .05;
PhaseAngle1 = 0;
PhaseAngle2 = 0;
PhaseAngle3 = 0;

scalar1 = 1/((L2-Mw)*(L2.^2-4*M.^2+2*L2*Mw+Mw.^2)); %Scale factor in front of the entries to the A matrix that are affected by mutual inductance

x3a =  (-L2.^2)*R2+(2*M.^2)*R2-L2*Mw*R2;
x3b = (-L2.^2)+2*M.^2-L2*Mw;
x3c = (L2.^2)*R2-(2*M.^2)*R2+L2*Mw*R2+(L2.^2)*R3-(2*M.^2)*R3+L2*Mw*R3;
x3d = L2*M*R2 - M*Mw*R2;
x3e = L2*M-M*Mw;
x3f = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3g = L2*M*R2-M*Mw*R2;
x3h = L2*M-M*Mw;
x3i = -L2*M*R2+M*Mw*R2-L2*M*R3+M*Mw*R3;
x3j = -2*(M.^2)*R2+L2*Mw*R2+(Mw.^2)*R2;
x3k = (-2*(M.^2)+L2*Mw+Mw.^2);
x3l = 2*(M.^2)*R2-L2*Mw*R2-(Mw.^2)*R2+2*R3*M.^2-L2*Mw*R3-R3*Mw.^2;

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
R = diag(.1*ones(1,size(B,2))); % Noise covariance
time = (0:SampleTime:RunTime);
backwards_vals = (Amplitude*sin(time*Frequency*2*pi));

scaling_amplitude = 1; %(sin(time*2*pi*Frequency/20)); This is a test that the method works for sin waves that scale up and down in amplitude



%%

% backwards_circuit_simin.time = (time)';
% time = 0:SampleTime:RunTime;
% backwards_circuit_simin.signals.values = (scaling_amplitude.*backwards_vals)';
% sim("RLC_Sin_To_Square_Backwards.slx", "StopTime", "RunTime")
% % sim("RLC_Sin_To_Square_Backwards.slx")
% time = ans.tout;
% voltage = ans.simout.signals.values;

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
[kalmf, L, P] = kalman(sys_d, Q, R, 0);
J = @(ref_signal, output, current_switch, next_switch) ((ref_signal-output).^2) + penalty* abs(current_switch-next_switch);
syskf = ss(Ad-L*Cd, [Bd L],eye(12), 0*[Bd L], dT);
[y,t] = lsim(sys_d, [voltage voltage voltage voltage], time);
[yk, tk] = lsim(syskf, [voltage voltage voltage voltage, y], time);



%% MPC Loop
clc;
ref_signal = voltage; %reference signal to track
horizon = 29*dT; %size of step horizon in terms of time samples. 29 seems to be the magic number? Not really sure why...
step_size = horizon/dT; %size of one horizon in terms of indices 
total_steps = RunTime/dT;  % number of time samples
switch_position = 0;
current_states = diag(zeros(size(A,1)));
allInputs = zeros(size(t));
num_iters = 0;
reference = L2_Current_Flux_1;
%The change I want to make is for the MPC to look forward in the reference
%signal farther than the update time. As in, if the rate at which I perform
%a control input is every 24dt's I want the MPC to look at the reference
%signal 36dt's in the future instead of 24dt's. Turns out that shit don't
%work :)))))
for i = round(step_size:step_size:(total_steps-step_size)) % I is always one horizon ahead of where the simulation is
    current_time = (i-step_size)*dT; %get starting time for upcmoming lsim command 
    [next_input, switch_position, results] = MPC(current_states, reference(i), sys_d, switch_position, dT, horizon,current_time, Amplitude, J, RunTime);

        % the reference point at the next time step, state space model in
    % discrete time, current switch position, sample rate in seconds, the
    % horizon (how far out you are trying to predict), the current time,
    % amplitude of SPA, and the cost function
    current_states = lsim(syskf, [next_input, next_input, next_input, next_input, results],current_time:dT:(i)*dT-dT, current_states); % update value of current states for next iteration
    current_states = current_states(end,:); %update current states
    allInputs(1+(length(next_input)*num_iters):length(next_input)*(num_iters+1),1) = next_input; % collect all inputs
    num_iters = num_iters+1;
end
ouch = lsim(sys_d, [allInputs allInputs allInputs allInputs], time);
plot(ouch(:,1))
hold on
plot(L2_Current_Flux_1)
title("Output from Simulation using MPC vs. Reference Signal");
xlabel('Time (s)');
ylabel('Current (Amps)')
legend('Output from Simulation', 'Reference Signal');
ylim([-1.5e4, 1.5e4]);
function [ts] = locsToTimes(locs, time)
    for i = 1:length(locs)
        t = locs(i);
        locs(i) = time(t);
    end
    ts = locs;
end

function [newValues] = toSquare(values, nada, troughs, peaks, nada_times, trough_times, peak_times, Amplitude, SampleTime)
    Magic_number = 4/pi; % Magic number for getting area under sin curve without integrating
    num_waves = length(nada); % number of individual peaks/troughs
    Areas = zeros(num_waves,1); % pre-allocate area array
    first_maxima = 0;
    second_maxima = 0;
    first_maxima_times = 0;
    second_maxima_times = 0;
    last_extrema = 0;
    last_extrema_time = 0;
    if trough_times(1) < peak_times(1) % decide if first extrema is a peak or a trough
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
    if trough_times(end) < peak_times(end)
        last_extrema = peaks(end);
        last_extrema_time = peak_times(end);
    else
        last_extrema = troughs(end);
        last_extrema_time = trough_times(end);
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
    extrema(end) = last_extrema;
    extrema_times(end) = last_extrema_time;
    for i = 1:length(nada_times)-1
        base = nada_times(i+1) - nada_times(i);
        Areas(i,1) = .5*base*extrema(i)*Magic_number;
    end
    if length(nada_times) < length(extrema)
        Areas(end) = -.5*(nada_times(end) - nada_times(length(nada_times)-1))*extrema(end)*Magic_number;

    else
        Areas(end) = .5*(nada_times(end) - nada_times(length(nada_times)-1))*extrema(end)*Magic_number;
    end
    %Make Square Waves here
    % Make new voltage 
    newValues = zeros(length(values),1);
    Amplitude = 600;
    for i = 1:length(Areas)
        width = abs(Areas(i)/Amplitude);
        b1 = round((extrema_times(i)-width/2)/SampleTime);
        b2 = round((extrema_times(i) + width/2)/SampleTime);
        if Areas(i) < 0
            newValues(b1:b2,1) = -Amplitude;
        else 
            newValues(b1:b2,1) = Amplitude;
        end
    end
end



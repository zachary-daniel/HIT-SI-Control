clear;close all;clc;
% 
 
MODEL = "RLC_model_2_input_finder";
% sim("RLC_model_2_input_finder.slx")

%t = ans.ScopeDataVoltages.time;
%V = ans.ScopeDataVoltages.signals.values;

[MaxCapVoltages, MaxLCurrent, Frequencies] = frequencyScan(19000, 600, 19010, 2, MODEL);
plot(Frequencies, MaxLCurrent)
xlabel("Frequencies")
ylabel("Peak Current Through the Inductor")
% figure(1);
% plot(t,V)
% title("Desired Output")
% xlabel("time")
% ylabel("voltage")




%  
%  CapVoltage = ans.OutputC.signals.values;
%  LVoltage = ans.OutputL.signals.values; 
  %MaxCapVoltage = trapz(ans.MaxCapValues.signals.values);

% figure(2)
% plot(t, CapVoltage)
% title("Voltage Across the Capacitor")
% xlabel("time")
% ylabel("voltage")


% figure(3)
% plot(t, LVoltage)
% title("Voltage Across the Inductor")
% xlabel("time")
% ylabel("voltage")


%Functions below this Line
function [MaxCapVoltageArray, MaxLCurrentArray, frequencies] = frequencyScan(initFreq, Amplitude, endFreq, steps, model)
    open(model)
    stepSize = round((endFreq-initFreq)/steps);
    Amp = Amplitude;
    frequency = initFreq;
    assignin("base", "frequency", frequency);
    assignin("base", "Amp", Amp);
    MaxCapVoltageArray = zeros(steps,1);
    MaxLCurrentArray = zeros(steps,1);
    frequencies = zeros(steps,1);
    index = 1;
    for i = initFreq:stepSize:endFreq
        frequencies(index,1) = frequency;
        frequency = i;
        assignin("base", "frequency", frequency);
        assignin("base", "Amp", Amp);
        set_param("RLC_model_2_input_finder/Desired Output", "Frequency", "frequency");
        set_param("RLC_model_2_input_finder/Desired Output", "Amplitude", "Amp");
        sim("RLC_model_2_input_finder.slx")
        MaxLCurrentArray(index,1) = ans.L2CurrentMax.signals.values(end);
        index = index+1;
    end
end

function [width, height] = squareInput(frequency, amplitude)
    Period = 1/frequency;
    Area = int(amplitude*sin(2*pi*frequency),0,Period/2);
    height = Area/(Period/2);
    width = Period/2;
end

    

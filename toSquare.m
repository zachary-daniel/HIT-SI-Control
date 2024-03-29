function [newValues] = toSquare(values, Amplitude, SampleTime, time) %Inputs are the waveform you want to convert to a square wave, The amplitude of the square wave you want
%Sample time, and a time vector for the length of the signal passed in.                                                                  
                                                                     
Magic_number = 4/pi; % Magic number for getting area under sin curve without integrating

    % Get the zero times, peak times, and trough times, as well as the
    % trough and peak values.
    [peaks, loc_peaks] = findpeaks(values); % Peaks
    [troughs, loc_troughs] = findpeaks(-values); %Troughs
    [peak_times] = locsToTimes(loc_peaks, time); %Peak times
    [trough_times] = locsToTimes(loc_troughs, time); %Trough times
    [nada, loc_nada] = findpeaks(-(abs(values))); % nada
    [nada_times] = locsToTimes(loc_nada, time); % nada times
    troughs = -troughs;

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
    % Make new values 
    newValues = zeros(length(values),1);

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

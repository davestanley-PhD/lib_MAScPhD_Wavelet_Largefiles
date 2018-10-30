


function [datafilt frequencies_removed] = smartfilter (datatimes, data, plot_filter_comparison, notch_bandwidth, window_hertz, spike_threshold_multiplier, start_freq, raw_interval)

%   plot_filter_comparison;     Set to 1 to automatically generate a comparison plot
%   start_freq;                 Ignore frequencies below this value.
%   window_hertz;               Window size measured in measured in hertz
%   spike_threshold_multiplier; Consider it a spike if it is this many times larger than the
%                               average surrounding power
%   notch_bandwidth;            After spikes are identified, use this bandwidth of
%                               the notch filter to remove them.
%
%   ** Units **
%   All units of time are expected to be in seconds and frequencies are in hertzs
%
%   --Optional Parameters--
%   raw_interval - Allows you to hard code in some specific frequencies you
%                   would like removed. Expected format is as a column vector of frequencies
%   == Output ==
%   datafilt            Filtered data
%   frequencies_removed List of frequencies that have been removed




% % % % % Some default values to try
%     raw_interval = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
%     notch_bandwidth = 0.1;
%     window_hertz = 60;  %Hz
%     spike_threshold_multiplier = 10;
%     start_freq = 100; %Hz --> Start searching at this frequency


    %Use smart filtering to remove excess mechanical noise
    [f fft_val] = daveFFT(datatimes, data, 1);
    temp = round(length(f)/2); f = f(1:temp); fft_val = fft_val(1:temp);
    
    if ~exist('raw_interval','var')
        raw_interval = [];
    end

    raw_interval = [raw_interval identify_fft_noise(f, fft_val, start_freq, window_hertz, spike_threshold_multiplier)];
    df = f(2) - f(1);

    
    interval = zeros(length(raw_interval), 2);
    for i = 1:length(raw_interval)
%        interval(i,:) =  [max((raw_interval(i) - 2*df), 0) (raw_interval(i) + 2*df)];   %no negative intervals
       interval(i,:) =  [max((raw_interval(i) - notch_bandwidth), 0) (raw_interval(i) + notch_bandwidth)];   %no negative intervals
    end

    datafilt = qif(datatimes, data, interval);
    frequencies_removed = raw_interval;
    
    if (plot_filter_comparison)
        plot_comparison (datatimes, data, datafilt);
    end
    
end



function plot_comparison (datatimes, data, datafiltered)
    [f fft_val] = daveFFT(datatimes, data, 1);
    temp = round(length(f)/2); f = f(1:temp); fft_val = fft_val(1:temp);
    
    [f2 fft_val2] = daveFFT(datatimes, datafiltered, 1);
    temp = round(length(f2)/2); f2 = f2(1:temp); fft_val2 = fft_val2(1:temp);
    
    ufiltered_power = davePower(data);
    filtered_power = davePower (datafiltered);
    percent_of_original = filtered_power / ufiltered_power * 100;
    
    figure; subplot(211)
    plot (datatimes, data - mean(data)); hold on
    plot (datatimes, datafiltered - mean(datafiltered),'r');
%     legend (['Unfiltered Power=' num2str(ufiltered_power, '%e')], ['Filtered Power=' num2str(filtered_power, '%e')]);   
    legend (['Unfiltered Power=' num2str(ufiltered_power, '%e')], ['Filtered Power is ' num2str(percent_of_original) '% of original']);
    xlabel ('Time (s)');
    
    subplot (212)
    plot (f, abs(fft_val).^2, 'r'); hold on;
    plot (f2, abs(fft_val2).^2, 'b');
%     axis ([min(f) max(f) 0 var(fft_val)]);
    
    legend ('Unfiltered FFT', 'Filtered FFT');
    xlabel ('Frequency (Hz)');
    
end
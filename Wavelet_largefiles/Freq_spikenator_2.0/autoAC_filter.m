


function [datafilt frequencies_removed] = autoac (datatimes, data, thr, windw, notch_bandwidth, hardcode_bands, start_freq, plot_filter_comparison)

%   plot_filter_comparison;     Set to 1 to automatically generate a comparison plot
%   start_freq;                 Ignore frequencies below this value.
%   windw;                      Width of the scanning window measured in hertz
%   thr;                        Consider it a spike if it is this many times larger than the
%                               ambient signal in windw
%   notch_bandwidth;            After spikes are identified, use this bandwidth of
%                               the notch filter to remove them.
%   hardcode_bands;             Allows you to hard code in some specific frequencies you
%                               would like removed. Expected format is as a column
%                               vector of frequencies
%
%   ** Units **
%   All units of time are to be in seconds and frequencies are in hertzs
%
%   == Output ==
%   datafilt            Filtered data
%   frequencies_removed List of frequencies that have been removed

%   David Stanley, May 2009
%   *** For academic/research use. Please acknowledge author. ***


if nargin < 2
    fprintf ('Must enter data series and time series together. Assuming integer time series.');
    data = datatimes;
    datatimes = 1:length(data);
    thr = 10;
    windw = 60;  %Hz
    notch_bandwidth = 0.1; %Hz
    hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 3
    thr = 10;
    windw = 60;  %Hz
    notch_bandwidth = 0.1; %Hz
    hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 4
    windw = 60;  %Hz
    notch_bandwidth = 0.1; %Hz
    hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 5
    notch_bandwidth = 0.1; %Hz
    hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 6
    hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 7
    start_freq = 100; %Hz --> Disregard spikes below this frequency
    plot_filter_comparison = 0;
elseif nargin < 8
    plot_filter_comparison = 0;
end



    %Use smart filtering to remove excess mechanical noise
    [f fft_val] = daveFFT(datatimes, data, 1);
    temp = round(length(f)/2); f = f(1:temp); fft_val = fft_val(1:temp);
    
    if ~exist('hardcode_bands','var')
        hardcode_bands = [];
    end

    hardcode_bands = [hardcode_bands identify_fft_noise(f, fft_val, start_freq, windw, thr)];
    df = f(2) - f(1);

    
    interval = zeros(length(hardcode_bands), 2);
    for i = 1:length(hardcode_bands)
%        interval(i,:) =  [max((hardcode_bands(i) - 2*df), 0) (hardcode_bands(i) + 2*df)];   %no negative intervals
       interval(i,:) =  [max((hardcode_bands(i) - notch_bandwidth), 0) (hardcode_bands(i) + notch_bandwidth)];   %no negative intervals
    end

    datafilt = qif(datatimes, data, interval);
    frequencies_removed = hardcode_bands;
    
    if (plot_filter_comparison)
        plot_comparison (datatimes, data, datafilt);
    end
    
end



function plot_comparison (datatimes, data, datafiltered)
    [f fft_val] = daveFFT(datatimes, data, 1);
    temp = round(length(f)/2); f = f(1:temp); fft_val = fft_val(1:temp);
    
    [f2 fft_val2] = daveFFT(datatimes, datafiltered, 1);
    temp = round(length(f2)/2); f2 = f2(1:temp); fft_val2 = fft_val2(1:temp);
    
    ufiltered_power = davePower(data - mean(data));
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



%Quick ideal filter

function dout = qif (datatimes, data, interval)

    ts1 = timeseries(data,datatimes);
    ts1filt = idealfilter (ts1, interval, 'notch');
    dout = ts1filt.data;
    
end



function noise_bands = identify_fft_noise (f, fft, start_freq, window_hertz, spike_threshold_multiplier)

start_freq; % Only search for noise at frequencies above this one
window_hertz; % window size measured in measured in hertz
spike_threshold_multiplier; % Consider it a spike if it is this many times larger than the average surrounding power

%Convert this to the indicies of f
df = f(2) - f(1);
window = round(window_hertz / df);

%Change FFT coefficients to power spectrum
powerSpec = (abs(fft).^2);

noise_bands = [];
start_loc = find (f >= start_freq, 1, 'first'); % Identify location of starting frequency
for i = start_loc:window:(length(f)-window)
    meanpower = mean(abs(powerSpec(i:(i+window-1))));  %Approxixmate power of power spectrum for this particular window
    noise_bands = [noise_bands (find(powerSpec(i:(i+window-1)) > (spike_threshold_multiplier*meanpower)) + (i-1))];   %Find if there are any spikes where the power is many times greater than the background noise
    
end

noise_bands = noise_bands * df;

end


function y = davePower (x)     %Take 2nd (non-central) moment

N = length(x);
y = sum(x.^2) / N;


end


function [f X] = daveFFT (t_input, x_input, scale_freq)

tstep = t_input(2)-t_input(1);
min_t = min(t_input);
max_t = max(t_input);

t = min_t:tstep/scale_freq:max_t;
N = length(t);
dt = tstep/scale_freq;
df = 1/(N*dt);
%t = (0:N-1)*dt;
f = df * (0:N-1);

x = interp1(t_input,x_input,t);
X = fft (x)/N;      %Not entirely sure why we need to divide by N, but the amplitudes seem correct when we run this on a standard exponential function

end




function [f X] = dave_bin_FFT (t_input, x_input, bin_duration)

tstep = t_input(2)-t_input(1);

length_t = length(t_input);
length_bin = round(bin_duration/tstep);

if (length_bin > length_t)
    fprintf ('FFT bin duration is longer than dataset. Decreasing bin size');
	length_bin = length_t;
    
end

nbins = floor(length_t/length_bin)

X = zeros(1,length_bin);
for i=1:nbins
    t_bin = (0:(length_bin-1))*tstep;
    x_bin = x_input( ((i-1)*(length_bin) + 1):(i*length_bin) );
    
    [f_bin X_bin] = daveFFT (t_bin, x_bin, 1);
    X = X + X_bin;
end

X = X / nbins;  % Return the average power spectrum
f = f_bin;

end














load Demitre_testdata.mat

raw_interval = [60 180 300 420 540 660 780 900 1500];    % Some default electrical AC filtering
notch_bandwidth = 0.1; %Hz
window_hertz = 60;  %Hz
spike_threshold_multiplier = 10;
start_freq = 100; %Hz --> Start searching at this frequency
plot_filter_comparison = 1;

 
[datafilt freqremoved]= smartfilter (datatimes, data, plot_filter_comparison, notch_bandwidth, window_hertz, spike_threshold_multiplier, start_freq, raw_interval);


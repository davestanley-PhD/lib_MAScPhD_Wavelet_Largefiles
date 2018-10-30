

load datac.mat

hardcode_bands = [60 180 300 420 540 660 780 900 1500];    % We hard code in some common AC harmonics to remove. Add your own here!
notch_bandwidth = 0.1; %Hz    Size of the notch to use
windw = 60;  %Hz      Width of the scanning window. Generally you want to make the window less than the spacing between frequency spikes
thr = 20;   % threshold above the ambient standard deviation
start_freq = 130; %Hz --> Start searching at this frequency (use hardcoded values to remove the 60 Hz noise below)
plot_filter_comparison = 1;         % Plots a before/after comparison 

 
[datafilt freqremoved]= autoAC_filter (datatimes, data, thr, windw, notch_bandwidth, hardcode_bands, start_freq, plot_filter_comparison);

subplot(212)
axis ([0 800 -1e-6 4e-5]);

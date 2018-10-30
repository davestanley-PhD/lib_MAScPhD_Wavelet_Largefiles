



% function [y,spk] = spikinator(x,windw,thr,amp,smth,cleanthr,frame)
%
%[y,spk] = spikinator(x,windw,thr,amp,smth,cleanthr,frame);
%
%x: spikey signal
%windw: spike window size in samples (... not time units); set equal to 1
%       to plot detection threshold before running program in standard 
%       mode        
%thresh: threshold for spike detection (greater or equal to 1); default=1.4
%amp: spike removal filter amplitude; default=2.5e-4
%    ****If amp = 0, THIS SETS PROGRAM TO EXCISION OF SPIKES***
%smth: magnitude of extra smoothing of spike stumps (0 to 1); default=0
%cleanthr: degree of spike cleaning (between 0 and 1); default=0.5
%          higher = greater cleaning, but may incur signal distortion
%frame: window frame margin (fraction of window size); default=0.2
%y: cleaned signal
%spk: time series containing removed spikes (only for non-excision mode)
%

% s = downsample (fb7s_bdet_PSPs(:,:),4);
% s = s(400000:end,:);
load s.mat

[jcs2 spk2] = spikinator(s(:,2),2000,5,0, 0, 0.5, 0.2);
[jcs3 spk3] = spikinator_dav3(s(:,2),2000,5,0, 0, 0.5, 0.2, 2500);

figure; plot (s(:,2));
hold on; plot(jcs2,'r')
hold on; plot(jcs3,'g')
hold on; plot(spk3,'k')
legend ('original signal','after excision','after excision/fade');
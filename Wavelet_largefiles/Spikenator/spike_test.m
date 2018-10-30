



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

% s = downsample (fb8s_bstochPSP_m(:,:),10);

% figure;
% s = s(10000:18000,:);

load s.mat
figure; plot (s(:,2));
% [jcs spk] = spikinator_dav2(s(:,2),2000,10,0, 0, 0.5, 0.2);
% hold on; plot(jcs,'r')
[jcs spk] = spikinator_dav3(s(:,2),2000,10,0, 0, 0.5, 0.2);
hold on; plot(jcs,'g')
hold on; plot(spk,'k')
legend ('original signal','after excision/fade','spike locations');
% File: Sulas21.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code to the fetal ECG extraction algorithm
% by Sulas et al. 2021. The original publication can be found under this 
% DOI: 10.1038/s41597-021-00811-3

% The original code was adapted to fit the benchmarking criteria, see
% publication.

function [fPeaks] = Sulas21(signal,fs)
% Fetal and maternal ECG extraction and heart rate calculation
% adapted for the NInFEA Dataset (https://physionet.org/content/ninfea/1.0.0/)
%
% By: Reza Sameni
% Email: reza.sameni@gmail.com
% May 2019
%
% Note: Use the most recent updates of the Open-Source
% Electrophysiological Toolbox (OSET) online available at: 
% https://gitlab.com/rsameni/OSET
%
% [fPeaks] = Sulas21(signal,fs)
%
% inputs:
%   signal:    [num_samples x num_channels] matrix of
%           abdominal ECG channels.
%   fs :    sampling frequency
%
% outputs:
%   fPeaks: fetal peak positions in seconds. Each marker indicates the 
%           position of one of the FQRS detected by the algorithm.

%% Algorithm parameters
f0 = 50; % notch frequency
Q = 15; % notch filter Q-factor
fm = 1.3; % approximate maternal heart rate in Hz
ITR = 4; % number of deflation iterations used for mECG cancellation
MM = 1; % number of top channels denoised by WDEN per deflation iteration
TPTR = 'rigrsure'; % WDEN threshold selection rule
SORH = 's'; % WDEN thresholding approach
SCAL = 'one'; % WDEN  multiplicative threshold rescaling
NDEN = 1; % Number of WDEN levels decrease to follow the average as it is, increase to make smoother
WNAME = 'coif5'; % WDEN mother wavelat
mpeak_detection_smoothing_len = round(fs*0.040); % maternal HR detector window length (not always used)
fpeak_detection_smoothing_len = round(fs*0.025); % fetal HR detector window length
MultiNotchNum = 1; % number of channels to remove using multichannel notch filter
wlenL = round(0.030*fs); % before R-peak
wlenR = round(0.030*fs); % after R-peak

ff = 2.4; % approximate fetal heart rate in Hz

%% Load data

% remove DC channels (if any). DC channels are later problematic, resulting in zer eigenvalues
signal = signal(std(signal, 1, 2)~=0, :);

data = signal;

T = size(data, 2);

%% Preprocessing
% preprocess the data (wide-band mode)
x_abd = signal - LPFilter(signal, 7.0/fs);
x_abd = LPFilter(x_abd, 150.0/fs);

% preprocess the data (narrow-band mode)
x = data - LPFilter(data, 15.0/fs);
x = LPFilter(x, 80.0/fs);
mref = x(1, :);


% A simple Notch filter:
Wc = f0/(fs/2);
BW = Wc/Q;
[b_notch_filter, a_notch_filter] = iirnotch(Wc, BW);
y = zeros(size(x));
for i = 1:size(x, 1)
    y(i, :) = filtfilt(b_notch_filter, a_notch_filter, x(i, :));
end
y_abd = zeros(size(x_abd));
for i = 1:size(x_abd, 1)
    y_abd(i, :) = filtfilt(b_notch_filter, a_notch_filter, x_abd(i, :));
end


%% Maternal R-peak detection

[mHR, mpeaks] = HRCalculation2(mref, ff, fs, 1, mpeak_detection_smoothing_len, 'trmean');

I_mpeaks = find(mpeaks);
mHR = 60*fs./diff(I_mpeaks);

%% mECG cancellation by deflation
[tt0, tt1] = SynchPhaseTimes2(mpeaks);
yy = PeriodicDeflDecompositionWDEN(y, ITR, MM, tt0, tt1, TPTR, SORH, SCAL, NDEN, WNAME);

% BSS
% JADE
W = jadeR(yy);
s = W*yy;

%% Fetal ECG detection
% A two-run fetal ECG detector
for run = 1:2 % repeat twice to make the fetal heart rate guess more accurate
    % ind = 1./ChannelIndex9(s, ff, fs); % Channel selection based on variance around average beat
    ind = ChannelIndex10(s, ff, fs); % Channel selection based on fixed template matched filter
    [sorted_indexes, II] = sort(ind, 1, 'ascend');
    s_sorted = s(II, :);
    fref0 = s_sorted(1, :);
    % fpeaks = PeakDetection(fref, ff/fs);
    [fhr0, fpeaks0] = HRCalculation2(fref0, ff, fs, 1, fpeak_detection_smoothing_len, 'trmean');
    I_fpeaks0 = find(fpeaks0);
    ff = fs/median(diff(I_fpeaks0));
end
% disp([D(subject).name ', ff = ' num2str(ff)]);

% fetal ECG template extraction
segment_width = wlenL + wlenR + 1;
X = zeros(length(I_fpeaks0), segment_width);
for k = 1 : length(I_fpeaks0)
    start = max(I_fpeaks0(k)-wlenL, 1);
    stop = min(I_fpeaks0(k)+wlenR, T);
    xx = [zeros(1, start - (I_fpeaks0(k)-wlenL)), fref0(start : stop), zeros(1, I_fpeaks0(k)+wlenR - T)];
    X(k, :) = xx;% - fref0(I_fpeaks0(k));
end
avg_fECG_beat = RWAverage(reshape(X, length(I_fpeaks0), segment_width));
fpeaks0 = PeakDetection4(fref0, fs, avg_fECG_beat, ff);
I_fpeaks0 = find(fpeaks0);

% PiCA using fetal R-peaks
[s, W, A] = PiCA(y_abd, fpeaks0);
fref1 = s(1, :);
c1 = corrcoef(fref0, fref1);
fref1 = sign(c1(1,2))*fref1;


% fetal ECG template extraction
segment_width = wlenL + wlenR + 1;
X = zeros(length(I_fpeaks0), segment_width);
t_segment = (0 : segment_width-1)/fs;
for k = 1 : length(I_fpeaks0)
    start = max(I_fpeaks0(k)-wlenL, 1);
    stop = min(I_fpeaks0(k)+wlenR, T);
    xx = [zeros(1, start - (I_fpeaks0(k)-wlenL)), fref1(start : stop), zeros(1, I_fpeaks0(k)+wlenR - T)];
    X(k, :) = xx;
end

avg_fECG_beat = RWAverage(reshape(X, length(I_fpeaks0), segment_width));
fpeaks1 = PeakDetection4(fref1, fs, avg_fECG_beat, ff);
fPeaks = find(fpeaks1);

%% HR refinement
wlen = 2; % # of beats before and after the current beat used for smoothing
th = 10; % HR error detection threshold in BPM
fHR0 = 60*fs./diff(I_fpeaks0);
fHR1 = 60*fs./diff(fPeaks);
fHR1_med = zeros(1, length(fHR1));
fHR1smoothed = fHR1;
hrlen = length(fHR1);
for i = 1 : hrlen
    index = max(i-wlen, 1) : min(i+wlen, hrlen);
    fHR1_med(i) = median(fHR1(index));
    if(abs(fHR1(i) - fHR1_med(i)) > th)
        k = find(I_fpeaks0 >= fPeaks(i), 1);
        if(~isempty(k) && k <= length(fHR0))
            fHR1smoothed(i) = fHR0(k);
        else
            fHR1smoothed(i) = fHR1_med(i);
        end
    end
    

end
end


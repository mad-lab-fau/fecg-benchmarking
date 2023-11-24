% File: PowerMF.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code to Power-MF algorithm, A Power Density
% Estimation and Matched Filter Based Fetal ECG Extraction Algorithm. The 
% algorithm was published in the publication with the DOI: xxxxxxxxx

% The Power-MF algorithm is a built on the algorithm 
% Varanini et al. 2014. Oritinal publication DOI: 
% 10.1088/0967-3334/35/8/1607

% The corresponding source code can be found here: 
% https://archive.physionet.org/challenge/2013/sources/


function [fPeaks] = PowerMF(ECG, fs, ms)
% Algorithm for Power-MF fetal ECG extraction algorithm
%
% [fPeaks] = PowerMF(ECG, fs, ms)
%
% inputs:
%   ECG:    [num_samples x num_channels] matrix of
%           abdominal ECG channels.
%   fs :    sampling frequency
%   ms :    minimum peak distance in milliseconds

% outputs:
%   fPeaks: fetal peak positions in seconds. Each marker indicates the 
%           position of one of the FQRS detected by the algorithm.

addpath(genpath('.'));
if nargin < 3
  ms = 340;
end

%% Preprocessing and maternal ECG cancellation according to 
%  Varanini et al. 2014

try

    cName='';
    qrsAf=[];
    dbFlag=0;                   % debug flag
    graph=0;                    % enable/disable graphical representation
    saveFig=0;                  % =1 => save figures of the processing phases


    % ---- check size of ECG ----
    if size(ECG,2)>size(ECG,1)
        ECG = ECG';
    end


    % ---- Artifact canceling ----
    % X=FecgFecgImpArtCanc(ECG,fs,cName,graph,dbFlag);
    X=FecgImpArtCanc(ECG,fs,cName,0,0);

    % ---- detrending  ----
    % Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag);
    Xd=FecgDetrFilt(X,fs,cName,0,0);

    % ---- Power line interference removal by notch filtering ----
    % Xf=FecgNotchFilt(Xd,fs,cName,graph,dbFlag);
    Xf=FecgNotchFilt(Xd,fs,cName,0,0);

    % ---- Independent Component Analysis ----
    % Xm=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig);
    Se=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig);

    % ---- Signal interpolation
    % Xi=FecgInterp(X,fs,interpFact,cName,graph);
    [Se,fs]=FecgInterp(Se,fs,4,cName,0);

    % ---- Channel selection and Mother QRS detection
    qrsM=FecgQRSmDet(Se,fs,cName,graph,dbFlag,saveFig,qrsAf);

    % ---- Mother QRS cancelling
    Xr=FecgQRSmCanc(Se,qrsM,fs,cName,graph,dbFlag,saveFig,qrsAf);

    % ---- Source separation by ICA on residual signals
    Ser=FecgICAf(Xr,fs,cName,graph,dbFlag,saveFig);
    sig = Ser';

    %% Channel Selection using Power Spectral Density (PSD)
    [~, ns]=size(Ser);

    for is=1:ns
        Ser(:,is)= (Ser(:,is)-mean(Ser(:,is)))/std(Ser(:,is));
    end
    
    % Raw derivative filter coefficients
    nu=ceil(0.005*fs); nz=floor(0.0030*fs/2)*2+1;  
    % nz = nearest odd value
    B=[ones(nu,1);zeros(nz,1);-ones(nu,1)];
    delay=floor(length(B)/2);
    
    % Compute the absolute derivative signal
    ecgfx=[repmat(Ser(1,:),delay,1);Ser;repmat(Ser(end,:),delay,1)];
    decgr=filter(B,1,ecgfx);   
    decgr= decgr(2*delay+1:end,:);
    adecg=abs(decgr);
    
    % Butterworth forward and backward bandpass filtered (1-6Hz)
    fmind=0.7; fmaxd=8;
    Wn = [fmind, fmaxd]/(fs/2);  % normalized cut-off frequency (0,1)
    [b,a]= butter(1,Wn);
    abs_dev=filtfilt(b,a,adecg);
    abs_dev = abs_dev';

    % PSD for each channel 
    for ch = 1:size(abs_dev,1)
        current = abs_dev(ch,:);
        nsc = uint16(fs *15*(60/110)); % window size 
        Nfft = max(256,2^nextpow2(nsc)); 
        if Nfft > length(current)/2 
            Nfft = length(current)/2;
        end
        [Pxx,f] = pwelch(current-mean(current),gausswin(Nfft), Nfft/2, ...
        Nfft, fs);

        % range of expected fetal heart rate from 1.8 Hz to 3.0 Hz
        min_freq = find(f>1.8);
        max_freq = find(f<3);
        min_freq = min_freq(1);
        max_freq = max_freq(end);

        wert = max(findpeaks(Pxx(min_freq:max_freq))); % finds peaks in 
        % fetal heart rate range
        if isempty(wert) % no peak in range
            max_peak(ch) = 0;
            continue
        end
        freq = f(find(ismember(Pxx,wert)));
        disp("channel " + ch + ": " + "freq. " + freq + "/" + freq*60);
        max_peak(ch) = wert; % largest peak per channel
    end

    [~,channel] = max(max_peak); % channel with maximum peak is selected 

    %% Fetal R-peak detection using matched filter
    
    % Preliminary peak detection
    distance = ms / 1000 * fs;
    [~,peaks] = findpeaks(abs_dev(channel,:),'MINPEAKDISTANCE',distance); 

    % template size: median RR interval 
    med_size= round(median(diff(peaks))/2); 

    % template generation
    for j = 3:length(peaks)-2 % exclude peaks in the beginning and end to 
        % avoid reaching end of array
        % template: from (peak - median/2) to (peak + median/2)
        template(j,:)=sig(channel,peaks(j)-med_size:peaks(j)+med_size);
    end
    template_med = median(template,1);
    r = matched_filter(sig(channel,:),template_med); % apply matched filter
    [~,fPeaks] = findpeaks(r,'MinPeakDistance',distance); % find peaks
    fPeaks = fPeaks / 4;

catch e
    fprintf(1,'There was an error! The message was:\n%s',e.message);
    fPeaks    = [];
end
end

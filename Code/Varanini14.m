% File: Varanini14.m
%
% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023
%
% This file contains the source code to the fetal ECG extraction algorithm
% by Varanini et al. 2014. The original publication can be found under this 
% DOI: 10.1088/0967-3334/35/8/1607

% The original code was adapted to fit the benchmarking criteria, see
% publication.


function [fetal_QRSAnn_est,Ser,ecgs] = Varanini14(ECG,fs,cName,qrsAf,graphDflags)
% Algorithm for Physionet/CinC competition 2013.
%
% [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG,cName,qrsAf)
%
% inputs :
%   ECG :   60000x4 (4 channels and 1min of signal at 1000Hz) matrix of
%           abdominal ECG channels.
%   Fs :    Sampling frequency
%   cName : record name (optional)
%   qrsAf : QRS markers in seconds (optional, learning set only)

% outputs :
%   fetal_QRSAnn_est :  FQRS markers. Each marker indicates the position of one
%                       of the FQRS detected by the algorithm.
%   Ser : Signal after mECG suppression and fetal ECG enhancement
%   ecgs : Derivative signal, signal on which peaks are detected
%

%% Parameter definition
addpath(genpath('.'));
Fs = fs;

try
    if(nargin<3), cName=''; end
    if(nargin<4), qrsAf=[]; end
    if(nargin<5)
        dbFlag=0;                   % debug flag
        graph=0;                    % enable/disable graphical representation
        saveFig=0;                  % =1 => save figures of the processing phases
        saveFigRRf=0;               % =1 => save estimated fetal RR figures
    else
        dbFlag=graphDflags.dbFlag;
        graph=graphDflags.graph;
        saveFig=graphDflags.saveFig;
        saveFigRRf=graphDflags.saveFigRRf;
    end
           
    % ---- check size of ECG ----
    if size(ECG,2)>size(ECG,1)
        ECG = ECG';
    end
        
    %% Preprocessing   
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
    
    %% Maternal ECG cancellation
    % ---- Channel selection and Mother QRS detection
    qrsM=FecgQRSmDet(Se,fs,cName,graph,dbFlag,saveFig,qrsAf);
    
    % ---- Mother QRS cancelling
    Xr=FecgQRSmCanc(Se,qrsM,fs,cName,graph,dbFlag,saveFig,qrsAf);


    % ---- Source separation by ICA on residual signals
    Ser=FecgICAf(Xr,fs,cName,graph,dbFlag,saveFig);

    %% Fetal peak detection
    % ---- Channel selection and Fetal QRS detection
    % qrsF=FecgQRSfDniAdf(Ser,fs,cName,qrsM,graph,dbFlag,saveFig,saveFigRRf,qrsAf);
    [qrsF,ics, ecgs]=FecgQRSfDet(Ser,fs,cName,qrsM,graph,dbFlag,saveFig,saveFigRRf,qrsAf);
    
    fetal_QRSAnn_est=qrsF.'*Fs;
    QT_Interval         = [];
    
catch e
    fprintf(1,'There was an error! The message was:\n%s',e.message);
    fetal_QRSAnn_est    = [];
    QT_Interval         = [];
end

end %== function ================================================================
%

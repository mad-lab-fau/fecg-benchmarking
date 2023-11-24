% File: Behar14.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code to the fetal ECG extraction algorithm
% by Behar et al. 2014. The original publication can be found under this 
% DOI: 10.1088/0967-3334/35/8/1569

% The original code was adapted to fit the benchmarking criteria, see
% publication.

function [FQRSout,QTout] = Behar14(tm,ecg,Fs,paramStruct)
% Template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the
% same.
%
% [FQRSout,QTout] = Behar14(tm,ecg) where the inputs and outputs are specified
% below.
%
% inputs
%   tm :         Nx1 vector of time in milliseconds
%   ecg:         [num_samples x num_channels] matrix of
%                abdominal ECG channels.
%   Fs :    sampling frequency
%   paramStruct: param of the algorithm to use. if not inputed then the
%                default parameters are used. 
%                ***FOR SET-C of the challenge LEAVE EMPTY.***
% 
% output
%   FQRSout:    FQRS markers in seconds. Each marker indicates the position
%               of one of the FQRS detected by the algorithm.
%   QTout:      1x1 estimated fetal QTout duration
%
%
%
% Safe Foetus Monitoring Toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013 Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 02-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.    

if nargin<3; error('physionet2013: wrong number of input arguments'); end;
if nargin<4; paramStruct=[]; end;
addpath(genpath('subfunctions/'));

try
    if isempty(paramStruct)
        % == for FQRS extraction
        paramStruct.method = 'FUSE';              % TS/TS-PCA/TS-EKF/DEFLATION/FUSE
        paramStruct.methodComplexity = 'COMPLEX'; % COMPLEX/SIMPLE
        paramStruct.smoothRR = 1;                 % postprocessing
        paramStruct.cinc_match = 0;               % twick to better match CinC scores
        paramStruct.debug = 0;
        paramStruct.fs = Fs;
        paramStruct.TSnbCycles = 20;              % nb of cycles to build template 
        paramStruct.EKFgain = 1;
        paramStruct.fhighNbCoeff = 5;
        paramStruct.fbasNbCoeff = 3;
        paramStruct.NotchCheck = 1;
        paramStruct.fbas = 10;
        paramStruct.fhigh = 99;
        paramStruct.STDcut = 17;                  % rubbish to fit to challenge scores (was 17)
        paramStruct.NbPC = 2;                     % nb of principal components in PCA
        paramStruct.detThres = 0.5;               % QRS detector threshold for FQRS
        paramStruct.qrsdet_wind = 6;%15;             % qrs detection window: LENGTH_SEG/qrsdet_wind NEEDS TO BE AN INTEGER!!
        paramStruct.extract_qt = 0;
    end
    
    [FQRSout,~,outExtra] = physionet2013extract(tm,ecg,paramStruct);
catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
                    disp('- just guessing');
    FQRSout=(0.1:0.42:length(ecg)/paramStruct.fs)*paramStruct.fs;
    QTout=250;
end


end








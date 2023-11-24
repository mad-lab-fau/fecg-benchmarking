% File: benchmark_ninfea.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file computes performance metrics based on ground truth fetal V-peak
% positions from PWD traces and alorithmically computed R-peak positions. 
% For more information see original publication DOI:
% 10.1038/s41597-021-00811-3.


function[F1,PPV,SE,TP,FN,FP]= benchmark_ninfea(labels,results)
% Algorithm benchmarking function, adapted for NInFEA dataset.

% [F1,a,PPV,SE,TP,FN,FP]= benchmark_ninfea(labels,results)
%
% inputs:
%   labels  :  Ground truth fetal peak location from PWD traces, in samples
%   results :  Algorithmically computed fetal peak locations in samples.  
%
% outputs:
%   F1  : F1-score
%   PPV : Positive-predictive value
%   SE  : Sensitivity
%   TP  : True positives, peaks that are closest to label
%   FN  : False negatives, not found peaks
%   FP  : False positives, found peaks, that are not true positives

%% Parameter definition
TP = [];
FN = [];
acceptint = 200*2048/1000; % 200 ms acceptance interval

%% Metric calculation
for i = 1:length(labels)
    a = intersect(find(results < labels(i)),find(results>(labels(i)-acceptint))); % find qrs 200ms before label
    if ~isempty(a)
        if size(a,2) > 1
            TP(end+1)=a(end); % True Positive: qrs closest to label
        else
            TP(end+1)=a;
        end
    else
        FN(end+1) = i; % False Negative: no qrs found
    end
end

results_ind = 1:length(results);
FP = setdiff(results_ind,TP); % False Positive: qrs that are not True Positives

FN = length(FN);
FP = length(FP);

% compute Se, PPV, F1
if (FN + TP) == 0
    SE = 0;
else
    SE  = TP/(TP+FN);
end
if(FP+TP) == 0
    PPV = 0;
else
    PPV = TP/(FP+TP);
end
if (SE+PPV) == 0
    F1 = 0;
else
F1 = 2*SE*PPV/(SE+PPV);
end
end
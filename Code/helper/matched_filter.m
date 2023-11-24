% File: matched_filter.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code of the matched filter of the Power-MF 
% fetal ECG extraction algorithm corresponding to the publication with the 
% DOI: xxxxxxxxx

function [r] = matched_filter(signal,template)
%
% r = matched_filter(signal, template),
% Fetal QRS detector based on Matched Filter
%
% inputs:
%   signal: vector of input data
%   template: template waveform
%
% outputs:
%   r: filtered signal, has peaks at location where fetal QRS is expected


N = length(signal);
L = length(template);
w = floor(L/2);

template = template(end:-1:1);

r = filter(template,1,[signal zeros(1,w-1)]);

r = r(w:N+w-1);

% File: benchmark_algorithms.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code for the paramerer optimization of the 
% Power-MF algorithm corresponding to the publication with the DOI: xxxxxxxxx

% Optimized parameter is minimum peak distance. For the minimum peak 
% distance optimization, physiologically meaningful values from 290 ms to 
% 360 ms in 10 ms steps were set. For each value, the local maxima were
% computed on all recordings of CinC2013 and compared with the ground truth 
% fetal R-peak annotations. The minimum peak distance is set to the value 
% leading to the highest F1.

clear; close all;
addpath(genpath('..'));

% check if subfunctions are installed
foldercontent = dir('./subfunctions/CinC_Behar');
if numel(foldercontent) < 5
    warning('You need to download the code for Behar algorithm here: https://archive.physionet.org/challenge/2013/sources')
    return
end
foldercontent = dir('./subfunctions/CinC_Varanini');
if numel(foldercontent) < 5
    warning('You need to download the code for Varanini algorithm here: https://archive.physionet.org/challenge/2013/sources')
    return
end
foldercontent = dir('./subfunctions/fecgsyn');
if numel(foldercontent) < 5
    warning('You need to download the fecgsyn toolbox here: https://github.com/fernandoandreotti/fecgsyn')
    return
end
foldercontent = dir('./subfunctions/OSET-master');
if numel(foldercontent) < 5
    warning('You need to download the OSET toolbox for Sulas algorithm here: https://github.com/alphanumericslab/OSET')
    return
end

% check if Challenge2013 dataset is downloaded
path = fullfile('../Data/','Challenge2013','SetA');
if ~exist(path, 'dir')
    warning('You need to download the Challenge2013 dataset to the folder ./Data/Challenge2013. Download link: https://physionet.org/content/challenge-2013/1.0.0/')
    return
end

%% Parameter definitions
% Define which datasets / algorithms to benchmark
datasets = "challenge";
algorithms="PowerMF";

pth = pwd;
window = [260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360];


for mss = 1:length(window)
    for i = 1:length(datasets)
        for alg = 1:length(algorithms)
            
            %% Load and prepare data
            ms = window(mss);
            dataset = datasets(i);
            algorithm = algorithms(alg);
            
            path = '../Data/Challenge2013/SetA/';

            
            D = dir([path '/*.mat']);
            score=[];
            for record = 1:length(D)
                load([path,D(record).name]);
                acceptint = 50*Fs/1000; % 50 ms acceptance interval interval
                
                % Delete NaN values from signal
                idx = find(isnan(signal));
                for a = 1:length(idx)
                    signal(idx(a)) = mean([signal(idx(a)-1) signal(idx(a)+1)]);
                end
                
                %% Compute performance of algorithm
                % Signals with length>500000 are split up
                % (there were problems for some algorithms)
                if length(signal)>500000
                    bis = floor(length(signal)/2);
     
                    % Compute F1 of PowerMF with parameter set
                    results1 = PowerMF(signal(:,1:bis),Fs, ms);
                    results2 = PowerMF(signal(:,bis+1:end),Fs, ms);
      
                    results = [results1 results2+bis+1];
                    cd(pth)
                    
                else
                    % Compute F1 of PowerMF with parameter set
                    results = PowerMF(signal,Fs,ms);
                end
                
                %% Compute performance metrics
                
                [F1,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(fqrs.', results, acceptint);
                score(record).F1 = F1;
                score(record).ACC = TP / (TP + FP + FN);
                score(record).MAE = MAE;
                score(record).PPV = PPV;
                score(record).SE = SE;
                score(record).TP = TP;
                score(record).FN = FN;
                score(record).FP = FP;
                disp("done: " + dataset + " record " + record + ". " + algorithm+ ": F1=" + F1);
                clearvars results results1 results2 signal_n
                cd(pth)
            end
            %% Save results
            name = strcat(dataset,"_",algorithm,'_',string(ms),'ms.mat');
            save_path = '../Results/paramoptimization/'+name;
            save(save_path,'score');
            clearvars score
        end
    end
end

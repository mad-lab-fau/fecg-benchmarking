% File: benchmark_algorithms.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023

% This file contains the source code to the fetal ECG algorithm
% benchmarking corresponding to the publication with the DOI: xxxxxxxxx

% The benchmarking includes the algorithms by Varanini et al. 2014, Behar
% et al. 2014 and Sulas et al. 2021. The corresponding source code is 
% open-source available:

% https://archive.physionet.org/challenge/2013/sources/
% https://sameni.org/OSET/
% https://github.com/rsameni/NInFEADataset

% The benchmarking is performed on the NInFEA and ADFECG datasets which are
% available in online repositories:

% CinC2013 (https://physionet.org/content/challenge-2013/1.0.0/)
% NInFEA (https://doi.org/10.13026/c4n5-3b04)
% ADFECG (https://doi.org/10.6084/m9.figshare.c.4740794.v1)

clear; close all;
addpath(genpath('.'));
datapath = '../Data';


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

%% Parameter definitions
% define which datasets / algorithms to benchmark
datasets = ["b1" "b2" "ninfea" "challenge"];

algorithms = ["powermf" "varanini" "sulas" "behar"];
configs = ["a" "b" "c" "d" "e" "all"];

pth = pwd;
for i = 1:length(datasets)
    for alg = 1:length(algorithms)
        for c = 1:length(configs)
            config = configs(c);
            
            dataset = datasets(i);
            algorithm = algorithms(alg);
            
            % load signals + labels
              switch dataset
                case "ninfea"
                    path = fullfile(datapath,'NInFEA','all_channels');
                    if ~exist(path, 'dir')
                        warning('You need to download the NInFEA dataset to the folder ./Data/NInFEA. Download link: https://physionet.org/content/ninfea/1.0.0/')
                        return
                    end
                        
                case "b1"
                    path = fullfile(datapath,'ADFECGDB','B1_Pregnancy_dataset');
                    if ~exist(path, 'dir')
                        warning('You need to download the ADFECG dataset to the folder ./Data/ADFECG. Download link: https://springernature.figshare.com/collections/Fetal_electrocardiograms_direct_and_abdominal_with_reference_heart_beats_annotations/4740794/1')
                        return
                    end
                case "b2"
                    path = fullfile(datapath,'ADFECGDB','B2_Labour_dataset');
                    if ~exist(path, 'dir')
                        warning('You need to download the ADFECG dataset to the folder ./Data/ADFECG. Download link: https://springernature.figshare.com/collections/Fetal_electrocardiograms_direct_and_abdominal_with_reference_heart_beats_annotations/4740794/1')
                        return
                    end
                case "challenge"
                    path = fullfile(datapath,'Challenge2013','SetA');
                    if ~exist(path, 'dir')
                        warning('You need to download the Challenge2013 dataset to the folder ./Data/Challenge2013. Download link: https://physionet.org/content/challenge-2013/1.0.0/')
                        return
                    end
               end
            
            D = dir(fullfile(path, '*.mat'));
            score=[];
            timeconsumption=[];
            for record = 1:length(D)              
                load(fullfile(path,D(record).name));
                               
                acceptint = 50*Fs/1000; % 50 ms acceptance interval interval
                
                if dataset == "b1" || dataset == "b2"
                    signal = signal(1:4,:);
                elseif dataset == "ninfea"
                    switch config
                        case "a"
                            signal_n(1,:) = signal(5,:);
                            signal_n(2,:) = signal(8,:);
                            signal_n(3,:) = signal(12,:);
                            signal_n(4,:) = signal(15,:);
                            signal_n(5,:) = signal(18,:);
                            signal = signal_n;
                        case "b"
                            signal_n(1,:) = signal(5,:);
                            signal_n(2,:) = signal(12,:);
                            signal_n(3,:) = signal(18,:);
                            signal_n(4,:) = signal(21,:);
                            signal = signal_n;
                        case "c"
                            signal_n(1,:) = signal(2,:);
                            signal_n(2,:) = signal(3,:);
                            signal_n(3,:) = signal(14,:);
                            signal_n(4,:) = signal(15,:);
                            signal = signal_n;
                        case "d"
                            signal_n(1,:) = signal(1,:);
                            signal_n(2,:) = signal(2,:);
                            signal_n(3,:) = signal(4,:);
                            signal_n(4,:) = signal(18,:);
                            signal = signal_n;
                        case "e"
                            signal_n(1,:) = signal(1,:);
                            signal_n(2,:) = signal(2,:);
                            signal_n(3,:) = signal(4,:);
                            signal_n(4,:) = signal(17,:);
                            signal_n(5,:) = signal(18,:);
                            signal = signal_n;
                        case "all"
                            signal = signal(1:22,:);
                            if record == 34
                                signal(9,:) = [];
                            end
      
                    end
                    
                    acceptint = 200*Fs/1000;                    
                end
                
                % delete NaN values from signal
                idx = find(isnan(signal));
                for a = 1:length(idx)
                    signal(idx(a)) = mean([signal(idx(a)-1) signal(idx(a)+1)]);
                end
                
                % signals with length>500000 are split up
                % (there were problems for some algorithms)
                if length(signal)>500000
                    bis = floor(length(signal)/2);
                    signallength = floor(length(signal)/Fs);
                    switch algorithm
                        case "varanini"
                            results1 = Varanini14(signal(:,1:bis), Fs);
                            results2 = Varanini14(signal(:,bis+1:end), Fs);
                            
                        case "sulas"
                            results1 = Sulas20(signal(:,1:bis), Fs);
                            results2 = Sulas20(signal(:,bis+1:end), Fs);
                        case "behar"
                            results1 = Behar14(tm, signal(:,1:bis), Fs);
                            results2 = Behar14(tm, signal(:,bis+1:end), Fs);
                        case "powermf"
                            results1 = PowerMF(signal(:,1:bis),Fs);
                            results2 = PowerMF(signal(:,bis+1:end),Fs);
                                                      
                    end
                    results = [results1 results2+bis+1];
                    cd(pth)
                    
                else
                    signallength = floor(length(signal)/Fs);
                    switch algorithm
                                                                      
                        case "powermf" 
                            results = PowerMF(signal,Fs);
                        case "varanini"
                            results = Varanini14(signal, Fs);
                        case "sulas"
                            results = Sulas20(signal, Fs);
                        case "behar"
                            results = Behar14(tm, signal, Fs);
                    end
                    cd(pth)
                end
                
                %% benchmark
                if dataset == "ninfea"
                    [F1,PPV,SE,TP,FN,FP] = benchmark_ninfea(fqrs.', results);
                    if isempty(F1)
                        score(record).F1 = 0;
                        score(record).ACC = 0;
                        score(record).PPV = 0;
                        score(record).SE = 0;
                        score(record).TP = 0;
                        score(record).FN = 0;
                        score(record).FP = 0;
                    else                      
                        score(record).F1 = F1;
                        score(record).ACC = TP / (TP + FP + FN);
                        score(record).PPV = PPV;
                        score(record).SE = SE;
                        score(record).TP = TP;
                        score(record).FN = FN;
                        score(record).FP = FP;
                    end
                    disp("done: " + dataset + " record " + record + ". " + algorithm+ ": F1=" + F1);
                else
                    [F1,MAE,PPV,SE,TP,FN,FP] = Bxb_compare(fqrs.', results, acceptint);
                    if isempty(F1)
                        score(record).F1 = 0;
                        score(record).ACC = 0;
                        score(record).MAE = 0;
                        score(record).PPV = 0;
                        score(record).SE = 0;
                        score(record).TP = 0;
                        score(record).FN = 0;
                        score(record).FP = 0;
                    else
                        score(record).F1 = F1;
                        score(record).ACC = TP / (TP + FP + FN);
                        score(record).MAE = MAE;
                        score(record).PPV = PPV;
                        score(record).SE = SE;
                        score(record).TP = TP;
                        score(record).FN = FN;
                        score(record).FP = FP;
                    end
                    
                    disp("done: " + dataset + " record " + record + ". " + algorithm+ ": F1=" + F1);
                end
                clearvars results results1 results2 signal_n
                cd(pth)
            end
            % save results
            if dataset == "ninfea"
                dataset = strcat(dataset,"_",config);
            end
            name = fullfile('../Results/', strcat(dataset,"_", algorithm,'.mat'));
            save(name,'score');
            clearvars score
           
        end
    end
end
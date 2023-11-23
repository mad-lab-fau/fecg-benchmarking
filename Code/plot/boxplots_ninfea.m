% File: boxplots_ninfea.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023
%
% This file is a script to create boxplots of the algorithm benchmarking
% for the NInFEA dataset
%

close all
addpath(genpath('..'));

path = '../../Results/';

Configs = {'a)' 'b)' 'c)' 'd)' 'e)' 'All abdominal electrodes'};

Colors = {[0.255 0.451 0.553], [0.851 0.491 0.227], [0.663 0.761 0.816], [0 0.184 0.424]};
Colors_light = {[0.643 0.753 0.804], [0.949 0.765 0.627], [0.788 0.835 0.855], [0.514 0.604 0.741]};

legendEntries = {'Power-MF' 'Varanini et al. 2014' 'Behar et al. 2014' 'Sulas et al. 2021'};

data = ["PowerMF";"Varanini";"Behar";"Sulas"];
for i=1:length(data)
    dataset = data(i);
    switch dataset
        case 'PowerMF'
            load(fullfile(path, 'ninfea_a_powermf.mat'))
            A = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_b_powermf.mat'))
            B = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_c_powermf.mat'))
            C = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_d_powermf.mat'))
            D = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_e_powermf.mat'))
            E = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_all_powermf.mat'))
            all = [score.F1]' * 100;
            PowerMF = [A, B, C, D, E, all];
        case 'Varanini'
            load(fullfile(path,  'ninfea_a_varanini.mat'))
            A = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_b_varanini.mat'))
            B = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_c_varanini.mat'))
            C = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_d_varanini.mat'))
            D = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_e_varanini.mat'))
            E = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_all_varanini.mat'))
            all = [score.F1]' * 100;
            Varanini = [A, B, C, D, E, all];
        case 'Sulas'
            load(fullfile(path,  'ninfea_a_sulas.mat'))
            A = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_b_sulas.mat'))
            B = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_c_sulas.mat'))
            C = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_d_sulas.mat'))
            D = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_e_sulas.mat'))
            E = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_all_sulas.mat'))
            all = [score.F1]' * 100;
            Sulas = [A, B, C, D, E, all];
        case 'Behar'
            load(fullfile(path,  'ninfea_a_behar.mat'))
            A = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_b_behar.mat'))
            B = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_c_behar.mat'))
            C = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_d_behar.mat'))
            D = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_e_behar.mat'))
            E = [score.F1]' * 100;
            load(fullfile(path,  'ninfea_all_behar.mat'))
            all = [score.F1]' * 100;
            Behar = [A, B, C, D, E, all];
    end
end
GroupedData = {PowerMF Varanini Behar Sulas};

N = numel(GroupedData);
delta = linspace(-.1,.4,N); %// define offsets to distinguish plots
width = .15; %// small width to avoid overlap
cmap = hsv(N); %// colormap
legWidth = 0; %// make room for legend

figure;
hold on;

fsz = 13;
set(gca,'FontSize',fsz)

for ii=1:N
    labels = Configs;
    boxplot(GroupedData{ii},'Color', Colors{ii}, 'boxstyle','filled', ...
        'position',(1:numel(labels))+delta(ii), 'widths',width, 'labels',labels, 'OutlierSize',4, 'Symbol', 'k+')
    plot(NaN,1,'color',Colors{ii}); %// dummy plot for legend
end
xlabel('Electrode configuration', 'FontSize', fsz+4); ylabel('F1 in %', 'FontSize', fsz+4); grid on;
xlim([1+2*delta(1) numel(labels)+legWidth+2*delta(N)]) %// adjust x limits, with room for legend

[~, hobj, ~, ~] = legend(legendEntries,'Location', 'northoutside', 'Orientation', 'horizontal','FontSize', fsz+4);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',10);

for i=2:2:8
    a=get(get(gca,'children'),'children');
    t=get(a{i},'tag');
    idx=strcmpi(t,'box');
    boxes=a{i}(idx);
    set(boxes,'linewidth',10);
    idx=strcmpi(t,'median');
    median=a{i}(idx);
    set(median,'linewidth',2);
end

for ii=1:N
    plot((1:numel(labels))+delta(ii), mean(GroupedData{ii}, 'omitnan'), 'diamond', 'MarkerFaceColor',Colors_light{ii}, 'MarkerEdgeColor','black')
end

x=gcf;
pos = x.Position;
sizenew = [1200 500];
x.Position = [50 50 sizenew];
%exportgraphics(gcf,'../../Plots/boxplots_ninfea.pdf')

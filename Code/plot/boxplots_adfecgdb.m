% File: boxplots_adfecgdb.m

% Author: Katharina M. Jaeger, katharina.jaeger@fau.de
% Created: January 2023
%
% This file is a script to create boxplots of the algorithm benchmarking
% for the ADFECG dataset
%

close all
addpath(genpath('..'));

path = '../../Results/';

Configs = {'B1' 'B2'};

Colors = {[0.255 0.451 0.553], [0.851 0.491 0.227], [0.663 0.761 0.816], [0 0.184 0.424]};
Colors_light = {[0.643 0.753 0.804], [0.949 0.765 0.627], [0.788 0.835 0.855], [0.514 0.604 0.741]};

legendEntries = {'Power-MF' 'Varanini et al. 2014' 'Behar et al. 2014' 'Sulas et al. 2021'};

data = ["PowerMF";"Varanini";"Behar";"Sulas"];
for i=1:length(data)
    dataset = data(i);
    switch dataset
        case 'PowerMF'
            load(fullfile(path, 'b1_powermf.mat'))
            A = [score.F1]' * 100;
            A(11)=NaN;
            A(12)=NaN;
            load(fullfile(path, 'b2_powermf.mat'))
            B = [score.F1]' * 100;
            PowerMF = [A, B];
        case 'Varanini'
            load(fullfile(path, 'b1_varanini.mat'))
            A = [score.F1]' * 100;
            A(11)=NaN;
            A(12)=NaN;
            load(fullfile(path, 'b2_varanini.mat'))
            B = [score.F1]' * 100;
            Varanini = [A, B];
        case 'Sulas'
            load(fullfile(path, 'b1_sulas.mat'))
            A = [score.F1]' * 100;
            A(11)=NaN;
            A(12)=NaN;
            load(fullfile(path, 'b2_sulas.mat'))
            B = [score.F1]' * 100;
            Sulas = [A, B];
        case 'Behar'
            load(fullfile(path, 'b1_behar.mat'))
            A = [score.F1]' * 100;
            A(11)=NaN;
            A(12)=NaN;
            load(fullfile(path, 'b2_behar.mat'))
            B = [score.F1]' * 100;
            Behar = [A, B];
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

fsz = 10;
set(gca,'FontSize',fsz)

for ii=1:N
    labels = Configs;
    boxplot(GroupedData{ii},'Color', Colors{ii}, 'boxstyle','filled', ...
        'position',(1:numel(labels))+delta(ii), 'widths',width, 'labels',labels, 'OutlierSize',4, 'Symbol', 'k+')
    plot(NaN,1,'color',Colors{ii}); %// dummy plot for legend
end



xlabel('Subset','FontSize', fsz+4); ylabel('F1 in %','FontSize', fsz+4); grid on;
xlim([1+2*delta(1) numel(labels)+legWidth+2*delta(N)]) %// adjust x limits, with room for legend
ylim([-4 107])
[~, hobj, ~, ~] =legend(legendEntries,'Location','south','FontSize', fsz);
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
sizenew = [400 300];
x.Position = [50 50 sizenew];

%exportgraphics(gcf,'../../Plots/boxplots_adfecg.pdf')
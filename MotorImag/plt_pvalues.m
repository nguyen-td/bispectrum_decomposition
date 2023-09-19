% Script to plot p-values to determine significant interactions. Given a significance level, 
% forward vectors (of the forward mixing matrix A) corresponding to the significance
% interactions are plotted. 

%% Set paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition'))
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab')
eeglab

%% Set parameters and load data
% set alpha level
alpha = 0.05;

% load data
sub = 14;
load(sprintf('MotorImag/data/P_vp%d.mat', sub))
load(sprintf('MotorImag/data/A_vp%d.mat', sub))
load('MotorImag/data/chanlocs.mat')

%% Plotting
% plot p-values
n = size(P, 1);
tiledlayout(1, n)
for i = 1:n
    nexttile
    imagesc(-log10(P(:, :, i))); 
    c = colorbar;
    c.Label.String = '-log10(p)';
    title(sprintf('m:n:%d', i))
end

% plot significant interactions as topoplot
[x, y, z] = ind2sub(size(P), find(P < alpha));

% make topoplot
data = [A(:,x) A(:,y) A(:,z)];
for i = 1:size(data, 2)
    figure; topoplot(data(:,i), chanlocs, 'electrodes', 'on'); colorbar; clim([min(data, [], 'all') max(data, [], 'all')])
end

% manual plot
figure; topoplot(A(:,5), chanlocs, 'electrodes', 'on'); colorbar; clim([min(data, [], 'all') max(data, [], 'all')])

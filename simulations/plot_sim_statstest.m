% Plot resutlts of the main_sim_statstest.m pipeline.

%% Setup
clear
DIRIN = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/simulations/sim_stats/results/';
f_name = 'P_fdr_traintest_off';
sim_case = 1;
n_iter = 100;
alpha = 0.05;

%% Load .mat files containing FPRs
snr = [0.2 0.4 0.5 0.6 0.8];
fpr = zeros(length(snr), n_iter);

isnr = 3;
P = load([DIRIN f_name '_snr' num2str(snr(isnr)) '_case' num2str(sim_case) '.mat']);
FPR = FPR.fpr_iter;

t = randi(100, 2, 100);
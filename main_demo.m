% Demo script that shows the decomposition method on a simulated dataset without
% any complex analyses.

%% Simulate time series (one univariate PAC interaction, no bivariate interaction)
sim_case = 1;
n_univ = 1;
n_biv = 0;
isnr = 0.8;

[data, fs, source, filt, L, D] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr);

%% Compute sensor cross-bispectrum
% sampling frequency
fres = fs; 
frqs = sfreqs(fres, fs); % freqs in Hz

% set parameter values for (cross)-bispectrum estimation
freqinds = [10 10]; % in Hz
freqpairs = get_freqindices(round_to_05(freqinds(1)), round_to_05(freqinds(2)), frqs);   
len_epochs = 2; % 2-second epochs
segleng = fs * len_epochs; 
segshift = floor(segleng/2);
epleng = fs * len_epochs; % create epochs of [e.epleng] seconds

% compute cross-bispectrum
para.nrun = 1;
[~, bs_orig, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs, para); % compute surrogate and true cross-bispectra

%% Run decomposition
n = 3; % number of fitted source interactions
[A_hat, D_hat] = bsfit_freqbands(bs_orig, n); 

%% Comments
% There are different ways how to go from here. For ideas, read the thesis
% or/and have a look at the main simulation script 'simulations/main_sim_pac' 
% or the data analysis in 'main.m'.
%
% Various plotting scripts are available in the 'plotting/' folder. Again,
% have a look at the main scripts to see how they can be used to visualize
% results.
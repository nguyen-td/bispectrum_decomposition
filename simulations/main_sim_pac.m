% Pipeline to run the decomposition method on simulated PAC data with a
% specified number of univariate and bivariate interactions. 
%
% Input:
%   n_shuf - [integer] number of shuffles
%
% Optional inputs:
%   n      - [integer] model order/number of fitted sources, default is 5. Can also be an array of n's, e.g., [3 4 5]
%   n_univ - [integer] number of univariate interactions, default is 1
%   n_biv  - [integer] number of bivariate interactions, default is 2

function main_sim_pac(n_shuf)

    % setup
    eeglab
    
    % simulate one bivariate PAC interaction between an occopital and a parietal region
    roi_idx1 = 27; %lingual L
    roi_idx2 = 59; % superiorparietal L
    [signal_sensor, fs, source] = sim_wholebrain_pac(2, 0, 1, 0.5, [roi_idx1 roi_idx2]);
    psd = pwelch(signal_sensor(:,:)', 100, 50, 2*fs, fs);

    % plot sources
end
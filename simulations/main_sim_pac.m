% Pipeline to run the decomposition method on simulated PAC data with a
% specified number of univariate and bivariate interactions. 
%
% Input:
%   n_shuf - [integer] number of shuffles
%
% Optional inputs:
%   alpha  - [float] significance level, default is 0.05.

function main_sim_pac(n_shuf, varargin)
    
%     local_path = '/data/tdnguyen/git_repos/';
    local_path = '/Users/nguyentiendung/GitHub/';
    DIROUT = [local_path 'bispectrum_decomposition/simulations/sim_pac/figures/'];

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        });
    if ischar(g), error(g); end
    
    %% Simulate one univariate PAC interaction at a random region
%     [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(2, 0, 1, 0.5);

    %% Simulate one bivariate PAC interaction between an occopital and a parietal region
    roi_idx1 = 27; %lingual L
    roi_idx2 = 59; % superiorparietal L
    [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(2, 0, 1, 0.5, [roi_idx1 roi_idx2]);
    
    % sampling frequency
    fres = fs; 
    frqs = sfreqs(fres, fs); % freqs in Hz

    % set parameter values for (cross)-bispectrum estimation
    freqinds = [mean(filt.low) mean(filt.high)]; % in Hz
    len_epochs = 2; % 2-second epochs
    segleng = fs * len_epochs; 
    segshift = floor(segleng/2);
    epleng = fs * len_epochs; % create epochs of [e.epleng] seconds 

    % run decomposition without antisymmetrization
    n = [1 2 3];
    err_colors = ['r', 'b', 'k'];
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L);
    plot_error(errors, 1, n, err_colors, '', '', DIROUT)

end
% Pipeline to run the decomposition method on simulated PAC data for f1 = [9 11] Hz and f2 = [58 62] Hz with a
% specified number of univariate and bivariate interactions. 
%
% Input:
%   n_shuf   - [integer] number of shuffles
%   poolsize - [integer] number of workers in the parellel pool (check parpool documentation) for parallel computing
%
% Optional inputs:
%   alpha    - [float] significance level, default is 0.05.

function main_sim_pac(n_shuf, varargin)

    %% Setup
%     local_path = '/data/tdnguyen/git_repos/';
    local_path = '/Users/nguyentiendung/GitHub/';
    DIROUT = [local_path 'bispectrum_decomposition/simulations/sim_pac/figures/'];

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end

    % setup plotting
    load cm17
    
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

    % analyze univariate sensor cross-bispectrum
    [~, ~, P_sens_fdr_uni] = freq_preselection(signal_sensor, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    p_cmap = cmap_pvalues(P_sens_fdr_uni, cm17, cm17a);
    plot_pvalues_univ(P_sens_fdr_uni, frqs, '', p_cmap, DIROUT, 'f_ext', '.fig')
    
    % analyze normal sensor cross-bispectrum
    freqpairs = get_freqindices(round_to_05(freqinds(1)), round_to_05(freqinds(2)), frqs);   
    para.nrun = n_shuf;
    [bs_all, bs_orig, ~] = data2bs_event_surro_final(signal_sensor(:, :)', segleng, segshift, epleng, freqpairs, para); % compute surrogate and true cross-bispectra
    ichan = 3;
    [~, P_sens_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig), ichan)), squeeze(mean(abs(bs_all), ichan)), n_shuf, g.alpha);
    p_cmap = cmap_pvalues(P_sens_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_fdr, frqs, '', p_cmap, DIROUT, 'bispec_type', ['_cross_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false)

    % analyze antisymmetrized sensor cross-bispectrum
    bs_orig_anti = bs_orig - permute(bs_orig, [2, 1, 3]); % B_ijk - B_jik
    bs_all_anti = bs_all - permute(bs_all, [2, 1, 3, 4]); % B_ijk - B_jik
    [~, P_sens_anti_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig_anti), ichan)), squeeze(mean(abs(bs_all_anti), ichan)), n_shuf, g.alpha);
    P_sens_anti_fdr(logical(eye(size(P_sens_anti_fdr)))) = 1; % exclude diagonal elements from analysis, 1 will become zero after log transformation
    p_cmap = cmap_pvalues(P_sens_anti_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_anti_fdr, frqs, '', p_cmap, DIROUT, 'bispec_type', ['_cross_anti_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false)

    % run decomposition on normal cross-bispectrum
    n = [1 2 3]; % number of fitted source interactions
    err_colors = ['r', 'b', 'k'];
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L);

    % plotting (normal cross-bispectrum)
    plot_error(errors, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log', 'f_ext', '.fig')
    plot_error(errors, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear', 'f_ext', '.fig')
    m_order = 3; % only plot for n = 3 because it's easier to visualize
    p_cmap = cmap_pvalues(P_source_fdr{m_order}, cm17, cm17a);
    plot_bispec_slices(-log10(P_source_fdr{m_order}), [1 1 1], p_cmap, '' , DIROUT, 'cbar_label', '-log10(p)', 'f_ext', '.fig');

    % run decomposition on partially antisymmetrized cross-bispectrum
    antisymm = [2 1 3]; % will correspond to B_ijk - B_jik
    [P_source_anti_fdr, P_source_anti, F_anti, F_moca_anti, A_hat_anti, A_demixed_anti, D_hat_anti, D_demixed_anti, errors_anti] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L, 'antisymm', antisymm);

    % plotting (partially antisymmetrized cross-bispectrum)
    plot_error(errors_anti, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log_anti', 'f_ext', '.fig')
    plot_error(errors_anti, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear_anti', 'f_ext', '.fig')
    p_cmap = cmap_pvalues(P_source_anti_fdr{m_order}, cm17, cm17a);
    plot_bispec_slices(-log10(P_source_anti_fdr{m_order}), [1 2 2], p_cmap, '', DIROUT, 'cbar_label', '-log10(p)', 'f_name', '_anti', 'f_ext', '.fig');
    
    % run decomposition on totally antisymmetrized cross-bispectrum 
    [P_source_total_fdr, ~, ~, ~, ~, ~, ~, ~, errors_total] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L, 'total_antisymm', 'on');

    % plotting (totally antisymmetrized cross-bispectrum)
    plot_error(errors_total, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log_total', 'f_ext', '.fig')
    plot_error(errors_total, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear_total', 'f_ext', '.fig')
    p_cmap = cmap_pvalues(P_source_total_fdr{m_order}, cm17, cm17a);
    plot_pvalues_bispec_source(freqinds(1), freqinds(2), '', DIROUT, p_cmap, P_source_total_fdr{m_order}, P_source_total_fdr{m_order}, 'bispec_type', '_total', 'istitle', false, 'f_ext', '.fig')

    % localize interacting sources and plot topographies
    plot_topomaps_patterns(A_hat_anti{m_order}, m_order, EEG.chanlocs, cm17, isub, 'estimated', DIROUT, 'bispec_type', file_name) 

end
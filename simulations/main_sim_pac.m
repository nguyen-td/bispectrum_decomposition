% Pipeline to run the decomposition method on simulated PAC data for f1 = [5 7] Hz and f2 = [30 34] Hz with a
% specified number of univariate and bivariate interactions. 
%
% Input:
%   n_shuf         - [integer] number of shuffles
%   poolsize       - [integer] number of workers in the parellel pool (check parpool documentation) for parallel computing
%   n_univ         - [integer] number of univariate interactions, default is 0
%   n_biv          - [integer] number of bivariate interactions, default is 1
%   isnr           - [float] total signal-to-noise ratio which controls the ratio of signal (PAC) and 1/f noise on sensor level, 
%                     default is 0.5 (0 dB). 
%   isrand_simroi  - [boolean] whether to simulate interactions at randomly chosen ROIs. If true, regions are selected randomly. If false,
%                     interactions will be simulated in the occipital region (univariate PAC) and parietal-occipital region (bivariate PAC).
%                     See the 'iroi_pac' parameter in sim_wholebrain_pac.m for more details. Default is false
%
% Optional inputs:
%   alpha    - [float] significance level, default is 0.05.

function main_sim_pac(n_shuf, varargin)

    %% Setup
    local_path = '/data/tdnguyen/git_repos/';
    % local_path = '/Users/nguyentiendung/GitHub/';
    DIROUT = [local_path 'bispectrum_decomposition/simulations/sim_pac/figures/'];

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }              1;
        'n_univ'         'integer'       { }              0;
        'n_biv'          'integer'       { }              1;
        'isnr'           'float'         { }              0.5; 
        'isrand_simroi'  'boolean'       { }              false; 
        });
    if ischar(g), error(g); end

    % setup plotting
    load cm17

    %% Simulate PAC interactions 
    % check case
    if ~g.n_univ == 0 && g.n_biv == 0
        sim_case = 1;
    elseif g.n_univ == 0 && ~g.n_biv == 0
        sim_case = 2;
    else
        sim_case = 3;
    end

    rng('default')
    if g.isrand_simroi
        [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv, g.isnr);
    else
        if sim_case == 3
            roi_idx1 = 27; % lingual L
            roi_idx2 = 59; % superiorparietal L
            roi_idx3 = 28; % lingual R
            roi_idx4 = 60; % superiorparietal R
            roi_idx5 = 5; % 'caudalmiddlefrontal L'
            iroi_pac = [[roi_idx1 roi_idx2]; [roi_idx3 roi_idx4]; [roi_idx5 roi_idx5]];
            [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv, g.isnr, iroi_pac);
        else
            roi_idx1 = 27; % lingual L
            roi_idx2 = 59; % superiorparietal L, will be ignored for sim_case == 1
            [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv, g.isnr, [roi_idx1 roi_idx2]);
        end
    end
    
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
    plot_pvalues_univ(P_sens_fdr_uni, frqs, '', p_cmap, DIROUT, 'f_ext', '.fig', 'istitle', false)
    
    % analyze normal sensor cross-bispectrum
    freqpairs = get_freqindices(round_to_05(freqinds(1)), round_to_05(freqinds(2)), frqs);   
    para.nrun = n_shuf;
    [bs_all, bs_orig, ~] = data2bs_event_surro_final(signal_sensor(:, :)', segleng, segshift, epleng, freqpairs, para); % compute surrogate and true cross-bispectra
    [~, P_sens_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig), 3)), squeeze(mean(abs(bs_all), 3)), n_shuf, g.alpha);
    p_cmap = cmap_pvalues(P_sens_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_fdr, frqs, '', p_cmap, DIROUT, 'bispec_type', ['_cross_chan' int2str(3)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false, 'istitle', false)

    % analyze antisymmetrized sensor cross-bispectrum
    bs_orig_anti = bs_orig - permute(bs_orig, [2, 1, 3]); % B_ijk - B_jik
    bs_all_anti = bs_all - permute(bs_all, [2, 1, 3, 4]); % B_ijk - B_jik
    [~, P_sens_anti_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig_anti), 3)), squeeze(mean(abs(bs_all_anti), 3)), n_shuf, g.alpha);
    P_sens_anti_fdr(logical(eye(size(P_sens_anti_fdr)))) = 1; % exclude diagonal elements from analysis, 1 will become zero after log transformation
    p_cmap = cmap_pvalues(P_sens_anti_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_anti_fdr, frqs, '', p_cmap, DIROUT, 'bispec_type', ['_cross_anti_chan' int2str(3)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false, 'istitle', false)

    % run decomposition on normal cross-bispectrum
    n = [1 2 3]; % number of fitted source interactions
    err_colors = ['r', 'b', 'k'];
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L, 'bs_orig', bs_orig, 'bs_all', bs_all);

    % plotting (normal cross-bispectrum)
    plot_error(errors, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log', 'f_ext', '.fig')
    plot_error(errors, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear', 'f_ext', '.fig')
    m_order = 3; % only plot for n = 3 because it's easier to visualize
    p_cmap = cmap_pvalues(P_source_fdr{m_order}, cm17, cm17a);
    plot_bispec_slices(-log10(P_source_fdr{m_order}), [1 1 1], p_cmap, '' , DIROUT, 'cbar_label', '-log10(p)', 'f_ext', '.fig'); % only diagonals are relevant because it cannot be proven if off-diagonals result from interactions involving 2 or 3 sources
    plot_bispectra(D_demixed{m_order}, '', '', '', 'demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig')

    % run decomposition on partially antisymmetrized cross-bispectrum
    antisymm = [2 1 3]; % will correspond to B_ijk - B_jik
    [P_source_anti_fdr, P_source_anti, F_anti, F_moca_anti, A_hat_anti, A_demixed_anti, D_hat_anti, D_demixed_anti, errors_anti] = bsfit_stats(signal_sensor, freqinds(1), ...
        freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L, 'antisymm', antisymm, 'bs_orig', bs_orig, 'bs_all', bs_all);

    % plotting (partially antisymmetrized source cross-bispectrum)
    plot_error(errors_anti, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log_anti', 'f_ext', '.fig')
    plot_error(errors_anti, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear_anti', 'f_ext', '.fig')
    p_cmap = cmap_pvalues(P_source_anti_fdr{m_order}, cm17, cm17a);
    plot_bispec_slices(-log10(P_source_anti_fdr{m_order}), [1 2 2], p_cmap, '', DIROUT, 'cbar_label', '-log10(p)', 'f_name', '_anti', 'f_ext', '.fig');
    plot_bispectra(D_demixed_anti{m_order}, '', '', '', 'anti_demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig')
    
    if sim_case == 2 || sim_case == 3
        % run decomposition on totally antisymmetrized source cross-bispectrum 
        [P_source_total_fdr, ~, ~, ~, ~, ~, ~, D_demixed_total, errors_total] = bsfit_stats(signal_sensor, freqinds(1), ...
            freqinds(2), n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L, 'total_antisymm', 'on', 'bs_orig', bs_orig, 'bs_all', bs_all);

        % plotting (totally antisymmetrized source cross-bispectrum)
        plot_error(errors_total, 1, n, err_colors, '', '', DIROUT, 'islog', true, 'f_name', '_log_total', 'f_ext', '.fig')
        plot_error(errors_total, 1, n, err_colors, '', '', DIROUT, 'islog', false, 'f_name', '_linear_total', 'f_ext', '.fig')
        p_cmap = cmap_pvalues(P_source_total_fdr{m_order}, cm17, cm17a);
        plot_pvalues_bispec_source(freqinds(1), freqinds(2), '', DIROUT, p_cmap, P_source_total_fdr{m_order}, P_source_total_fdr{m_order}, 'bispec_type', '_total', 'istitle', false, 'f_ext', '.fig')
        plot_bispectra(D_demixed_total{m_order}, '', '', '', 'total_demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig')
    end

    % plot topographies and localized sources 
    chanlocs = readlocs('channel_BrainProducts_ActiCap_97.mat');
    figure; topoplot([], chanlocs, 'style', 'blank',  'electrodes', 'labelpoint'); 
    load bs_results
    if sim_case == 1 % plot results for normal (not antisymmetrized) bispectra
        plot_topomaps_patterns(A_hat{m_order}, m_order, chanlocs, cm17, '', 'estimated', DIROUT, 'f_ext', '.fig') 
        plot_topomaps_patterns(A_demixed{m_order}, m_order, chanlocs, cm17, '', 'demixed', DIROUT, 'f_ext', '.fig') 
        plot_sources(F_moca{m_order}, m_order, [], [], [], p_cmap, '', DIROUT, 'bispec_type', '', 'cortex_BS', cortex_highres, 'in_normal_to_high', in_normal_to_high)
    else % plot results for partially antisymmetrized bispectra
        plot_topomaps_patterns(A_hat_anti{m_order}, m_order, chanlocs, cm17, '', 'estimated', DIROUT, 'f_ext', '.fig') 
        plot_topomaps_patterns(A_demixed_anti{m_order}, m_order, chanlocs, cm17, '', 'demixed', DIROUT, 'f_ext', '.fig') 
        plot_sources(F_moca_anti{m_order}, m_order, [], [], [], p_cmap, '', DIROUT, 'bispec_type', '_anti', 'cortex_BS', cortex_highres, 'in_normal_to_high', in_normal_to_high)
    end
end
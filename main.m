% Wrapper method to run the pipeline that estimates the statistical significance of the
% bispectrum decomposition method for a single subject. Check bsfit_stats.m for the whole pipeline.
%
% Inputs:
%   n_shuf - [integer] number of shuffles
%   isub   - [integer] index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   n              - [integer] model order/number of fitted sources, default is 5. Can also be an array of n's, e.g., [3 4 5]
%   alpha          - [float] significance level, default is 0.05.
%   freq_manual    - [string] manual frequency selection, default is 'on', i.e., frequencies will not be selected automatically and have to be passed (see f1, f2).
%   f1             - [integer] phase frequency, single frequency in Hz or frequencyband, e.g., [9 11], default is 11.
%   f2             - [integer] amplitude frequency single frequency in Hz or frequency band, e.g., [20 22], default is 22.
%   run_ica        - [string] run ICA decomposition and save the first n components, default is 'off'.
%   poolsize       - [integer] number of workers in the parellel pool (check parpool documentation) for parallel computing, default is 2.
%   epleng         - [integer] length of epochs in seconds, default is 2 seconds
%   freq_down      - [integer] if 'downsample' is activated, the data will be downsampled to <freq_down> Hz. Default is 125 Hz.
%   downsample     - [string] check whether to downsample data to <freq_down> Hz, default is 'on'.
%   bispec_type    - [string] type of bispectrum (for file name), default is '_cross'.
%   antisymm       - [int, int, int] array of integer combinations to partially antisymmetrize, e.g. [1 3 2] would compute B_xyz - B_xzy. Default is [1 2 3] (no antisymmetrization).
%   train_test     - ['on' | 'off'] whether A should be fitted on train data and D on test data. If yes, the train-test split is 80-20. Default is 'off'.
%   dim_chan       - [integer] channel over which the bispectrum will becollapsed. Default is 2.

function main(n_shuf, isub, varargin)

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'                'integer'       { }                5;
        'alpha'            'float'         { }              0.05;
        'freq_manual'      'string'        { 'on' 'off' }   'on';
        'f1'               'integer'       { }              11; 
        'f2'               'integer'       { }              22; 
        'run_ica'          'string'        { 'on' 'off' }   'off';
        'poolsize'         'integer'       { }              1;
        'epleng'           'integer'       { }              2;
        'downsample'       'string'        { 'on' 'off'}    'on';
        'freq_down'        'integer'       { }              125;
        'bispec_type'      'string'        { }              ''; 
        'antisymm'         'integer'       { }              [1 2 3];
        'train_test'       'string'        { 'on' 'off'}    'off';
        'dim_chan'         'integer'       { 1 2 3}         2;
        });
    if ischar(g), error(g); end

    % set directory paths
    dict_name = ['fa' int2str(g.f1) '_fb' int2str(g.f2) '_' num2str(isub) '/'];
    DIROUT = ['/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/main_' dict_name];
    f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';
%     DIROUT = ['/data/tdnguyen/git_repos/bispectrum_decomposition/Lemon/figures/main_' dict_name];
%     f_path = '/data/tdnguyen/data/lemon/data/';

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % load data
    sub = ['sub-032' num2str(isub)];
    f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);

%     % epoching
%     epoch = [1 3]; % from 1 second after stimulus onset to 3 seconds after stimulus onset
%     EEG = eeg_checkset(EEG);
%     EEG = pop_epoch(EEG, { }, epoch, 'epochinfo', 'yes');

    % compute and plot ICA components (optional)
    if strcmpi(g.run_ica, 'on')
        run_ica(EEG, g.n, isub, DIROUT)
    end
    
    % downsample data to 100 Hz and plot
    if strcmpi(g.downsample, 'on')
        EEG = downsampling(EEG, g.freq_down);
        plot_spectra(EEG, 'EC', ['First peak: ' num2str(g.f1), 'Hz, Second peak: ' num2str(g.f2) ' Hz'], DIROUT, ...
            'title_str', ['psd_downsampled' num2str(isub)], 'f1', g.f1)
    end

    % set parameter values for (cross)-bispectrum estimation
    data = EEG.data;
    segleng = EEG.srate * g.epleng; 
    segshift = floor(segleng/2);
    epleng = EEG.srate * g.epleng; % create epochs of [e.epleng] seconds 
    fres = EEG.srate;
    frqs = sfreqs(fres, EEG.srate);

    f1 = g.f1;
    f2 = g.f2;
    
    % make frequency pre-selection by assessing the significance of frequency pairs of the univariate sensor bispectrum if 'freq_manual' = 'off'
    if strcmpi(g.freq_manual, 'off')
        [f1, f2, P_sens_fdr, P_sens] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    else
        f1 = g.f1;
        f2 = g.f2;
    end
    
    % analyze univariate sensor cross-bispectrum
    load cm17.mat 
    % [~, ~, P_sens_fdr_uni, ~] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    % p_cmap = cmap_pvalues(P_sens_fdr_uni, cm17, cm17a);
    % plot_pvalues_univ(P_sens_fdr_uni, frqs, isub, p_cmap, DIROUT, 'f_ext', '.fig', 'istitle', false)

    % analyze normal sensor cross-bispectrum
    dim_chan = g.dim_chan;
    freqpairs = get_freqindices(round_to_05(f1), round_to_05(f2), frqs); 
    para.nrun = n_shuf;
    [bs_all, bs_orig, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs, para); % compute surrogate and true cross-bispectra
    [~, P_sens_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig), dim_chan)), squeeze(mean(abs(bs_all), dim_chan)), n_shuf, g.alpha); % collapse over 2nd dimension
    p_cmap = cmap_pvalues(P_sens_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_fdr, frqs, isub, p_cmap, DIROUT, 'bispec_type', ['_cross_chan' int2str(dim_chan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false, 'istitle', false)
    
    % analyze antisymmetrized sensor cross-bispectrum
    bs_orig_anti = bs_orig - permute(bs_orig, g.antisymm); % B_ijk - B_jik
    bs_all_anti = bs_all - permute(bs_all, [g.antisymm, 4]); % B_ijk - B_jik
    [~, P_sens_anti_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig_anti), dim_chan)), squeeze(mean(abs(bs_all_anti), dim_chan)), n_shuf, g.alpha);
    P_sens_anti_fdr(logical(eye(size(P_sens_anti_fdr)))) = 1; % exclude diagonal elements from analysis, 1 will become zero after log transformation
    p_cmap = cmap_pvalues(P_sens_anti_fdr, cm17, cm17a);
    plot_pvalues_univ(P_sens_anti_fdr, frqs, isub, p_cmap, DIROUT, 'bispec_type', ['_cross_anti_chan' int2str(dim_chan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false, 'istitle', false)
    
    if ~(f1 == f2)
        % analyze totally antisymmetrized sensor cross-bispectrum
        bs_orig_total = bs_orig + permute(bs_orig, [3, 1, 2]) + permute(bs_orig, [2, 3, 1]) - permute(bs_orig, [3, 2, 1]) - permute(bs_orig, [2, 1, 3]) - permute(bs_orig, [1, 3, 2]); % TACB
        bs_all_total = bs_all + permute(bs_all, [3, 1, 2, 4]) + permute(bs_all, [2, 3, 1, 4]) - permute(bs_all, [3, 2, 1, 4]) - permute(bs_all, [2, 1, 3, 4]) - permute(bs_all, [1, 3, 2, 4]); % TACB
        [~, P_sens_total_fdr] = compute_pvalues(squeeze(mean(abs(bs_orig_total), dim_chan)), squeeze(mean(abs(bs_all_total), dim_chan)), n_shuf, g.alpha);
        P_sens_total_fdr(logical(eye(size(P_sens_total_fdr)))) = 1; % exclude diagonal elements from analysis, 1 will become zero after log transformation
        p_cmap = cmap_pvalues(P_sens_total_fdr, cm17, cm17a);
        plot_pvalues_univ(P_sens_total_fdr, frqs, isub, p_cmap, DIROUT, 'bispec_type', ['_cross_total_chan' int2str(dim_chan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'f_ext', '.fig', 'label_latex', false, 'istitle', false)
    end

    % test significance of the fitted normal, partially antisymmetrized and totally antisymmetrized source cross-bispectrum within subjects  
    [L_3D, cortex75k, cortex2k] = reduce_leadfield_nyhead(EEG); 
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(data, f1, f2, g.n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L_3D, 'train_test', g.train_test, 'bs_orig', bs_orig, 'bs_all', bs_all);    
    [P_source_anti_fdr, P_source_anti, F_anti, F_moca_anti, A_hat_anti, A_demixed_anti, D_hat_anti, D_demixed_anti, errors_anti] = bsfit_stats(data, f1, f2, g.n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L_3D, 'antisymm', g.antisymm, 'train_test', g.train_test, 'bs_orig', bs_orig, 'bs_all', bs_all);    
    if ~(f1 == f2)
        [P_source_total_fdr, P_source_total, F_total, F_moca_total, A_hat_total, A_demixed_total, D_hat_total, D_demixed_total, errors_total] = bsfit_stats(data, f1, f2, g.n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L_3D, 'total_antisymm', 'on', 'train_test', g.train_test, 'bs_orig', bs_orig, 'bs_all', bs_all);    
    end

    % save values
    save([DIROUT 'P_source_fdr.mat'], 'P_source_fdr', '-v7.3')
    save([DIROUT 'P_source_anti_fdr.mat'], 'P_source_anti_fdr', '-v7.3')
    if ~(f1 == f2)
        save([DIROUT 'P_source_total_fdr.mat'], 'P_source_total_fdr', '-v7.3')
        save([DIROUT 'D_demixed_total.mat'], 'D_demixed_total', '-v7.3')
        save([DIROUT 'F_moca_total.mat'], 'F_moca_total', '-v7.3')
        save([DIROUT 'A_demixed_total.mat'], 'A_demixed_anti', '-v7.3')
    end
    save([DIROUT 'bs_orig.mat'], 'bs_orig', '-v7.3')
    save([DIROUT 'bs_all.mat'], 'bs_all', '-v7.3')
    save([DIROUT 'D_demixed.mat'], 'D_demixed', '-v7.3')
    save([DIROUT 'D_demixed_anti.mat'], 'D_demixed_anti', '-v7.3')
    save([DIROUT 'F_moca.mat'], 'F_moca', '-v7.3')
    save([DIROUT 'F_moca_anti.mat'], 'F_moca_anti', '-v7.3')
    save([DIROUT 'A_demixed.mat'], 'A_demixed', '-v7.3')
    save([DIROUT 'A_demixed_anti.mat'], 'A_demixed_anti', '-v7.3')

    % plotting
    err_colors = ['r', 'b', 'k', 'g', 'o'];
    if max(g.n) > length(err_colors)
        err_colors = 0;
    end
        
    for n_idx = 1:length(g.n)
        m_order = g.n(n_idx);
        file_name = ['_n' num2str(m_order) g.bispec_type];

        if strcmpi(g.freq_manual, 'off')
            p_cmap = cmap_pvalues(P_sens_fdr, cm17, cm17a);
            plot_pvalues_univ(P_sens_fdr, frqs, isub, p_cmap, DIROUT, 'bispec_type', g.bispec_type)
        end

        % plot p-values and source cross-bispectra for normal, partially antisymmetrized and totally antisymmetrized source cross-bispectra
        p_cmap = cmap_pvalues(P_source_fdr{n_idx}, cm17, cm17a); % normal 
        plot_bispec_slices(-log10(P_source_fdr{n_idx}), [1 1 1], p_cmap, isub , DIROUT, 'cbar_label', '-log10(p)', 'f_ext', '.fig', 'f_name', file_name); % only diagonals are relevant because it cannot be proven if off-diagonals result from interactions involving 2 or 3 sources
        p_cmap = cmap_pvalues(D_demixed{n_idx}, cm17, cm17a); % normal 
        plot_bispectra(D_demixed{n_idx}, '', '', isub, 'demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig', 'bispec_type', file_name, 'dim_chan', g.dim_chan)

        p_cmap = cmap_pvalues(P_source_anti_fdr{n_idx}, cm17, cm17a); % partially  antisymmetrized
        plot_bispec_slices(-log10(P_source_anti_fdr{n_idx}), [1 2 2], p_cmap, isub , DIROUT, 'cbar_label', '-log10(p)', 'f_name', ['_anti' file_name], 'f_ext', '.fig'); 
        p_cmap = cmap_pvalues(D_demixed_anti{n_idx}, cm17, cm17a); % partially  antisymmetrized
        plot_bispectra(D_demixed_anti{n_idx}, '', '', isub, 'anti_demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig', 'bispec_type', file_name, 'dim_chan', g.dim_chan)

        if ~(f1 == f2)
            p_cmap = cmap_pvalues(P_source_total_fdr{n_idx}, cm17, cm17a); % totally  antisymmetrized
            plot_pvalues_bispec_source(f1, f2, isub, DIROUT, p_cmap, P_source_total_fdr{n_idx}, P_source_total{n_idx}, 'bispec_type', '_total', 'istitle', false, 'f_ext', '.fig', 'dim_chan', g.dim_chan)
            p_cmap = cmap_pvalues(D_demixed_anti{n_idx}, cm17, cm17a); % partially  antisymmetrized
            plot_bispectra(D_demixed_total{n_idx}, '', '', isub, 'total_demixed', DIROUT, p_cmap, 'istitle', false, 'f_ext', '.fig')
        end

        % plot estimated and demixed spatial patterns and demixed sources
        plot_topomaps_patterns(A_hat{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'estimated', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig') 
        plot_topomaps_patterns(A_demixed{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'demixed', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig')
        plot_sources(F_moca{n_idx}, m_order, cortex75k, cortex2k, [], cm17a, isub, DIROUT, 'bispec_type', file_name)

        plot_topomaps_patterns(A_hat_anti{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'anti_estimated', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig') 
        plot_topomaps_patterns(A_demixed_anti{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'anti_demixed', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig')
        plot_sources(F_moca_anti{n_idx}, m_order, cortex75k, cortex2k, [], cm17a, isub, DIROUT, 'bispec_type', ['_anti' file_name])

        if ~(f1 == f2)
            plot_topomaps_patterns(A_hat_total{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'estimated', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig') 
            plot_topomaps_patterns(A_demixed_total{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'demixed', DIROUT, 'bispec_type', file_name, 'f_ext', '.fig')
            plot_sources(F_moca_total{n_idx}, m_order, cortex75k, cortex2k, [], cm17a, isub, DIROUT, 'bispec_type', ['_total' file_name])
        end
    end

    % plot fitting errors
    plot_error(errors, 1, g.n, err_colors, '', isub, DIROUT, 'f_ext', '.fig')
    plot_error(errors, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_log', 'islog', true, 'f_ext', '.fig')

    plot_error(errors_anti, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_anti', 'f_ext', '.fig')
    plot_error(errors_anti, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_anti_log', 'islog', true, 'f_ext', '.fig')
    
    if ~(f1 == f2)
        plot_error(errors_total, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_total', 'f_ext', '.fig')
        plot_error(errors_total, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_total_log', 'islog', true, 'f_ext', '.fig')
    end

    % % store metrics
    % Errors = cellfun(@(x) x(end), errors);
    % Errors = Errors(:);
    % ModelOrder = g.n(:);
    % T = table(ModelOrder, Errors);
    % 
    % t_name = [DIROUT 'metrics_' num2str(isub) '.xlsx']; % name
    % writetable(T, t_name, 'Sheet', 1)

end
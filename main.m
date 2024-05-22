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
%   antisymm       - [int, int, int] array of integer combinations to antisymmetrize, e.g. [1 3 2] would compute B_xyz - B_xzy. Default is [1 2 3] (no antisymmetrization).
%   total_antisymm - ['on' | 'off'] perform total antisymmetrization Bartz et al. (2020), default is 'off'.
%   train_test     - ['on' | 'off'] whether A should be fitted on train data and D on test data. If yes, the train-test split is 80-20. Default is 'off'.

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
        'bispec_type'      'string'        { }              '_cross'; 
        'antisymm'         'integer'       { }              [1 2 3];
        'total_antisymm'   'string'        { 'on' 'off'}    'off';
        'train_test'       'string'        { 'on' 'off'}    'off';
        });
    if ischar(g), error(g); end

    if strcmpi(g.total_antisymm, 'on')
        anti_label = 'anti_total';
    elseif ~isequal(g.antisymm, [1, 2, 3])
        anti_label = ['anti_' strrep(num2str(g.antisymm), ' ', '')];
    else
        anti_label = '';
    end

    % set directory paths
    dict_name = [anti_label 'fa' int2str(g.f1) '_fb' int2str(g.f2) '_' num2str(isub) '/'];
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
    
    % make frequency pre-selection by assessing the significance of frequency pairs of the univariate sensor bispectrum if 'freq_manual' = 'off'
    if strcmpi(g.freq_manual, 'off')
        [f1, f2, P_sens_fdr, P_sens] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    else
        f1 = g.f1;
        f2 = g.f2;
    end
        
    % test significance of the fitted source cross-bispectrum within subjects
    [L_3D, cortex75k, cortex2k] = reduce_leadfield_nyhead(EEG); 
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, errors] = bsfit_stats(data, f1, f2, g.n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L_3D, 'antisymm', g.antisymm, 'total_antisymm', g.total_antisymm, 'train_test', g.train_test);    

    % plotting
    load cm17.mat 
    err_colors = ['r', 'b', 'k', 'g', 'o'];
    if max(g.n) > length(err_colors)
        err_colors = 0;
    end
        
    for n_idx = 1:g.n
        m_order = g.n(n_idx);
        file_name = ['_n' num2str(m_order) g.bispec_type];

        if strcmpi(g.freq_manual, 'off')
            p_cmap = cmap_pvalues(P_sens_fdr, cm17, cm17a);
            plot_pvalues_univ(P_sens_fdr, frqs, isub, p_cmap, DIROUT, 'bispec_type', g.bispec_type)
        end
        p_cmap = cmap_pvalues(P_source_fdr{n_idx}, cm17, cm17a);
        plot_pvalues_bispec_source(f1, f2, isub, DIROUT, p_cmap, P_source_fdr{n_idx}, P_source{n_idx}, 'bispec_type', file_name)
    
        % plot estimated and demixed spatial patterns
        plot_topomaps_patterns(A_hat{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'estimated', DIROUT, 'bispec_type', file_name) 
        plot_topomaps_patterns(A_demixed{n_idx}, m_order, EEG.chanlocs, cm17, isub, 'demixed', DIROUT, 'bispec_type', file_name)
    
        % plot sources
        plot_sources(F_moca{n_idx}, m_order, cortex75k, cortex2k, [], cm17a, isub, DIROUT, 'bispec_type', file_name)
        
        % plot D_hat and D_demixed
        plot_bispectra(D_hat{n_idx}, f1, f2, isub, 'estimated', DIROUT, cm17a, 'bispec_type', file_name)
        plot_bispectra(D_demixed{n_idx}, f1, f2, isub, 'demixed', DIROUT, cm17a, 'bispec_type', file_name)
    end
    % plot fitting error
    plot_error(errors, 1, g.n, err_colors, '', isub, DIROUT)
    plot_error(errors, 1, g.n, err_colors, '', isub, DIROUT, 'f_name', '_log', 'islog', true)

    % store metrics
    Errors = cellfun(@(x) x(end), errors);
    Errors = Errors(:);
    ModelOrder = g.n(:);
    T = table(ModelOrder, Errors);
    
    t_name = [DIROUT 'metrics_' num2str(isub) '.xlsx']; % name
    writetable(T, t_name, 'Sheet', 1)

end
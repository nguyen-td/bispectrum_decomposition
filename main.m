% Wrapper method to run the pipeline that estimates the statistical significance of the
% bispectrum decomposition method for a single subject. Check bsfit_stats.m for the whole pipeline.
%
% Inputs:
%   n_shuf - [integer] number of shuffles
%   isub   - [integer] index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   n            - [integer] model order/number of fitted sources, default is 5.
%   alpha        - [float] significance level, default is 0.05.
%   freq_manual  - [string] manual frequency selection, default is 'on', i.e., frequencies will not be selected automatically and have to be passed (see f1, f2).
%   f1           - [integer] phase frequency, single frequency in Hz or frequencyband, e.g., [9 11], default is 11.
%   f2           - [integer] amplitude frequency single frequency in Hz or frequency band, e.g., [20 22], default is 22.
%   run_ica      - [string] run ICA decomposition and save the first n components, default is 'off'.
%   poolsize     - [integer] number of workers in the parellel pool (check parpool documentation) for parallel computing, default is 2.
%   epleng       - [integer] length of epochs in seconds, default is 2 seconds
%   freq_down    - [integer] if 'downsample' is activated, the data will be downsampled to <freq_down> Hz. Default is 125 Hz.
%   downsample   - [string] check whether to downsample data to <freq_down> Hz, default is 'on'.
%   bispec_type  - [string] type of bispectrum (for file name), default is '_cross'.

function main(n_shuf, isub, varargin)
    
    % set directory paths
%     DIROUT = '/data/tdnguyen/data/p_imag'; % save directory
%     DIROUT = '/data/tdnguyen/data/p_carracer'; % save directory
%     DIROUT = ['/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/' num2str(isub) '/'];
%     f_path = '/data/tdnguyen/data/imag_data'; % change if necessary
%     f_path = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/EmergencyBreaking/preprocessing/analysis_output/preprocessing/data';
%     f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';
    DIROUT = ['/data/tdnguyen/git_repos/bispectrum_decomposition/Lemon/figures/main_' num2str(isub) '/'];
    f_path = '/data/tdnguyen/data/lemon/data/';

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'alpha'          'float'         { }              0.05;
        'freq_manual'    'string'        { 'on' 'off' }   'on';
        'f1'             'integer'       { }              11; 
        'f2'             'integer'       { }              22; 
        'run_ica'        'string'        { 'on' 'off' }   'off';
        'poolsize'       'integer'       { }              1;
        'epleng'         'integer'       { }              2;
        'downsample'     'string'        { 'on' 'off'}    'on';
        'freq_down'      'integer'       { }              125;
        'bispec_type'    'string'        { }              '_cross'; 
        });
    if ischar(g), error(g); end

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
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed] = bsfit_stats(data, f1, f2, g.n, n_shuf, frqs, segleng, segshift, epleng, g.alpha, L_3D);    

    % plot p-values
    load cm17.mat 
    if strcmpi(g.freq_manual, 'off')
        plot_pvalues_univ(P_sens_fdr, frqs, isub, cm17a, DIROUT, 'bispec_type', g.bispec_type)
    end
    plot_pvalues_bispec_source(f1, f2, isub, DIROUT, cm17a, P_source_fdr, P_source, 'bispec_type', g.bispec_type)

    % plot estimated and demixed spatial patterns
    plot_topomaps_patterns(A_hat, g.n, EEG.chanlocs, cm17, isub, 'estimated', DIROUT, 'bispec_type', g.bispec_type) 
    plot_topomaps_patterns(A_demixed, g.n, EEG.chanlocs, cm17, isub, 'demixed', DIROUT, 'bispec_type', g.bispec_type)

    % plot sources
    plot_sources(F_moca, F, g.n, cortex75k, cortex2k, [], cm17a, isub, DIROUT, 'bispec_type', g.bispec_type)
    
    % plot D_hat and D_demixed
    plot_bispectra(D_hat, f1, f2, isub, 'estimated', DIROUT, cm17a, 'bispec_type', g.bispec_type)
    plot_bispectra(D_demixed, f1, f2, isub, 'demixed', DIROUT, cm17a, 'bispec_type', g.bispec_type)

end
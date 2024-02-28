% Wrapper method to run the pipeline that estimates the statistical significance of the
% bispectrum decomposition method on motor imagery data for a single
% subject.
%
% Inputs:
%   nshuf - number of shuffles
%   isub  - index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   n           - model order/number of fitted sources, default is 5.
%   alpha       - significance level, default is 0.05.
%   freq_manual - manual frequency selection, default is 'off', i.e., frequencies will be selected automatically.
%   f1          - phase frequency, single frequency in Hz or frequencyband, e.g., [9 11]. Default is 11.
%   f2          - amplitude frequency single frequency in Hz or frequency band, e.g., [20 22]. Default is 22.
%   run_ica     - run ICA decomposition and save the first n components, default is 'off'
%   poolsize    - number of workers in the parellel pool (check parpool documentation) for parallel computing

function main(nshuf, isub, varargin)
    
    % set directory paths
%     DIROUT = '/data/tdnguyen/data/p_imag'; % save directory
%     DIROUT = '/data/tdnguyen/data/p_carracer'; % save directory
    DIROUT = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/';
%     f_path = '/data/tdnguyen/data/imag_data'; % change if necessary
%     f_path = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/EmergencyBreaking/preprocessing/analysis_output/preprocessing/data';
%     f_path = '/data/tdnguyen/git_repos/bispectrum_decomposition/EmergencyBreaking/preprocessing/analysis_output/preprocessing/data';
    f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';

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
        });
    if ischar(g), error(g); end

    % load data
    sub = 
    f_name = ['sub-032' isub '_EC.set']; % load LEMON eyes-closed data
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);

    % compute and plot ICA components (optional)
    if strcmpi(g.run_ica, 'on')
        run_ica(EEG, g.n, isub, DIROUT)
    end
    
    % set parameter values for (cross)-bispectrum estimation
    data = EEG.data;
    segleng = EEG.pnts;
    segshift = floor(segleng/2);
    epleng = EEG.pnts; 
    fres = EEG.srate;
    
    % make frequency pre-selection by assessing the significance of frequency pairs of the univariate sensor bispectrum if 'freq_manual' = 'off'
    frqs = sfreqs(fres, EEG.srate);
    if strcmpi(g.freq_manual, 'off')
        [f1, f2, P_sens_fdr, P_sens] = freq_preselection(data, nshuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    else
        f1 = g.f1;
        f2 = g.f2;
    end
        
    % test significance of the fitted source cross-bispectrum within subjects
    [L_3D, cortex75k, cortex2k] = reduce_leadfield(EEG); 
    [P_source_fdr, P_source, F, F_moca, A_hat, A_demixed, D_hat, D_demixed] = bsfit_stats(data, f1, f2, g.n, nshuf, frqs, segleng, segshift, epleng, g.alpha, L_3D);    

    % plot p-values
    if strcmpi(g.freq_manual, 'off')
        plot_pvalues(f1, f2, frqs, isub, DIROUT, P_source_fdr, P_source, P_sens_fdr, P_sens)
    else
        plot_pvalues(f1, f2, frqs, isub, DIROUT, P_source_fdr, P_source)
    end

    % plot estimated and demixed spatial patterns
    plot_topomaps(A_hat, g.n, EEG.chanlocs, isub, 'estimated', DIROUT) 
    plot_topomaps(A_demixed, g.n, EEG.chanlocs, isub, 'demixed', DIROUT)

    % plot sources
    load cm17
    plot_sources(F_moca, F, g.n, cortex75k, cortex2k, [], cm17, isub, DIROUT)
    
    % plot D_hat and D_demixed
    plot_bispectra(D_hat, f1, f2, isub, 'estimated', DIROUT)
    plot_bispectra(D_demixed, f1, f2, isub, 'demixed', DIROUT)

end
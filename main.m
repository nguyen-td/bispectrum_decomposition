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
%   f1          - phase frequency, default is 11 (mu rhythm)
%   f2          - amplitude frequency, default is 22 (beta rhythm)
%   run_ica     - run ICA decomposition and save the first n components, default is 'off'
%   poolsize    - number of workers in the parellel pool (check parpool documentation) for parallel computing

function main(nshuf, isub, varargin)
    
    % set paths
    DIROUT = '/data/tdnguyen/data/p_imag'; % save directory
    f_path = '/data/tdnguyen/data/imag_data'; % change if necessary
%     f_path = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/MotorImag/data';
%     f_chanlocs = 'MotorImag/data/chanlocs.mat';
    f_chanlocs = '/data/tdnguyen/data/imag_data/chanlocs.mat'; % change if necessary)

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'alpha'          'float'         { }              0.05;
        'freq_manual'    'string'        { 'on' 'off' }   'off';
        'f1'             'integer'       { }              11; 
        'f2'             'integer'       { }              22; 
        'run_ica'        'string'        { 'on' 'off' }   'off';
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end

    % load data
    sub = ['vp' num2str(isub)];
    f_name = ['prep_' sub '.set'];
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);

    % compute and plot ICA components (optional)
    if strcmpi(g.run_ica, 'on')
        run_ica(EEG, g.n, f_chanlocs, DIROUT)
    end
    
    % epoching
    epoch = [1 3]; % from 1 second after stimulus onset to 3 seconds after stimulus onset
    EEG = eeg_checkset(EEG);
    EEG = pop_epoch(EEG, { }, epoch, 'epochinfo', 'yes');
    
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
    L_3D = reduce_leadfield(EEG); 
    [P_source_fdr, P_source, A] = bsfit_stats(data, f1, f2, g.n, nshuf, frqs, segleng, segshift, epleng, g.alpha, L_3D);

    % create plots
    if strcmpi(g.freq_manual, 'off')
        plot_pvalues(A, f1, f2, frqs, isub, f_chanlocs, DIROUT, P_source_fdr, P_source, P_sens_fdr, P_sens)
    else
        plot_pvalues(A, f1, f2, frqs, isub, f_chanlocs, DIROUT, P_source_fdr, P_source)
    end

%     save_P = ['/P_' sub '.mat'];
%     save_A = ['/A_' sub '.mat'];
% 
%     save(strcat(DIROUT, save_P), '/P_source_fdr', '-v7.3');
%     save(strcat(DIROUT, save_A), '/A', '-v7.3');
    
end
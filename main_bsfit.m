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
%   freq_manual - manual frequency selection. Default is 'off', i.e., frequencies will be selected automatically.
%   f1          - phase frequency, default is 11 (mu rhythm)
%   f2          - amplitude frequency, default is 22 (beta rhythm)
%   poolsize    - number of workers in the parellel pool (check parpool documentation) for parallel computing

function main_bsfit(nshuf, isub, varargin)
    
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'alpha'          'float'         { }              0.05;
        'freq_manual'    'string'        { 'on' 'off' }   'off';
        'f1'             'integer'       { }              11; 
        'f2'             'integer'       { }              22; 
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end
    DIROUT = '/data/tdnguyen/data/p_imag'; % save directory

    % load data
    sub = ['vp' num2str(isub)];
    f_name = ['prep_' sub '.set'];
    f_path = '/data/tdnguyen/data/imag_data'; % change if necessary
%     f_path = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/MotorImag/data';
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);
    
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
    
    if strcmpi(g.freq_manual, 'off')
        % make frequency pre-selection by assessing the significane of frequency pairs of the univariate sensor bispectrum
        [f1, f2, P_sens_fdr, P_sens, frqs] = freq_preselection(data, nshuf, fres, EEG.srate, segleng, segshift, epleng, g.alpha, g.poolsize);
    else
        f1 = g.f1;
        f2 = g.f2;
    end
        
    % test significance of the fitted source cross-bispectrum within subjects
    [P_source_fdr, P_source, A] = bsfit_stats(data, f1, f2, g.n, nshuf, fres, EEG.srate, segleng, segshift, epleng, g.alpha);

    % create plots
    plot_pvalues(P_sens_fdr, P_sens, P_source_fdr, P_source, A, f1, f2, frqs, isub, DIROUT)

%     save_P = ['/P_' sub '.mat'];
%     save_A = ['/A_' sub '.mat'];
% 
%     save(strcat(DIROUT, save_P), '/P_source_fdr', '-v7.3');
%     save(strcat(DIROUT, save_A), '/A', '-v7.3');
    
end
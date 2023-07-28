% Wrapper method to run the pipeline that estimates the statistical significance of the
% bispectrum decomposition method on motor imagery data for a single
% subject.
%
% Inputs:
%   nshuf - number of shuffles
%   isub  - index of subject (the pipeline works for a single subject)

function [bs_all, bs_orig, P] = main_bsfit(nshuf, isub)
    DIROUT = '/data/tdnguyen/data/p_imag/';
    eeglab

    % load data
    sub = ['vp' num2str(isub)];
    f_name = ['prep_' sub '.set'];
    f_path = '/data/tdnguyen/data/imag_data'; % change if necessary
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);
    
    % epoching
    epoch = [1 3]; % from 1 second after stimulus onset to 3 seconds after stimulus onset
    EEG = eeg_checkset(EEG);
    EEG = pop_epoch(EEG, { }, epoch, 'epochinfo', 'yes');
    
    % set parameter values for cross-bispectrum estimation
    data = EEG.data;
    segleng = EEG.pnts;
    segshift = floor(segleng/2);
    epleng = EEG.pnts; 
    
    % test significance of the fitted source cross-bispectrum within subjects
    f1 = 11; %  mu
    f2 = 22; % beta
    n = 5;
    fres = EEG.srate;
    
    [bs_all, bs_orig, P] = run_bsfit(data, f1, f2, n, nshuf, fres, EEG.srate, segleng, segshift, epleng);

%     save_bsall = ['/bsall_' sub '.mat'];
%     save_bsorig = ['/bsorig_' sub '.mat'];
    save_P = ['/P_' sub '.mat'];

%     save(strcat(DIROUT, save_bsall), 'bs_all', '-v7.3');
%     save(strcat(DIROUT, save_bsorig), 'bs_orig', '-v7.3');
    save(strcat(DIROUT, save_P), 'P', '-v7.3');
    
end
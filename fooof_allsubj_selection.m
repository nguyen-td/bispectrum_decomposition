% Wrapper script to make a subject preselection for the subsequent analysis. 
% The metrics (SNR and difference between fundamental frequency and first harmonic) 
% are then stored in an Excel file which can then be used to sort by subject.

function fooof_allsubj_selection
    % setup
    f_path = '/Volumes/PortableSSD/LEMON/';
    sbjs = 301:1:528; % subject IDs
    
    snrs = [];
    diffs = [];
    eeglab
    for isub = sbjs
        disp(isub)
        try
            [snr, diff] = fooof_subj_selection(isub, f_path);
            disp(['Subject: ' int2str(isub)])

            % store values for all subjects
            snrs(end+1) = snr;
            diffs(end+1) = diff;
        catch
            continue 
        end
    end


end
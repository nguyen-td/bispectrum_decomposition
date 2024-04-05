% Wrapper script to make a subject preselection for the subsequent analysis. 
% The metrics (SNR and difference between fundamental frequency and first harmonic) 
% are then stored in an Excel file which can then be used to sort by subject.

function fooof_allsubj_selection
    % setup
    f_path = '/Volumes/PortableSSD/LEMON/';
    sbjs = 301:1:528; % subject IDs
    
    snrs = [];
    diffs = [];
    inds = [];
    eeglab
    for isub = sbjs
        try
            disp(['Subject: ' int2str(isub)])
            [snr, diff] = fooof_subj_selection(isub, f_path);

            % store values for all subjects
            snrs(end+1) = snr;
            diffs(end+1) = diff;
            inds(end+1) = isub;
        catch
            continue % skip if subject ID cannot be found
        end
    end

    results = vertcat(inds, snrs, diffs)';
    writematrix(results, 'subjects_snr.xlsx')
end
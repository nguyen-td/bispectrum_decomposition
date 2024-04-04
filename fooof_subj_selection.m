% Script to:
%   - Quantify the first harmonic by computing the SNR using FOOOF.
%   - Quantify the difference between the fundamental frequency and the first harmonic using FOOOF.
%
% This is to select subjects for the subsequent analysis.
%
% Input:
%   isub   - [integer] index of subject 
%   f_path - [String] path to data
%
% Outputs:
%   snr    - [double] signal-to-noise ratio (spectral value evaluated at the first harmonic)
%   diff   - [double] difference between alpha (fundamental frequency) and its first harmonic

function [snr, diff] = fooof_subj_selection(isub, f_path)

    % load data
    sub = ['sub-032' num2str(isub)];
    f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data

    % load preprocessed EEG
%     eeglab
    EEG = pop_loadset(f_name, f_path);

    % downsampling
    EEG = downsampling(EEG, 125);
    
    % run FOOOF
    [first_peak, ~, fooof_results] = find_peak_fooof(EEG);
%     fooof_plot(fooof_results)

    % correct spectrum
    ps_corrected = fooof_results.power_spectrum - fooof_results.ap_fit;
%     figure; plot(fooof_results.power_spectrum); hold on; plot(ps_corrected); plot(fooof_results.ap_fit); legend({'old', 'new', 'ap_fit'})

    % get metrics 
    snr = ps_corrected(fooof_results.freqs == floor(first_peak) * 2);
    diff = abs(ps_corrected(fooof_results.freqs == floor(first_peak)) - ps_corrected(fooof_results.freqs == floor(first_peak) * 2));
end

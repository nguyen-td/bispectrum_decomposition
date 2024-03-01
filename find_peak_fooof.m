% Find the peaks (fundametal frequency and its first harmonic) of the power spectrum using FOOOF.
%
% Input:
%   isub  - index of subject (the pipeline works for a single subject)
%
% Output:
%   first_peak  - mean (over channels) of the first peaks in the power spectrum (fundamental frequency)
%   second_peak - mean (over channels) of the first peaks in the power spectrum (first harmonic)

function [mean_first_peak, mean_second_peak] = find_peak_fooof(isub)

    % add paths and set directory
    addpath(genpath('/Users/nguyentiendung/GitHub/fooof_mat/'))
    addpath('/Users/nguyentiendung/GitHub/eeglab')
    f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';

    % load data
    eeglab
    sub = ['sub-032' num2str(isub)];
    f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data
    EEG = pop_loadset(f_name, f_path);

    % compute power spectrum
    [psd, freqs] = pwelch(EEG.data', 100, 50, [], EEG.srate);
    
    % FOOOF settings
    settings = struct();  
    settings.max_n_peaks = 2; % find fundamental frequency and its first harmonic
    f_range = [2, 40];

    % Run FOOOF
    n_chan = size(psd, 2);
    first_peak = zeros(1, n_chan);
    second_peak = zeros(1, n_chan);
    for ichan = 1:n_chan
        fooof_results = fooof(freqs', psd(:, ichan)', f_range, settings, true);
%         fooof_plot(fooof_results)

        % save peaks
        first_peak(ichan) = fooof_results.peak_params(1);
        second_peak(ichan) = fooof_results.peak_params(2);
    end

    % output mean values of peaks over channels
    mean_first_peak = mean(first_peak);
    mean_second_peak = mean(second_peak);
end

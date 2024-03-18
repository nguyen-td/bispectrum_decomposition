% Find the peaks (fundametal frequency and its first harmonic) of the power spectrum using FOOOF.
%
% Input:
%   EEG - EEG struct containing the sensor-level dataset
%
% Optional inputs:
%   alpha     - significance level, default is 0.05.
%   poolsize  - number of workers in the parellel pool (check parpool documentation) for parallel computing
%
% Output:
%   first_peak    - mean (over channels) of the first peaks in the power spectrum (fundamental frequency)
%   second_peak   - mean (over channels) of the first peaks in the power spectrum (first harmonic)
%   fooof_results - last FOOOF fit for visualization

function [first_peak, second_peak, fooof_results] = find_peak_fooof(EEG)

    % compute power spectrum
    [psd, freqs] = pwelch(EEG.data', 100, 50, 2*EEG.srate, EEG.srate);
    psd = mean(psd, 2);
    
    % FOOOF settings
    settings = struct();  
%     settings.max_n_peaks = 2; % find fundamental frequency and its first harmonic
    f_range = [1, 40];

    % run FOOOF
    py.importlib.import_module('fooof')
    fooof_results = fooof(freqs', psd, f_range, settings, true);
%         fooof_plot(fooof_results)

    % save peaks
    first_peak = fooof_results.peak_params(1);
    second_peak = fooof_results.peak_params(2);
end

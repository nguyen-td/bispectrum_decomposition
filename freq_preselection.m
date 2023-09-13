% Make a frequency pre-selection by evaluating the frequency pairs of the
% univariate sensor bispectrum. 
%
% Notations:
%   n_chans  - number of EEG channels (sensors)
%   n        - model order (number of estimated sources)
%   epleng   - length of a single epoch/trial
%   n_epochs - number of epochs/trials
%   n_freq   - number of frequencies
%
% Inputs:
%   data     - (n_chans x epleng x n_epochs) time series used to estimate the sensor cross-bispectrum
%   nshuf    - number of shuffles
%   fres     - frequency resolution
%   srate    - sampling rate
%   segleng  - segment length (see METH toolbox documentation)
%   segshift - overlap of segments (see METH toolbox documentation)
%   epleng   - epoch length (see METH toolbox documentation)
%   alpha    - significance level, default is 0.05.
%
% Outputs:
%   f1, f2 - frequencies in Hz
%   P_fdr  - (n_freq x n_freq) matrix of FDR-corrected p-values
%   P      - (n_freq x n_freq) matrix of p-values
%   frqs   - (n_frq x 1) array of frequencies

function [f1, f2, P_fdr, P, frqs] = freq_preselection(data, nshuf, fres, srate, segleng, segshift, epleng, alpha)
    
    % compute univariate sensor bispectrum
    clear para
    para.nrun = nshuf;
    frqs = sfreqs(fres, srate);
    
    disp('Start calculating surrogate univariate sensor bispectra for frequency selection...')
    [bsall, ~, ~] = data2bs_univar_stat(data(:, :)', segleng, segshift, epleng, length(frqs) - 1, para);
    
    % compute p-values, take mean over regions
    P = squeeze(sum(abs(mean(bsall(:, :, :, 1), 1)) < abs(mean(bsall(:, :, :, 2:end), 1)), 4) ./ nshuf);
    P(P==0) = 1 / nshuf;
    [~, argmax] = max(-log10(P(:)));
    [f1_bin, f2_bin] = ind2sub(size(P), argmax); 

    % convert to Hz
    f1 = frqs(f1_bin);
    f2 = frqs(f2_bin);

    % correct for multiple comparisons
    [p_fdr, ~] = fdr(P, alpha);
    P_fdr = P;
    P_fdr(P > p_fdr) = 1;

end
% Make a frequency pre-selection for the main pipeline by evaluating the frequency pairs of the
% univariate sensor bispectrum. Also serves as a wrapper function to compute univariate bispectra 
% and the corresponding p-values.
%
% Notations:
%   n_chans  - number of EEG channels (sensors)
%   n        - model order (number of estimated sources)
%   epleng   - length of a single epoch/trial
%   n_epochs - number of epochs/trials
%   n_freq   - number of frequencies
%   n_shufs  - number of shuffles
%
% Inputs:
%   data     - (n_chans x epleng x n_epochs) time series used to estimate the sensor cross-bispectrum
%   nshuf    - number of shuffles
%   frqs     - (n_frq x 1) array of frequencies
%   segleng  - segment length (see METH toolbox documentation)
%   segshift - overlap of segments (see METH toolbox documentation)
%   epleng   - epoch length (see METH toolbox documentation)
%   alpha    - significance level, default is 0.05.
%   poolsize - number of workers in the parellel pool (check parpool documentation) for parallel computing
%
% Outputs:
%   f1, f2       - frequencies in Hz
%   P_fdr        - (n_freq x n_freq) matrix of FDR-corrected p-values
%   P            - (n_freq x n_freq) matrix of p-values
%   bispec_orig  - (n_chans x n_freqs x n_freqs) surrogate univariate bispectral tensors (without normalization)
%   bicoh        - (n_chans x n_freqs x n_freqs) univariate bicoherence tensor

function [f1, f2, P_fdr, P, bispec_orig, bicoh] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, alpha, poolsize)
    
    % compute univariate sensor bispectrum
    clear para
    para.nrun = n_shuf;
    maxfreqbins = floor(segleng/2);
    
    disp('Start calculating surrogate univariate sensor bispectra for frequency selection...')
    parpool(poolsize)
    tic
    [bsall, bsallnr] = data2bs_univar_stat(data(:, :)', segleng, segshift, epleng, maxfreqbins, para);
    toc
    
    % shut down current parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    % compute bicoherence
    bispec_orig = squeeze(bsall(:, :, :, 1)); % original bispectral tensor
    bicoh = bispec_orig ./ bsallnr;

    % compute p-values
    [P, P_fdr] = compute_pvalues(squeeze(mean(abs(bsall(:, :, :, 1)), 1)), squeeze(mean(abs(bsall(:, :, :, 2:end)), 1)), n_shuf, alpha);
    
%     % compute p-values, take mean over regions
%     P = squeeze(sum(abs(mean(bsall(:, :, :, 1), 1)) < abs(mean(bsall(:, :, :, 2:end), 1)), 4) ./ nshuf);
%     P(P==0) = 1 / nshuf;
% 
%     % correct for multiple comparisons
%     [p_fdr, ~] = fdr(P, alpha);
%     P_fdr = P;
%     P_fdr(P > p_fdr) = 1;
%     
    % extract frequencies
    [maxval, ~] = max(-log10(P_fdr(:)));
    argmaxs = find(-log10(P_fdr(:)) == maxval);
    argmax = argmaxs(randi(length(argmaxs), 1));
    [f1_bin, f2_bin] = ind2sub(size(P_fdr), argmax); 

    % convert to Hz
    f1 = frqs(f1_bin);
    f2 = frqs(f2_bin);

end
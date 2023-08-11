% Test significance of the fitted source cross-bispectrum within 
% subjects based on the following steps:
%
% 1. Estimate the (true) sensor cross-bispectrum B_true.
% 2. Run the decomposition method on B_true to estimate A_hat (true mixing
%    matrix) and D_hat (true source cross-bispectrum)
% 3. Get a null distribution by computing surrogate sensor cross-bispectra
%    B_shuf 
% 4. Estimate surrogate source cross-bispectra D_shuf by holding the mixing 
%    matrtix A_hat fixed (use/fix A_hat and only fit D_shuf)
% 5. Compute p-values based on the absolute values of the estimate source cross-bispectra.
%    The result will be a (n x n x n) tensor of p-values.
% 6. TO-DO: Perform source localization using significant A_hat.
%
% Notations:
%   k - number of EEG channels (sensors)
%   n - model order (number of estimated sources)
%   T - length of a time series (number of data points)
%
% Inputs:
%   data        - (k, T) time series used to estimate the sensor cross-bispectrum
%   f1          - single frequency in Hz (low frequency)
%   f2          - single frequency in Hz (high frequency)
%   n           - number of estimated sources
%   nshuf       - number of shuffles
%   fres        - frequency resolution
%   srate       - sampling rate
%   segleng     - segment length (see METH toolbox documentation)
%   segshift    - overlap of segments (see METH toolbox documentation)
%   epleng      - epoch length (see METH toolbox documentation)

function [bs_all, bs_orig, P, A_hat] = bsfit_stats(data, f1, f2, n, nshuf, fres, srate, segleng, segshift, epleng, alpha)

    % estimate sensor cross-bispectrum
    clear para
    para.nrun = nshuf;
    frqs = sfreqs(fres, srate);
    freqpairs = [find(frqs == f1), find(frqs == f2)];

    disp('Start calculating surrogate sensor cross-bispectra...')
%     [bs_all, bs_orig, ~] = data2bs_event_surro_final(data', segleng, segshift, epleng, freqpairs, para);
    [bs_all, bs_orig, ~] = data2bs_event_surro_final(data(:,:)', segleng, segshift, epleng, freqpairs, para);

    % run decomposition on the original sensor cross-bispectrum 
    [A_hat, D_hat, ~, ~, ~] = bsfit(bs_orig, n);

    % fit surrogate source cross-bispectra with fixed mixing matrix 
    disp('Start calculating surrogate source cross-bispectra...')
    clear para
    para.a = A_hat;
    para.isfit_a = false;

    fprintf('Progress of %d:', nshuf);
    D_shuf = zeros(n, n, n, nshuf);
    for ishuf = 1:nshuf
        if mod(ishuf, 10) == 0
            fprintf('%d', ishuf);
        elseif mod(ishuf, 2) == 0
            fprintf('.');
        end

        [~, D_ishuf, ~, ~, ~] = bsfit(bs_all(:, :, :, ishuf), n, para);
        D_shuf(:, :, :, ishuf) = D_ishuf;
    end
    
    % compute p-values
    P = sum(abs(D_hat) < abs(D_shuf), 4) ./ nshuf;

    % correct for multiple comparisons
    [p_fdr, ~] = fdr(P, alpha);
    P(P > p_fdr) = 1;
    P(P==0) = 1 / nshuf;

end
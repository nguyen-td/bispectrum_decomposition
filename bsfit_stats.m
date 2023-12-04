% Test significance of the fitted source cross-bispectrum within 
% subjects based on the following steps:
%
% 1. Estimate the (true) sensor cross-bispectrum B_true.
% 2. Run the decomposition method on B_true to estimate A_hat (true mixing
%    matrix) and D_hat (true source cross-bispectrum)
% 3. Get a null distribution by computing surrogate sensor cross-bispectra
%    B_shuf 
% 4. Estimate surrogate source cross-bispectra D_shuf by holding the mixing 
%    matrix A_hat fixed (use/fix A_hat and only fit D_shuf)
% 5. Unmix the estimated source interactions using MOCA (applied on A_hat).
% 6. Compute p-values based on the absolute values of the estimate source cross-bispectra.
%    The result will be a (n x n x n) tensor of p-values.
% 7. TO-DO: Perform source localization using significant A_hat.
%
% Notations:
%   n_chans  - number of EEG channels (sensors)
%   n        - model order (number of estimated sources)
%   epleng   - length of a single epoch/trial
%   n_epochs - number of epochs/trials
%
% Inputs:
%   data        - (n_chans x epleng x n_epochs) time series used to estimate the sensor cross-bispectrum
%   f1          - single frequency in Hz (low frequency)
%   f2          - single frequency in Hz (high frequency)
%   n           - number of estimated sources
%   nshuf       - number of shuffles
%   frqs        - (n_frq x 1) array of frequencies
%   srate       - sampling rate
%   segleng     - segment length (see METH toolbox documentation)
%   segshift    - overlap of segments (see METH toolbox documentation)
%   epleng      - epoch length (see METH toolbox documentation)
%   alpha       - significance level, default is 0.05.
%   L_3D        - (n_chans x n_voxels x n_dum) leadfield tensor, dipole directions are typically 3 
%
% Outputs:
%   P_fdr  - (n x n x n) tensor of fdr-corrected p-values
%   P      - (n x n x n) tensor of p-values (before fdr correction)
%   A_sens - (n_chan x n) mixing matrix

function [P_fdr, P, A_sens] = bsfit_stats(data, f1, f2, n, nshuf, frqs, segleng, segshift, epleng, alpha, L_3D)

    % estimate sensor cross-bispectrum
    clear para
    para.nrun = nshuf;
    freqpairs = [find(frqs == f1), find(frqs == f2)];

    disp('Start calculating surrogate sensor cross-bispectra...')
    [bs_all, bs_orig, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs, para);

    % run decomposition on the original sensor cross-bispectrum 
    [A_hat, D_hat, ~, ~, ~] = bsfit(bs_orig, n);

    % unmix source interactions using MOCA
    A_moca = apply_moca(L_3D, A_hat, n);
    A_sens = A_hat * A_moca; % combined sensor patterns

    % fit surrogate source cross-bispectra with fixed mixing matrix 
    disp('Start calculating surrogate source cross-bispectra...')
    clear para d
    para.a = A_sens;
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
    P(P==0) = 1 / nshuf;

    % correct for multiple comparisons
    [p_fdr, ~] = fdr(P, alpha);
    P_fdr = P;
    P_fdr(P > p_fdr) = 1;

end
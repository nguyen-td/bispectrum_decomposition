% Test significance of the fitted source cross-bispectrum within 
% subjects based on the following steps:
%
% 1. Estimate the (true) sensor cross-bispectrum B_true.
% 2. Run the decomposition method on B_true to estimate A_hat (true mixing
%    matrix) and D_hat (true source cross-bispectrum)
% 3. Fit surrogate source cross-bispectra by holding the mixing matrix A_hat fix (i.e., only fit D_shuf).
% 4. Unmix the estimated pattern A_hat and cross-bispectrum D_hat using MOCA to get the demixed pattern 
%    A_demixed and demixed cross-bispectrum D_demixed.
% 5. Perform source localization using significant demixed pattern A_demixed.
% 6. Compute p-values based on the absolute values of the estimate source cross-bispectra.
%    The result will be a (n x n x n) tensor of p-values.
%
% Notations:
%   n_chans  - number of EEG channels (sensors)
%   n        - model order (number of estimated sources)
%   epleng   - length of a single epoch/trial
%   n_epochs - number of epochs/trials
%
% Inputs:
%   data           - (n_chans x epleng x n_epochs) time series used to estimate the sensor cross-bispectrum
%   f1             - single frequency in Hz or frequency band, e.g., [9 11](low frequency)
%   f2             - single frequency in Hz or frequency band, e.g., [22 24] (high frequency)
%   n              - number of estimated sources
%   nshuf          - number of shuffles
%   frqs           - (n_frq x 1) array of frequencies
%   srate          - sampling rate
%   segleng        - segment length (see METH toolbox documentation)
%   segshift       - overlap of segments (see METH toolbox documentation)
%   epleng         - epoch length (see METH toolbox documentation)
%   alpha          - significance level, default is 0.05.
%   L_3D           - (n_chans x n_voxels x n_dum) leadfield tensor, dipole directions are typically 3 
%   antisymm       - [idx, idx, idx] array containing indices to permute, will not perform antisymmetrization if [1, 2, 3]
%   total_antisymm - {'on' | 'off'} whether TACB should be compuuted
%
% Outputs:
%   P_fdr     - (1 x n) cell array with (n x n x n) tensors of fdr-corrected p-values
%   P         - (1 x n) cell array with (n x n x n) tensors of p-values (before fdr correction)
%   F         - (1 x n) cell array with (n_voxels x n_dum x n) mixed sources
%   F_moca    - (1 x n) cell array with (n_voxels x n_dum x n) demixed sources after MOCA
%   A_hat     - (1 x n) cell array with (n_chan x n) estimated spatial patterns (mixing matrix)
%   A_demixed - (1 x n) cell array with (n_chan x n) demixed spatial patterns 
%   D_hat     - (1 x n) cell array with (n x n x n) estimated source cross-bispectra
%   D_demixed - (1 x n) cell array with (n x n x n) demixed source cross-bispectra
%   err       - (1 x n) cell array with (n_freqcombs x n) errors over iterations

function [P_fdr, P, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, err] = bsfit_stats(data, f1, f2, n, nshuf, frqs, segleng, segshift, epleng, alpha, L_3D, antisymm, total_antisymm)
    
    % get frequency pairs (in bins)
    freqpairs = get_freqindices(round_to_05(f1), round_to_05(f2), frqs); 

    % estimate sensor cross-bispectrum
    clear para
    para.nrun = nshuf;
    disp('Start calculating surrogate sensor cross-bispectra...')
    [bs_all, bs_orig] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs, para);

    % antisymmetrization
    if strcmpi(total_antisymm, 'on')
        disp('Perform total antisymmetrization')
        bs_orig = bs_orig + permute(bs_orig, [3, 1, 2]) + permute(bs_orig, [2, 3, 1]) - permute(bs_orig, [3, 2, 1]) - permute(bs_orig, [2, 1, 3]) - permute(bs_orig, [1, 3, 2]);
    elseif ~isequal(antisymm, [1, 2, 3])
        disp('Perform partial antisymmetrization')
        bs_orig = bs_orig - permute(bs_orig, antisymm); 
    else
        disp('No antisymmetrization')
    end

    % run decomposition on the original sensor cross-bispectrum  
    A_hat = {};
    D_hat = {};
    err = {};
    for m_order = n
%        [A_hat, D_hat] = bsfit(bs_orig, n);
       [A_hat{end+1}, D_hat{end+1}, ~, ~, err{end+1}] = bsfit_freqbands(bs_orig, m_order); 
    end

    % fit surrogate source cross-bispectra with fixed mixing matrix 
    disp('Start calculating surrogate source cross-bispectra...')
    D_shuf = {};
    for n_idx = 1:length(n)
        clear para
        para.a = A_hat{n_idx};
        para.isfit_a = false;
        m_order = n(n_idx);
        disp(['Fitting for model order ' num2str(m_order) ':...'])

        fprintf('Progress of %d:', nshuf);
        D_shuf_n = zeros(m_order, m_order, m_order, nshuf);
        for ishuf = 1:nshuf
            if mod(ishuf, 10) == 0
                fprintf('%d', ishuf);
            elseif mod(ishuf, 2) == 0
                fprintf('.');
            end
            
    %         [~, D_ishuf] = bsfit(bs_all(:, :, :, ishuf), n, para);
            [~, D_ishuf] = bsfit_freqbands(bs_all(:, :, :, ishuf), m_order, para);
            D_shuf_n(:, :, :, ishuf) = D_ishuf;
        end
        D_shuf{end+1} = D_shuf_n;
        fprintf('\n');
    end
    fprintf('\n');

    % demix the source interactions using MOCA and compute p-values
    P = {};
    P_fdr = {};
    F = {};
    F_moca = {};
    A_demixed = {};
    D_demixed = {};
    for n_idx = 1:length(n)
        [A_moca, F_moca{end+1}, F{end+1}] = apply_moca(L_3D, A_hat{n_idx}, n(n_idx));
        A_demixed{end+1} = A_hat{n_idx} * A_moca'; % demix sensor pattern
        D_demixed{end+1} = calc_bsmodel(A_moca', D_hat{n_idx}); % demix source cross-bispectrum

        % compute p-values
        [P{end+1}, P_fdr{end+1}] = compute_pvalues(abs(D_hat{n_idx}), abs(D_shuf{n_idx}), nshuf, alpha, 'shuf_dim', 4);
    end

end
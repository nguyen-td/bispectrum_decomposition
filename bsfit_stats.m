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
%   segleng        - segment length (see METH toolbox documentation)
%   segshift       - overlap of segments (see METH toolbox documentation)
%   epleng         - epoch length (see METH toolbox documentation)
%   alpha          - significance level, default is 0.05.
%   L_3D           - (n_chans x n_voxels x n_dum) leadfield tensor, dipole directions are typically 3 
%
% Optional inputs:
%   antisymm       - [idx, idx, idx] array containing indices to permute, will not perform antisymmetrization if [1, 2, 3]. Default is [1, 2, 3].
%   total_antisymm - ['on' | 'off'] whether TACB should be compuuted. Default is 'off'.
%   train_test     - ['on' | 'off'] whether A should be fitted on train data and D on test data. If yes, the train-test split is 80-20. Default is 'off'.
%   bs_orig        - (n_chans x n_chans x n_chans) original sensor cross-bispectrum. Default is empty.
%   bs_all         - (n_chans x n_chans x n_chans x nshuf) surrogate sensor cross-bispectrum. Default is empty.

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
%   bs_orig   - (n_chans x n_chans x n_chans) original sensor cross-bispectrum
%   bs_all    - (n_chans x n_chans x n_chans x nshuf) surrogate sensor cross-bispectrum

function [P_fdr, P, F, F_moca, A_hat, A_demixed, D_hat, D_demixed, err, bs_orig, bs_all] = bsfit_stats(data, f1, f2, n, nshuf, frqs, segleng, segshift, epleng, alpha, L_3D, varargin)
    
    g = finputcheck(varargin, { ...
        'antisymm'         'integer'       { }               [1 2 3];
        'total_antisymm'   'string'        { 'on' 'off' }    'off';
        'train_test'       'string'        { 'on' 'off' }    'off';
        'bs_orig'          'float'         { }               [];
        'bs_all'           'float'         { }               []; 
        });
    if ischar(g), error(g); end
    
    % get frequency pairs (in bins)
    freqpairs = get_freqindices(round_to_05(f1), round_to_05(f2), frqs); 
    
    if isempty(g.bs_orig) && isempty(g.bs_all)
        % estimate sensor cross-bispectrum, either on full dataset or on train-test splots
        clear para
        disp('Start calculating surrogate sensor cross-bispectra...')
        if strcmpi(g.train_test, 'off')
            para.nrun = nshuf;
            [bs_all, bs_orig] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs, para);
        else
            % fit A_hat and D_hat on bs_orig that is computed using training data, later fit D_shuf on bs_all that is computed using test data
            cut = round(size(data, 2) * 0.8);
            train = data(:, 1:cut);
            test = data(:, cut+1:end);
    
            para.nrun = 1; % original bispectrum on training data
            [~, bs_orig] = data2bs_event_surro_final(train', segleng, segshift, epleng, freqpairs, para);
    
            para.nrun = nshuf; % surrogates on test data
            [bs_all, ~] = data2bs_event_surro_final(test', segleng, segshift, epleng, freqpairs, para);
        end
    else
        bs_orig = g.bs_orig;
        bs_all = g.bs_all;
    end

    % antisymmetrization
    if strcmpi(g.total_antisymm, 'on')
        disp('Perform total antisymmetrization')
        bs_orig = bs_orig + permute(bs_orig, [3, 1, 2]) + permute(bs_orig, [2, 3, 1]) - permute(bs_orig, [3, 2, 1]) - permute(bs_orig, [2, 1, 3]) - permute(bs_orig, [1, 3, 2]);
        bs_all = bs_all + permute(bs_all, [3, 1, 2, 4]) + permute(bs_all, [2, 3, 1, 4]) - permute(bs_all, [3, 2, 1, 4]) - permute(bs_all, [2, 1, 3, 4]) - permute(bs_all, [1, 3, 2, 4]);
    elseif ~isequal(g.antisymm, [1, 2, 3])
        disp('Perform partial antisymmetrization')
        bs_orig = bs_orig - permute(bs_orig, g.antisymm); 
        bs_all = bs_all - permute(bs_all, [g.antisymm, 4]); 
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
        
        % % demix shuffles
        % D_shuf_demixed = zeros(n_idx, n_idx, n_idx, nshuf);
        % for ishuf = 1:n_idx
        %     D_shuf_demixed(:, :, :, ishuf) = calc_bsmodel(A_moca', D_shuf{n_idx}(:, :, :, ishuf));
        % end

        % compute p-values
        % [P{end+1}, P_fdr{end+1}] = compute_pvalues(abs(D_demixed{n_idx}), abs(D_shuf_demixed), nshuf, alpha, 'shuf_dim', 4);
        [P{end+1}, P_fdr{end+1}] = compute_pvalues(abs(D_demixed{n_idx}), abs(D_shuf{n_idx}), nshuf, alpha, 'shuf_dim', 4);
    end

end
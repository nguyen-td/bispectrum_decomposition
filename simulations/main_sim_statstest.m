% Pipeline to test the FPR of the statistical test. To this end, a single univariate PAC interaction 
% is simulated and placed randomly in one of the 68 Desikan-Killiany atlas regions. Non-interacting 
% 1/f noise is placed in all other regions. The time series are then projected to sensor level. The 
% decompositin method is then applied to this data.
%
% Input:
%   n_shuf - [integer] number of shuffles (for computing surrogates)
%   n_iter - [integer] number of iterations (number of PAC scores to be computed)
%
% Optional inputs:
%   n          - [integer] model order/number of fitted sources, default is 1. Can also be an array of n's, e.g., [3 4 5]
%   n_univ     - [integer] number of univariate interactions, default is 1
%   n_biv      - [integer] number of bivariate interactions, default is 0
%   isnr       - [float] total signal-to-noise ratio which controls the ratio of signal (PAC) and 1/f noise on sensor level, 
%                default is 0.5 (0 dB). Can also be an array of SNRs, e.g., [0.2 0.4 0.5 0.6 0.8].
%   alpha      - [float] significance level, default is 0.05.
%   epleng     - [integer] epoch length (see METH toolbox documentation)
%   train_test - ['on' | 'off'] whether A should be fitted on train data and D on test data. If yes, the train-test split is 80-20. Default is 'off'.

function main_sim_statstest(n_shuf, n_iter, varargin)
    
    local_path = '/data/tdnguyen/git_repos/';
%     local_path = '/Users/nguyentiendung/GitHub/';
    DIROUT = [local_path 'bispectrum_decomposition/simulations/results/'];

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                1;
        'n_univ'         'integer'       { }                1;
        'n_biv'          'integer'       { }                0;
        'isnr'           'float'         { }               0.5; 
        'epleng'         'integer'       { }                2;
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }                1;
        'train_test'     'string'        { 'on' 'off' }    'off';
        });
    if ischar(g), error(g); end

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % check case
    if ~g.n_univ == 0 && g.n_biv == 0
        sim_case = 1;
    elseif g.n_univ == 0 && ~g.n_biv == 0
        sim_case = 2;
    else
        sim_case = 3;
    end
%     f_name = ['_snr' int2str(20 * log10(g.isnr / (1 - g.isnr))) '_case' int2str(sim_case)];
    
    for snr = g.isnr
        disp(['Test statistical test for SNR ' num2str(snr) '...'])
        % generate simulated data
        [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv, snr);
    
        % sampling frequency
        fres = fs; 
        frqs = sfreqs(fres, fs); % freqs in Hz
    
        % set parameter values for (cross)-bispectrum estimation
        freqinds = [mean(filt.low) mean(filt.high)]; % in Hz
        segleng = fs * g.epleng; 
        segshift = floor(segleng/2);
        epleng = fs * g.epleng; % create epochs of [e.epleng] seconds 
    
        % run statistics on decomposition and compute FPR
        disp(['Start computing FPR ' int2str(n_iter) ' times...'])
        parpool(g.poolsize)
        fpr_iter = zeros(length(g.n), n_iter); 
        P_fdr = {};
        tic
        parfor i_iter = 1:n_iter
            if strcmpi(g.train_test, 'off')
                P_fdr{i_iter} = bsfit_stats(signal_sensor, freqinds(1), freqinds(2), g.n, n_shuf, frqs, ...
                    segleng, segshift, epleng, g.alpha, L);
            else
                P_fdr{i_iter} = bsfit_stats(signal_sensor, freqinds(1), freqinds(2), g.n, n_shuf, frqs, ...
                    segleng, segshift, epleng, g.alpha, L, 'train_test', 'on');
            end
            fpr = cellfun(@(x) x < g.alpha, P_fdr{i_iter}, 'UniformOutput', false);
            fpr_iter(:, i_iter) = cell2mat(cellfun(@(x) sum(x, 'all'), fpr, 'UniformOutput', false));
        end
        toc
    
        % shut down current parallel pool
        poolobj = gcp('nocreate');
        delete(poolobj);
    
        % save structs
        save([DIROUT 'P_fdr_traintest_' g.train_test '_snr' num2str(snr) '_case' int2str(sim_case) '.mat'], 'P_fdr', '-v7.3')
        save([DIROUT 'FPR_traintest_' g.train_test '_snr' num2str(snr) '_case' int2str(sim_case) '.mat'], 'fpr_iter', '-v7.3')
    end

%     % raincloud/half-violin plot on linearly scaled y-axis
%     for i_source = 1:g.n
%         titles = {'Normal', 'Train-Test Split'};
%         colors = [[0 0 0.5]; [0.8 0 0.2]];
%         plot_metrics_raincloud(fpr_iter(i_source, :), colors, 1, 0.2, 'ks', titles, DIROUT, 'name', ['_n' int2str(i_source) f_name]);
%     end
    
end
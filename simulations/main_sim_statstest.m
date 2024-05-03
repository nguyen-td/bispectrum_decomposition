% Pipeline to test FPR of the statistical test. To this end, univariate PAC
% is simulated and decomposed.
%
% Input:
%   n_shuf - [integer] number of shuffles (for computing surrogates)
%   n_iter - [integer] number of iterations (number of PAC scores to be computed)
%
% Optional inputs:
%   n        - [integer] model order/number of fitted sources, default is 1. Can also be an array of n's, e.g., [3 4 5]
%   n_univ   - [integer] number of univariate interactions, default is 1
%   n_biv    - [integer] number of bivariate interactions, default is 0
%   alpha    - [float] significance level, default is 0.05.
%   epleng   - [integer] epoch length (see METH toolbox documentation)

function main_sim_statstest(n_shuf, n_iter, varargin)
    
    local_path = '/data/tdnguyen/git_repos/';
%     local_path = '/Users/nguyentiendung/GitHub/';
    DIROUT = [local_path 'bispectrum_decomposition/simulations/figures/stats_test/'];

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                1;
        'n_univ'         'integer'       { }                1;
        'n_biv'          'integer'       { }                0;
        'epleng'         'integer'       { }                2;
        'alpha'          'float'         { }              0.05;
        });
    if ischar(g), error(g); end

    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % check case
    if ~g.n_univ == 0 && g.n_biv == 0
        sim_case = 1;
        isnr = 0.8;
    elseif g.n_univ == 0 && ~g.n_biv == 0
        sim_case = 2;
        isnr = 0.5;
    else
        sim_case = 3;
        isnr = 0.5;
    end
    f_name = ['_snr' int2str(20 * log10(isnr / (1 - isnr))) '_case' int2str(sim_case)];

    % generate simulated data
    [signal_sensor, fs, source, filt, L] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv, isnr);
%     psd = pwelch(signal_sensor(:,:)', 100, 50, 2*fs, fs);
    psd = pwelch(source', 100, 50, 2*fs, fs);

    % sampling frequency
    fres = fs; 
    frqs = sfreqs(fres, fs); % freqs in Hz
    plot_psd_pac(psd, frqs, DIROUT, 'name', f_name)

    % set parameter values for (cross)-bispectrum estimation
    freqinds = [mean(filt.low) mean(filt.high)]; % in Hz
    segleng = fs * g.epleng; 
    segshift = floor(segleng/2);
    epleng = fs * g.epleng; % create epochs of [e.epleng] seconds 

    % compute univariate bispectrum as a sanity check
    maxfreqbins = floor(segleng/2);
    bs_source = data2bs_univar(source, segleng, segshift, epleng, maxfreqbins);
%     bs_sens = data2bs_univar(signal_sensor(:,:)', segleng, segshift, epleng, maxfreqbins);
    plot_bispectra_univ(squeeze(abs(bs_source)), frqs, 0, jet, DIROUT, 'bispec_type', '_univ_source')
%     plot_bispectra_univ(squeeze(mean(abs(bs_sens), 1)), frqs, 0, jet, DIROUT, 'bispec_type', '_univ_sens')

    % run statistics on decomposition and compute FPR
    fpr = zeros(g.n, n_iter); 
    for i_iter = 1:n_iter
        P_fdr = bsfit_stats(signal_sensor, freqinds(1), freqinds(2), g.n, n_shuf, frqs, ...
            segleng, segshift, epleng, g.alpha, L);
        for i_source = 1:g.n
            fpr(i_source, i_iter) = squeeze(sum(P_fdr{i_source}(:) < g.alpha) ./ length(P_fdr{i_source}(:))); 
        end
    end

    % raincloud/half-violin plot on linearly scaled y-axis
    for i_source = 1:g.n
        titles = {'Normal', 'Train-Test Split'};
        colors = [[0 0 0.5]; [0.8 0 0.2]];
        plot_metrics_raincloud(fpr(i_source, :), colors, 1, 0.2, 'ks', titles, DIROUT, 'name', ['_n' int2str(i_source) f_name]);
    end
    
end
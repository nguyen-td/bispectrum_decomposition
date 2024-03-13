% Script to run a pre-analysis before the actual analysis involving the following steps:
%   1. Find peaks using FOOOF, first peak will be used as the f1 frequency.
%   2. Downsampling to 100 Hz.
%   3. Computation of univariate (surrogate) sensor bispectra and its p-values.
%   4. Computation of (surrogate) sensor cross-bispectra for (f1, f1, f1+f1) and (f1, f2, f1+f2).
% 
% Each compuational step is followed by various plotting steps. 
% 
% Inputs:
%   n_shuf - number of shuffles
%   isub  - index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   alpha    - significance level, default is 0.05.
%   poolsize - number of workers in the parellel pool (check parpool documentation) for parallel computing
%   f1       - fundamental frequency, used to compute cross-bispectra; if not passed, the fundamental frequency will be estimated via FOOOF
%   epleng   - length of epochs in seconds, default is 2 seconds

function main_preanalysis(n_shuf, isub, varargin)

    % set and add paths
%     DIROUT = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/';
    DIROUT = '/data/tdnguyen/git_repos/bispectrum_decomposition/Lemon/figures/';
    f_path = '/data/tdnguyen/data/lemon/data/';
%     f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';
%     f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/MotorImag/data/';
    
    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }              1;
        'f1'             'integer'       { }              0;
        'epleng'         'integer'       { }              2;
        });
    if ischar(g), error(g); end

    % load data
    sub = ['sub-032' num2str(isub)];
    f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data
%     sub = ['prep_vp' num2str(isub)];
%     f_name = [sub '.set']; % load motor imagery data
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);

    % find peaks using FOOOF and plot last fit if no f1 is passed
    if g.f1 == 0 % 
        [first_peak, ~, fooof_results] = find_peak_fooof(EEG);
        fooof_plot(fooof_results); xlabel('Frequency (Hz)'); ylabel('Power');
            exportgraphics(gcf, [DIROUT 'fooof_' lower(int2str(isub)) '.png'])
    else
        first_peak = g.f1;
%         second_peak = 2 * g.f1;
    end

    % plot channel locations
    figure; topoplot([], EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); 
    exportgraphics(gcf, [DIROUT, 'chanlocs_' num2str(isub) '.png']);

    % plot PSD and indicate the FOOOF peaks
    plot_spectra(EEG, 'EC', ['First peak: ' int2str(first_peak), 'Hz, Second peak: ' int2str(2 * first_peak) ' Hz'], DIROUT, ...
        'title_str', ['psd_' int2str(isub)], 'f1', first_peak)

    % downsample data to 100 Hz and plot
    EEG = downsampling(EEG, 100);
    plot_spectra(EEG, 'EC', ['First peak: ' int2str(first_peak), 'Hz, Second peak: ' int2str(2 * first_peak) ' Hz'], DIROUT, ...
        'title_str', ['psd_downsampled' int2str(isub)], 'f1', first_peak)

    % set parameter values for (cross)-bispectrum estimation
    data = EEG.data;
    segleng = EEG.srate * g.epleng; 
    segshift = floor(segleng/2);
    epleng = EEG.srate * g.epleng; % create epochs of [e.epleng] seconds 
    fres = EEG.srate;
    
    % compute univariate bicoherence and get the p-values
    frqs = sfreqs(fres, EEG.srate);
    [~, ~, P_sens_fdr_uni, ~, bispec, bicoh] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
%     save('P_sens.mat', 'P_sens_fdr')

    % get 2D positions of sensors and set up plotting
    locs_2D = create_locs_2D(EEG);
    
    % plot single matrices
    plot_pvalues_univ(P_sens_fdr_uni, frqs, isub, DIROUT)
    plot_bispectra_univ(bispec, frqs, isub, DIROUT, 'bispec_type', '_univ_unnorm', 'title_str', 'Unnormalized mean univariate sensor bispectrum') 
    plot_bispectra_univ(bicoh, frqs, isub, DIROUT, 'bispec_type', '_univ_norm', 'title_str', 'Normalized mean univariate sensor bispectrum')
    
    % plot matrices over the head
    clear para;
    para.tunit = 'Hz';
    para.funit = 'Hz';
    para.timeaxis = frqs(1:end-1);
    para.freqaxis = frqs(1:end-1);

    figure; showtfinhead(abs(bispec), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_unnormalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_normalized_' int2str(isub) '.png'])

    % set up computation of sensor cross-bispectra for (f1, f1, f1+f1) and (f1, f2, f1+f2)
    freqpairs1 = get_freqindices(round(first_peak), round(first_peak), frqs); % (f1, f1)
    freqpairs2 = get_freqindices(round(first_peak), 2 * round(first_peak), frqs); % (f1, f2)

    % estimate sensor cross-bispectrum
    clear bispec_para
    bispec_para.nrun = n_shuf;
    disp('Start calculating surrogate sensor cross-bispectra...')
    [bs_all1, bs_orig1, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs1, bispec_para);
    [bs_all2, bs_orig2, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs2, bispec_para);
    
    % compute normalization factor (threenorm)
    rtp1 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs1);
    rtp2 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs2);
    bicoh1 = bs_orig1 ./ rtp1;
    bicoh2 = bs_orig2 ./ rtp2;

    % compute and plot p-values for cross-bispectra
    [~, P_sens_fdr1] = compute_pvalues(mean(bs_all1(:, :, :, 1), 1), mean(bs_all1(:, :, :, 2:end), 1), n_shuf, g.alpha);
    [~, P_sens_fdr2] = compute_pvalues(mean(bs_all2(:, :, :, 1), 1), mean(bs_all2(:, :, :, 2:end), 1), n_shuf, g.alpha);
    plot_pvalues_univ(P_sens_fdr1, frqs, isub, DIROUT, 'bispec_type', '1_cross', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (cross-bispectrum)')
    plot_pvalues_univ(P_sens_fdr2, frqs, isub, DIROUT, 'bispec_type', '2_cross', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (cross-bispectrum)')
    
    % plot single matrices of net bispectra (collapsed over one channel dimension)
    plot_bispectra_univ(bs_orig1, frqs, isub, DIROUT, 'bispec_type', '1_cross_unnorm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net cross-bispectrum (f1, f1, f1+f1)') 
    plot_bispectra_univ(bicoh1, frqs, isub, DIROUT, 'bispec_type', '1_cross_norm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net cross-bispectrum (f1, f1, f1+f1)')
    plot_bispectra_univ(bs_orig2, frqs, isub, DIROUT, 'bispec_type', '2_cross_unnorm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net cross-bispectrum (f1, f2, f1+f2)') 
    plot_bispectra_univ(bicoh2, frqs, isub, DIROUT, 'bispec_type', '2_cross_norm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net cross-bispectrum (f1, f2, f1+f2)')
    
    % plot matrices over head
    clear plt_bispec_para
    plt_bispec_para.tunit = 'channel';
    plt_bispec_para.funit = 'channel';
    figure; showtfinhead(abs(bs_orig1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_unnorm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_norm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bs_orig2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_unnorm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_norm_sensor_head_' int2str(isub) '.png'])

    % plot cross-bispectra as topomaps with a seed
    net_bicoh1 = squeeze(mean(abs(bicoh1), 1));
    [seed1, seed1_idx] = max(max(net_bicoh1, [], 2), [], 'all'); % get max. value in the row (1st dimension)
    plot_topomaps_seed(net_bicoh1, seed1_idx, EEG.chanlocs, '1', 'Seed net cross-bicoherence (f1, f1, f1+f1)', DIROUT)

    net_bicoh2 = squeeze(mean(abs(bicoh2), 1));
    [seed2, seed2_idx] = max(max(net_bicoh1, [], 2), [], 'all'); % get max. value in the row (1st dimension)
    plot_topomaps_seed(net_bicoh2, seed2_idx, EEG.chanlocs, '2', 'Seed net cross-bicoherence (f1, f2, f1+f2)', DIROUT)

    % do we see between-site coupling? if yes, no non-sinusoids?

end
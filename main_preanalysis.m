% Script to run a pre-analysis before the actual analysis involving the following steps:
%   1. Find peaks using FOOOF, first peak will be used as the f1 frequency.
%   2. Downsampling to 100 Hz.
%   3. Computation of univariate (surrogate) sensor bispectra and its p-values.
%   4. Computation of (surrogate) sensor cross-bispectra for (f1, f1, f1+f1) and (f1, f2, f1+f2).
% 
% Each compuational step is followed by various plotting steps. 
% 
% Inputs:
%   n_shuf - [integer] number of shuffles
%   isub   - [integer] index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   alpha      - [float] significance level, default is 0.05.
%   poolsize   - [integer] number of workers in the parellel pool (check parpool documentation) for parallel computing
%   f1         - [integer] fundamental frequency, used to compute cross-bispectra; if not passed, the fundamental frequency will be estimated via FOOOF. Default is 0 (f1 will be estimated)
%   epleng     - [integer] length of epochs in seconds, default is 2 seconds
%   freq_down  - [integer] if 'downsample' is activated, the data will be downsampled to <freq_down> Hz. Default is 125 Hz.
%   downsample - [string] check whether to downsample data to <freq_down> Hz, default is 'on'

function main_preanalysis(n_shuf, isub, varargin)

    %% Prepare dataset, including downsampling to 125 Hz

    % set and add paths
%     DIROUT = ['/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/pre_' num2str(isub) '/'];
    DIROUT = ['/data/tdnguyen/git_repos/bispectrum_decomposition/Lemon/figures/pre_' num2str(isub) '/'];
    f_path = '/data/tdnguyen/data/lemon/data/';
%     f_path = '/Volumes/PortableSSD/LEMON/';
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
        'f1'             'float'       { }                0;
        'epleng'         'integer'       { }              2;
        'downsample'     'string'        {'on' 'off'}     'on';
        'freq_down'      'integer'       { }              125;
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
    plot_spectra(EEG, 'EC', ['First peak: ' num2str(first_peak), 'Hz, Second peak: ' num2str(2 * first_peak) ' Hz'], DIROUT, ...
        'title_str', ['psd_' num2str(isub)], 'f1', first_peak)

    % downsample data to 100 Hz and plot
    if strcmpi(g.downsample, 'on')
        EEG = downsampling(EEG, g.freq_down);
        plot_spectra(EEG, 'EC', ['First peak: ' num2str(first_peak), 'Hz, Second peak: ' num2str(2 * first_peak) ' Hz'], DIROUT, ...
            'title_str', ['psd_downsampled' num2str(isub)], 'f1', first_peak)
    end

    %% Compute univariate sensor bispectra

    % set parameter values for (cross)-bispectrum estimation
    data = EEG.data;
    segleng = EEG.srate * g.epleng; 
    segshift = floor(segleng/2);
    epleng = EEG.srate * g.epleng; % create epochs of [e.epleng] seconds 
    fres = EEG.srate;
    frqs = sfreqs(fres, EEG.srate);
    
    % compute univariate bicoherence and get the p-values
    [~, ~, P_sens_fdr_uni, ~, bispec, bicoh] = freq_preselection(data, n_shuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
%     save('P_sens.mat', 'P_sens_fdr')

    % get 2D positions of sensors and set up plotting
    locs_2D = create_locs_2D(EEG);
    load cm17.mat

    % plot single matrices
    p_cmap = cmap_pvalues(P_sens_fdr_uni, cm17, cm17a);
    plot_pvalues_univ(P_sens_fdr_uni, frqs, isub, p_cmap, DIROUT)
    plot_bispectra_univ(squeeze(mean(abs(bispec), 1)), frqs, isub, cm17a, DIROUT, 'bispec_type', '_univ_unnorm', 'title_str', 'Unnormalized mean univariate sensor bispectrum') 
    plot_bispectra_univ(squeeze(mean(abs(bicoh), 1)), frqs, isub, cm17a, DIROUT, 'bispec_type', '_univ_norm', 'title_str', 'Normalized mean univariate sensor bispectrum', 'isnorm', true)

    % plot matrices over the head
    clear para;
    para.tunit = 'Hz';
    para.funit = 'Hz';
    para.timeaxis = frqs(1:end-1);
    para.freqaxis = frqs(1:end-1);
    para.cmap = cm17a;

    figure; showtfinhead(abs(bispec), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_unnormalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_normalized_' int2str(isub) '.png'])

    %% Compute sensor cross-bispectra

    % set up computation of sensor cross-bispectra for (f1, f1, f1+f1) and (f1, f2, f1+f2)
    freqpairs1 = get_freqindices(round_to_05(first_peak), round_to_05(first_peak), frqs); % (f1, f1)
    freqpairs2 = get_freqindices(round_to_05(first_peak), 2 * round_to_05(first_peak), frqs); % (f1, f2)

    % estimate sensor cross-bispectrum
    clear bispec_para
    bispec_para.nrun = n_shuf;
    disp('Start calculating surrogate sensor cross-bispectra...')
    [bs_all1, bs_orig1, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs1, bispec_para);
    save([DIROUT 'bs_all1.mat'], 'bs_all1', '-v7.3')
    save([DIROUT 'bs_orig1.mat'], 'bs_orig1')
    [bs_all2, bs_orig2, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs2, bispec_para);
    
    % compute normalization factor (threenorm)
    rtp1 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs1);
    rtp2 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs2);
    bicoh1 = bs_orig1 ./ rtp1;
    bicoh2 = bs_orig2 ./ rtp2;

    for ichan = 1:3
        % compute and plot p-values for cross-bispectra
        [~, P_sens_fdr1] = compute_pvalues(squeeze(mean(abs(bs_orig1), ichan)), squeeze(mean(abs(bs_all1), ichan)), n_shuf, g.alpha);
        [~, P_sens_fdr2] = compute_pvalues(squeeze(mean(abs(bs_orig2), ichan)), squeeze(mean(abs(bs_all2), ichan)), n_shuf, g.alpha);
        p_cmap1 = cmap_pvalues(P_sens_fdr1, cm17, cm17a);
        p_cmap2 = cmap_pvalues(P_sens_fdr2, cm17, cm17a);
        plot_pvalues_univ(P_sens_fdr1, frqs, isub, p_cmap1, DIROUT, 'bispec_type', ['1_cross_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (f1, f1,  f1+f1)')
        plot_pvalues_univ(P_sens_fdr2, frqs, isub, p_cmap2, DIROUT, 'bispec_type', ['2_cross_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (f1,  f2, f1+f2)')

        % plot single matrices of net bispectra (collapsed over one channel dimension)
        plot_bispectra_univ(squeeze(mean(abs(bs_orig1), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '1_cross_unnorm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net cross-bispectrum (f1, f1, f1+f1)', 'mean_chan', ichan) 
        plot_bispectra_univ(squeeze(mean(abs(bicoh1), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '1_cross_norm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net cross-bispectrum (f1, f1, f1+f1)', 'mean_chan', ichan, 'isnorm', true)
        plot_bispectra_univ(squeeze(mean(abs(bs_orig2), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '2_cross_unnorm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net cross-bispectrum (f1, f2, f1+f2)', 'mean_chan', ichan) 
        plot_bispectra_univ(squeeze(mean(abs(bicoh2), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '2_cross_norm', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net cross-bispectrum (f1, f2, f1+f2)', 'mean_chan', ichan, 'isnorm', true)
    end

    % plot matrices over head
    clear plt_bispec_para
    plt_bispec_para.tunit = 'channel';
    plt_bispec_para.funit = 'channel';
    plt_bispec_para.cmap = cm17a;
    figure; showtfinhead(abs(bs_orig1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_unnorm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_norm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bs_orig2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_unnorm_sensor_head_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_norm_sensor_head_' int2str(isub) '.png'])

    % plot cross-bicoherence as topomaps with a seed
    for ichan = 1:3
        net_bicoh1 = squeeze(mean(abs(bicoh1), ichan));
        net_bicoh2 = squeeze(mean(abs(bicoh2), ichan));
        max_val = max([net_bicoh1, net_bicoh2], [], 'all');
        
        plot_topomaps_seed(net_bicoh1, EEG.chanlocs, max_val, cm17a, ['1_chan' int2str(ichan)], 'Seed net cross-bicoherence (f1, f1, f1+f1)', DIROUT, 'isnorm', true)
        plot_topomaps_seed(net_bicoh2, EEG.chanlocs, max_val, cm17a, ['2_chan' int2str(ichan)], 'Seed net cross-bicoherence (f1, f2, f1+f2)', DIROUT, 'isnorm', true)
    end

    %% Compute antisymmetrized cross-bispectra

    % compute antisymmetrized cross-bispectra
    bs_orig1_anti = bs_orig1 - permute(bs_orig1, [3, 2, 1]); % B_ijk - B_kji
    bicoh1_anti = bicoh1 - permute(bicoh1, [3, 2, 1]); % B_ijk - B_kji
    bs_all1_anti = bs_all1 - permute(bs_all1, [3, 2, 1, 4]); % B_ijk - B_kji
    bs_orig2_anti = bs_orig2 + permute(bs_orig2, [3, 1, 2]) + permute(bs_orig2, [2, 3, 1]) - permute(bs_orig2, [3, 2, 1]) - permute(bs_orig2, [2, 1, 3]) - permute(bs_orig2, [1, 3, 2]); % TACB
    bicoh2_anti = bicoh2 + permute(bicoh2, [3, 1, 2]) + permute(bicoh2, [2, 3, 1]) - permute(bicoh2, [3, 2, 1]) - permute(bicoh2, [2, 1, 3]) - permute(bicoh2, [1, 3, 2]); % TACB
    bs_all2_anti = bs_all2 + permute(bs_all2, [3, 1, 2, 4]) + permute(bs_all2, [2, 3, 1, 4]) - permute(bs_all2, [3, 2, 1, 4]) - permute(bs_all2, [2, 1, 3, 4]) - permute(bs_all2, [1, 3, 2, 4]); % TACB
    
    for ichan = 1:3
        % compute and plot p-values for antisymmetrized cross-bispectra
        [~, P_sens_anti_fdr1] = compute_pvalues(squeeze(mean(abs(bs_orig1_anti), ichan)), squeeze(mean(abs(bs_all1_anti), ichan)), n_shuf, g.alpha);
        [~, P_sens_anti_fdr2] = compute_pvalues(squeeze(mean(abs(bs_orig2_anti), ichan)), squeeze(mean(abs(bs_all2_anti), ichan)), n_shuf, g.alpha);
        p_cmap1 = cmap_pvalues(P_sens_anti_fdr1, cm17, cm17a);
        p_cmap2 = cmap_pvalues(P_sens_anti_fdr2, cm17, cm17a);
        plot_pvalues_univ(P_sens_anti_fdr1, frqs, isub, p_cmap1, DIROUT, 'bispec_type', ['1_cross_anti_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (f1, f1,  f1+f1)')
        plot_pvalues_univ(P_sens_anti_fdr2, frqs, isub, p_cmap2, DIROUT, 'bispec_type', ['2_cross_anti_chan' int2str(ichan)], 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'p-values (f1,  f2, f1+f2)')
        
        % plot matrices of net antisymmetrized bispectra (collapsed over one channel dimension)
        plot_bispectra_univ(squeeze(mean(abs(bs_orig1_anti), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '1_cross_unnorm_anti', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net antisymmetrized cross-bispectrum (f1, f1, f1+f1)', 'mean_chan', ichan) 
        plot_bispectra_univ(squeeze(mean(abs(bicoh1_anti), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '1_cross_norm_anti', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net antisymmetrized cross-bispectrum (f1, f1, f1+f1)', 'mean_chan', ichan, 'isnorm', true) 
        plot_bispectra_univ(squeeze(mean(abs(bs_orig2_anti), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '2_cross_unnorm_anti', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Unnormalized net totally antisymmetrized cross-bispectrum (f1, f2, f1+f2)', 'mean_chan', ichan)
        plot_bispectra_univ(squeeze(mean(abs(bicoh2_anti), ichan)), frqs, isub, cm17a, DIROUT, 'bispec_type', '2_cross_norm_anti', 'label_x', 'channel', 'label_y', 'channel', 'custom_label', 0, 'title_str', 'Normalized net totally antisymmetrized cross-bispectrum (f1, f2, f1+f2)', 'mean_chan', ichan, 'isnorm', true) 
    
        % plot cross-biscoherence as topomaps with a seed
        net_bicoh1_anti = squeeze(mean(abs(bicoh1_anti), ichan));
        net_bicoh2_anti = squeeze(mean(abs(bicoh2_anti), ichan));
        max_val_anti = max([net_bicoh1_anti, net_bicoh2_anti], [], 'all');
    
        plot_topomaps_seed(net_bicoh1_anti, EEG.chanlocs, max_val_anti, cm17a, ['1_anti_chan' int2str(ichan)], 'Seed net antisymmetrized cross-bicoherence (f1, f1, f1+f1)', DIROUT, 'isnorm', true)
        plot_topomaps_seed(net_bicoh2_anti, EEG.chanlocs, max_val_anti, cm17a, ['2_anti_chan' int2str(ichan)], 'Seed net antisymmetrized cross-bicoherence (f1, f2, f1+f2)', DIROUT, 'isnorm', true)
    end
end
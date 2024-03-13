% Script to run a pre-analysis before the actual analysis. Here, we first
% determine the alpha frequency and its first higher harmonic. Then, we
% compute the univariate sensor bispectrum and see if we find genuine 
% between-site interactions or within-site interactions (potentially spurious interactions, 
% possibly arising from non-sinusoidal signals, leading to coupling of higher harmonics). 
% 
% Inputs:
%   nshuf - number of shuffles
%   isub  - index of subject (the pipeline works for a single subject)
%
% Optional inputs:
%   alpha    - significance level, default is 0.05.
%   poolsize - number of workers in the parellel pool (check parpool documentation) for parallel computing
%   f1       - fundamental frequency, used to compute cross-bispectra; if not passed, the fundamental frequency will be estimated via FOOOF
%   epleng   - length of epochs in seconds, default is 2 seconds

function main_preanalysis(nshuf, isub, varargin)

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
    [~, ~, P_sens_fdr, ~, bispec, bicoh] = freq_preselection(data, nshuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);
    save('P_sens.mat', 'P_sens_fdr')

    % get 2D positions of sensors and set up plotting
    locs_2D = create_locs_2D(EEG);
%     max_freq = 60; % in Hz
%     freq_res = frqs(2) - frqs(1);
%     freq_inds = 1:max_freq * 1/freq_res;
    
    % plot single matrices
    plot_pvalues_univ(frqs, isub, DIROUT, P_sens_fdr)
%     plot_bispectra_univ(frqs, bispec(:, freq_inds, freq_inds), isub, 'Unnormalized mean', DIROUT) % plot bispectra/bicoherence between 0 - 60 Hz
%     plot_bispectra_univ(frqs, bicoh(:, freq_inds, freq_inds), isub, 'Normalized mean', DIROUT)
    plot_bispectra_univ(frqs, bispec, isub, 'Unnormalized mean', DIROUT) % plot bispectra/bicoherence between 0 - 60 Hz
    plot_bispectra_univ(frqs, bicoh, isub, 'Normalized mean', DIROUT)
    
    
    % plot matrices over the head
    clear para;
    para.tunit = 'Hz';
    para.funit = 'Hz';
%     para.timeaxis = frqs(freq_inds);
%     para.freqaxis = frqs(freq_inds);
    para.timeaxis = frqs(1:end-1);
    para.freqaxis = frqs(1:end-1);

%     figure; showtfinhead(abs(bispec(:, freq_inds, freq_inds)), locs_2D, para); 
%         exportgraphics(gcf, [DIROUT 'B_sensor_head_unnormalized_' int2str(isub) '.png'])
%     figure; showtfinhead(abs(bicoh(:, freq_inds, freq_inds)), locs_2D, para); 
%         exportgraphics(gcf, [DIROUT 'B_sensor_head_normalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bispec), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_unnormalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'B_sensor_head_normalized_' int2str(isub) '.png'])

    % set up computation of sensor cross-bispectra for (f1, f1, f1+f1) and (f1, f2, f1+f2)
    freqpairs1 = get_freqindices(round(first_peak), round(first_peak), frqs); % (f1, f1)
    freqpairs2 = get_freqindices(round(first_peak), 2 * round(first_peak), frqs); % (f1, f2)

    % estimate sensor cross-bispectrum
    clear bispec_para
    bispec_para.nrun = nshuf;
    disp('Start calculating surrogate sensor cross-bispectra...')
    [bs_all1, bs_orig1, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs1, bispec_para);
    [bs_all2, bs_orig2, ~] = data2bs_event_surro_final(data(:, :)', segleng, segshift, epleng, freqpairs2, bispec_para);
    
    % compute normalization factor (threenorm)
    rtp1 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs1);
    rtp2 = data2bs_threenorm(data(:, :)', segleng, segshift, epleng, freqpairs2);
    bicoh1 = bs_orig1 ./ rtp1;
    bicoh2 = bs_orig2 ./ rtp2;
    
    % plot single matrices of net bispectra (collapsed over one channel dimension)
    plot_bispectra_univ(frqs, bs_orig1, isub, 'Unnormalized net cross-bispectrum (f1, f1, f1+f1)', DIROUT, 'bispec_type', '', 'label_x', 'channel', 'label_y', 'channel') 
    plot_bispectra_univ(frqs, bicoh1, isub, 'Normalized net cross-bispectrum (f1, f1, f1+f1)', DIROUT, 'bispec_type', '', 'label_x', 'channel', 'label_y', 'channel')
    plot_bispectra_univ(frqs, bs_orig2, isub, 'Unnormalized net cross-bispectrum (f1, f2, f1+f2)', DIROUT, 'bispec_type', '', 'label_x', 'channel', 'label_y', 'channel') 
    plot_bispectra_univ(frqs, bicoh2, isub, 'Normalized net cross-bispectrum (f1, f2, f1+f2)', DIROUT, 'bispec_type', '', 'label_x', 'channel', 'label_y', 'channel')
    
    % plot matrices over head
    clear plt_bispec_para
    plt_bispec_para.tunit = 'channel';
    plt_bispec_para.funit = 'channel';
    figure; showtfinhead(abs(bs_orig1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_sensor_head_unnormalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bs_orig2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_sensor_head_unnormalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh1), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B1_cross_sensor_head_normalized_' int2str(isub) '.png'])
    figure; showtfinhead(abs(bicoh2), locs_2D, plt_bispec_para); 
        exportgraphics(gcf, [DIROUT 'B2_cross_sensor_head_normalized_' int2str(isub) '.png'])
    % do we see between-site coupling? if yes, no non-sinusoids?

end
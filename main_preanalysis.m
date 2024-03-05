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
%   poolsize - number of workers in the parellel pool (check parpool documentation) for parallel computing

function main_preanalysis(nshuf, isub, varargin)

    % set and add paths
    DIROUT = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/figures/';
    f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';
    addpath('/Users/nguyentiendung/GitHub/eeglab')
    addpath(genpath('/Users/nguyentiendung/GitHub/fooof_mat/'))
    
    if ~exist(DIROUT, 'dir')
        mkdir(DIROUT)
    end

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end

    % load data
    sub = ['sub-032' num2str(isub)];
    f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data
    
    % load preprocessed EEG
    EEG = pop_loadset(f_name, f_path);

    % set parameter values for (cross)-bispectrum estimation
    data = EEG.data;
    segleng = EEG.pnts;
    segshift = floor(segleng/2);
    epleng = EEG.pnts; 
    fres = EEG.srate;
    
    % compute univariate bicoherence and get the p-values
    frqs = sfreqs(fres, EEG.srate);
    [~, ~, P_sens_fdr, ~, bispec, bicoh] = freq_preselection(data, nshuf, frqs, segleng, segshift, epleng, g.alpha, g.poolsize);

    % get 2D positions of sensors
    locs_2D = create_locs_2D(EEG);

    % plotting
    max_freq = 60; % in Hz
    freq_res = frqs(2) - frqs(1);
    freq_inds = 1:max_freq * 1/freq_res;
    
    % plot single matrices
    plot_pvalues_univ(frqs, isub, DIROUT, P_sens_fdr)
    plot_bispectra_univ(frqs, bispec(:, freq_inds, freq_inds), isub, 'Unnormalized mean', DIROUT) % plot bispectra/bicoherence between 0 - 60 Hz
    plot_bispectra_univ(frqs, bicoh(:, freq_inds, freq_inds), isub, 'Normalized mean', DIROUT)
    
    % plot matrices over the head
    clear para;
    para.tunit = 'Hz';
    para.funit = 'Hz';
    para.timeaxis = frqs(freq_inds);
    para.freqaxis = frqs(freq_inds);
    figure; showtfinhead(abs(bispec(:, freq_inds, freq_inds)), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'Unnormalized bispectrum_' lower(int2str(isub)) '.png'])
    figure; showtfinhead(abs(bicoh(:, freq_inds, freq_inds)), locs_2D, para); 
        exportgraphics(gcf, [DIROUT 'Normalized bispectrum_' lower(int2str(isub)) '.png'])
        
    % find peaks using FOOOF
    [first_peak, second_peak] = find_peak_fooof(EEG);
end
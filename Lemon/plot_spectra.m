function plot_spectra(EEG, condition, listing, result_dir, isub)
    EEG = eeg_checkset( EEG );
    figure; pop_spectopo(EEG, 1, [0 EEG.xmax * 1000], 'EEG' , 'percent', 100, 'freq', [10 20 30], 'freqrange',[0 EEG.srate/2], ...
        'electrodes', 'on', 'overlap', 50, 'title', [condition, ' ', listing(isub).name]);
    exportgraphics(gcf, [result_dir, listing(isub).name, '_', condition, '_psd.png']);
end
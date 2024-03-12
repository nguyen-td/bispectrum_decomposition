% Downsample dataset from [EEG.srate] Hz to [fs_new] Hz, where [EEG.srate] > [fs_new].
%
% Inputs:
%   EEG - EEG struct, where EEG.data has size (n_chans x n_pnts)
%   fs_new - new frequency in Hz

function EEG = downsampling(EEG, fs_new)

    % low-pass filter to avoid aliasing
    cutoff_freq = fs_new/2 * 0.9;
    [b_low, a_low] = butter(2, cutoff_freq/(EEG.srate/2), 'low');
    filtdata = filtfilt(b_low, a_low, double(EEG.data));
    
    % simple downsampling
    dsRatio = EEG.srate / fs_new;
    EEG.data = filtdata(:, 1:dsRatio:end);
    EEG.srate = fs_new;
    EEG.pnts = size(EEG.data, 2);
    EEG.xmax = EEG.xmin + (EEG.pnts-1)/EEG.srate;
    EEG.times = linspace(EEG.xmin * 1000, EEG.xmax * 1000, EEG.pnts);
    EEG = eeg_checkset( EEG );
end
% Downsample dataset from [EEG.srate] Hz to [fs_new] Hz, where [EEG.srate] > fs_new].
%
% Inputs:
%   EEG - EEG struct, where EEG.data has size (n_chans x n_pnts)
%   fs_new - new frequency in Hz

function EEG = downsampling(EEG, fs_new)
    dsRatio = EEG.srate / fs_new;
    EEG.data = EEG.data(:, 1:dsRatio:end);
%     EEG.times = EEG.times(1:dsRatio:end);
    EEG.srate = fs_new;
    EEG = eeg_checkset( EEG );
end
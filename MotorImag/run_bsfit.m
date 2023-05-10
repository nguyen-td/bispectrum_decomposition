% Compute the bispectrum on preprocessed motor imagery data, then fit the
% low-rank model on the estimated bispectrum. 
%% Add paths
addpath(genpath('/Volumes/Seagate/MotorImag/data/'))
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab')
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition'))
eeglab

%% Load data
isub = 3; % pick a single subject
sub = ['vp' num2str(isub)];
f_name = ['prep_' sub '.set'];
f_path = '/Volumes/Seagate/MotorImag/data/';
EEG = pop_loadset(f_name, f_path);

%% Set inputs for bispectrum estimation
data = EEG.data;
segleng = EEG.srate;
segshift = floor(segleng/2);
epleng = 2 * EEG.srate; 

%% Calculate bispectrum
% choose frequency pair in Hz
f1 = 11; % as a start, choose 11 Hz (mu rhythm) and 22 Hz (beta rhythm) bc it was mentioned a couple of times in Sanelli et al., 2019
f2 = 22;
freqpairs = [f1,f2];

[bs_orig, ~] = data2bs_event(data', segleng, segshift, epleng, freqpairs);

%% Make low dimensional fit
n = 3; % number of sources  
[a, d, err, err_all, bsmodel] = bsfit(bs_orig, n);


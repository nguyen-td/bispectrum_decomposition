% Compute the cross-bispectrum on preprocessed motor imagery data, then fit the
% low-rank model on the estimated cross-bispectrum. 
%% Add paths and load stuff
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1'
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab')
eeglab
addpath(genpath('/Volumes/PortableSSD/MotorImag/data/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition'))

% load cortex
% load cm17

%% Load data
isub = 3; % pick a single subject
sub = ['vp' num2str(isub)];
f_name = ['prep_' sub '.set'];
f_path = '/Volumes/PortableSSD/MotorImag/data/';
EEG = pop_loadset(f_name, f_path);

%% Set parameter values for cross-bispectrum estimation
data = EEG.data;
segleng = EEG.srate;
segshift = floor(segleng/2);
epleng = EEG.srate; 

%% Test significance of the fitted source cross-bispectrum within subjects
f1 = 11;
f2 = 22;
n = 3;
nshuf = 3;
fres = EEG.srate;

[bs_all, bs_orig] = run_bsfit(data, f1, f2, n, nshuf, fres, EEG.srate, segleng, segshift, epleng);

%% Make topoplot
% first perform source reconstruction to localize mixing matrix A and then
% make topoplot I think
min_val = min(a, [], 'all');
max_val = max(a, [], 'all');
for i = 1:n
    figure;
    title(sprintf('Source %d', i))
    topoplot(a(:, i), EEG.chanlocs)
    colorbar; clim([min_val max_val])
end

%% Statistics
% 

%% Sanity check
p = 10;
q = 12;
r = 15;
bs_est_pqr = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            bs_est_pqr = bs_est_pqr + (a(p,i) * a(q,j) * a(r,k) * d(i,j,k));
        end
    end
end
disp(bs_est_pqr)

% tensor notation
bs_est = tensorprod(a, tensorprod(a, tensorprod(a, d, 2, 3), 2, 3), 2, 3);
figure; imagesc(squeeze(abs(bs_orig(:,:,5))))
figure; imagesc(squeeze(abs(bs_est(:,:,5))))
disp(bs_est(p,q,r))

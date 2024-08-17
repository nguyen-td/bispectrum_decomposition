% Script to test the method extension to decompose over frequency bands.

%% Setup
clear

f_path = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/Lemon/data/';
DIROUT = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/tests/figures/f1_f2/';
isub = 317;
freq_down = 125; % new sampling frequency after downsampling (Hz)
n_shuf = 2;

%% Load data
eeglab
sub = ['sub-032' num2str(isub)];
f_name = [sub '/' sub '_EC.set']; % load LEMON eyes-closed data

EEG = pop_loadset(f_name, f_path);
EEG = downsampling(EEG, freq_down);
fres = EEG.srate;
frqs = sfreqs(fres, EEG.srate);

%% Set parameter values for (cross)-bispectrum estimation
data = EEG.data;
segleng = EEG.pnts;
segshift = floor(segleng/2);
epleng = EEG.pnts; 

%% Set frequency pairs
f1_a = 9.5;
% f2_a = f1_a * 2;
f2_a = f1_a;

f1_b = [8 10];
% f2_b = f1_b * 2;
f2_b = f1_b;

f1_c = [9 10];
% f2_c = f1_c * 2;
f2_c = f1_c;

f1_d = [9.5 10];
% f2_d = f1_d * 2;
f2_d = f1_d;

%% Get frequency bins
freqpairs_a = get_freqindices(round_to_05(f1_a), round_to_05(f2_a), frqs); 
freqpairs_b = get_freqindices(round_to_05(f1_b), round_to_05(f2_b), frqs); 
freqpairs_c = get_freqindices(round_to_05(f1_c), round_to_05(f2_c), frqs); 
freqpairs_d = get_freqindices(round_to_05(f1_d), round_to_05(f2_d), frqs); 

%% Estimate sensor cross-bispectrum
bs_orig_a = data2bs_event(data(:, :)', segleng, segshift, epleng, freqpairs_a);
bs_orig_b = data2bs_event(data(:, :)', segleng, segshift, epleng, freqpairs_b);
bs_orig_c = data2bs_event(data(:, :)', segleng, segshift, epleng, freqpairs_c);
bs_orig_d = data2bs_event(data(:, :)', segleng, segshift, epleng, freqpairs_d);

%% Sanity check if bsfit.m and bsfit_freqbands.m give the same results for single frequencies
n = 3;
[A_a, D_a] = bsfit(bs_orig_a, n);
[A_a1, D_a1] = bsfit_freqbands(bs_orig_a, n);
if isequal(A_a, A_a1) && isequal(D_a, D_a1)
    disp('Both functions return equal results for single frequencies.')
end

%% Fit model for all other frequency (pair) combinations with different model orders
n = [3, 5, 8, 10];
errors = {};

for ind = 1:length(n)
    tic
    disp(['Fitting model with order ' int2str(n(ind))])
    [~, ~, ~, ~, err1] = bsfit_freqbands(bs_orig_a, n(ind));
    [~, ~, ~, ~, err2] = bsfit_freqbands(bs_orig_b, n(ind));
    [~, ~, ~, ~, err3] = bsfit_freqbands(bs_orig_c, n(ind));
    [~, ~, ~, ~, err4] = bsfit_freqbands(bs_orig_d, n(ind));

    errors{1, ind} = err1;
    errors{2, ind} = err2;
    errors{3, ind} = err3;
    errors{4, ind} = err4;
    toc
end
% save('tests/errors.mat', 'errors')

%% Plot results
colors = ['r', 'b', 'k', 'g'];

plot_error(errors, 1, n, colors, ['f1 = ' num2str(f1_a) ' Hz , f2 = ' num2str(f2_a) ' Hz'], '', DIROUT, 'f_name', '_a', 'f_ext', '.fig');
plot_error(errors, 2, n, colors, ['f1 = ' num2str(f1_b(1)) '-' num2str(f1_b(2)) ' Hz , f2 = ' num2str(f2_b(1)) '-' num2str(f2_b(2)) ' Hz'], '', DIROUT, 'f_name', '_b', 'f_ext', '.fig'); 
plot_error(errors, 3, n, colors, ['f1 = ' num2str(f1_c(1)) '-' num2str(f1_c(2)) ' Hz , f2 = ' num2str(f2_c(1)) '-' num2str(f2_c(2)) ' Hz'], '', DIROUT, 'f_name', '_c', 'f_ext', '.fig'); 
plot_error(errors, 4, n, colors, ['f1 = ' num2str(f1_d(1)) '-' num2str(f1_d(2)) ' Hz , f2 = ' num2str(f2_d(1)) '-' num2str(f2_d(2)) ' Hz'], '', DIROUT, 'f_name', '_d', 'f_ext', '.fig'); 

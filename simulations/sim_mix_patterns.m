% Simulate mixed patterns and demix with MOCA to determine correct demixing.

%% Set output filename and paths
clear
DIROUT = '/Users/nguyentiendung/GitHub/bispectrum_decomposition/simulations/figures/';
addpath(genpath('/Users/nguyentiendung/GitHub/bispectrum_decomposition/'))
addpath('/Users/nguyentiendung/GitHub/eeglab')
if ~exist(DIROUT, 'dir')
   mkdir(DIROUT)
end

%% Load the leadfield
try
    load sa_nyhead
catch
    warning("Please download the leadfield first: https://www.parralab.org/nyhead/ and make sure the file exists.")
end
% L_3D = sa.cortex75K.V_fem(:, sa.cortex2K.in_from_cortex75K, :);
load L_3D % 90 channels

%% Select spatially distributed patterns at n voxels randomly
n = 4;
rng(100)
n_voxels = size(L_3D, 2);
n_chans = size(L_3D, 1);

vox_ind = randi(n_voxels, n, 1);
B_3D = L_3D(:, vox_ind, :); 
B = zeros(n_chans, n);
for b = 1:n
    ori = randn(3, 1);
    ori = ori / norm(ori);
    B(:, b) = squeeze(B_3D(:, b, :)) * ori;
end

%% Create mixed patterns
M = randn(n, n);
A = B * M;

%% Create sources and demix them with MOCA
[A_moca, F_moca, F] = apply_moca(L_3D, A, n);

%% Plot and save mixed and demixed sources
load cm17
plot_sources(F_moca, n, sa.cortex75K, sa.cortex2K, vox_ind, cm17a, '', DIROUT)

%% Plot and save topomaps of mixed and demixed patterns
load chanlocs 
eeglab
plot_topomaps_patterns(B, n, chanlocs, cm17, '', 'original', DIROUT) 
plot_topomaps_patterns(A, n, chanlocs, cm17, '', 'mixed', DIROUT)

B_hat = A * A_moca'; % demixed patterns
plot_topomaps_patterns(B_hat, n, chanlocs, cm17, '', 'demixed', DIROUT) 

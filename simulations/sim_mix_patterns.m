%% Simulate mixed patterns and demix with MOCA to determine correct demixing.

% Set output filename and paths
DIROUT = '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/simulations/figures';
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab')
if ~exist(DIROUT, 'dir')
   mkdir(DIROUT)
end

% Load the leadfield
try
    load sa_nyhead
catch
    warning("Please download the leadfield first: https://www.parralab.org/nyhead/ and make sure the file exists.")
end
% L_3D = sa.cortex75K.V_fem(:, sa.cortex2K.in_from_cortex75K, :);
load L_3D % 90 channels

% Select spatially distributed patterns at n voxels randomly
n = 4;
rng(42)
vox_ind = randi(1006, n, 1);
B = L_3D(:, vox_ind, 1);

% Create mixed patterns
M = randn(n, n);
A = B * M;

% Create sources and demix them with MOCA
[A_moca, F_moca, F] = apply_moca(L_3D, A, n);

% Plot and save mixed and demixed sources
load cm17
plot_sources(F_moca, F, n, sa.cortex75K, sa.cortex2K, vox_ind, cm17, '', DIROUT)

% Plot and save topomaps of mixed and demixed patterns
load chanlocs 
eeglab
plot_topomaps(B, n, chanlocs, 'original', DIROUT) 
plot_topomaps(A, n, chanlocs, 'mixed', DIROUT)

B_hat = A * A_moca'; % demixed patterns
plot_topomaps(B_hat, n, chanlocs, 'demixed', DIROUT) 

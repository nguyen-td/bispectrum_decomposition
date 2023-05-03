%% Set up directories and load classification scores
DIRIN = '/Volumes/Seagate/VitalBCI/mat/imag/';
load('/Volumes/Seagate/VitalBCI/mat/scores.mat')

%% Load subjects
nsub = 40;

sub = ['vp' num2str(2)];
load([DIRIN sub '.mat'])
plot(cnt.x')
for isub = 2:nsub
    sub = ['vp' num2str(isub)];
    disp(sub)
end
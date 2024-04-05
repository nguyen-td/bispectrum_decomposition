% Script to plot PSDs for the LEMON dataset
%% Set paths and directories
clear

addpath('/Volumes/PortableSSD/LEMON/')
addpath('/Users/nguyentiendung/GitHub/eeglab')
data_path = '/Volumes/PortableSSD/LEMON/';
result_dir = '/Volumes/PortableSSD/LEMON/figures/PSDs/';

%% Settings
eeglab

listing = dir([data_path, 'sub*']);

%% Load data and plot PSDs
for isub = 1:length(listing)
    try
        disp(['Subject: ' listing(isub).name])
        % load eyes-closed (EC) and eyes-open (EO) files
        f_name = [listing(isub).folder, '/', listing(isub).name, '/'];
        EEG_EC = pop_loadset([f_name, listing(isub).name, '_EC.set']);
%         EEG_EO = pop_loadset([f_name, listing(isub).name, '_EO.set']);

        % plot PSDs and save them as .png files
        plot_spectra(EEG_EC, 'EC', listing(isub).name, result_dir)
%         plot_spectra(EEG_EO, 'EO', listing(isub).name, result_dir)
    catch
        continue
    end
end

%% Make presentation
import mlreportgen.ppt.*
ppt = Presentation(strcat(result_dir, 'power_spectra.pptx'));

EC_files = dir([result_dir, 'sub*EC*.png']);
% EO_files = dir([result_dir, 'sub*EO*.png']);

make_slides(EC_files, result_dir, ppt);
% make_slides(EO_files, result_dir, ppt);

close(ppt)
%rptview(ppt) % uncomment to open .ppt file
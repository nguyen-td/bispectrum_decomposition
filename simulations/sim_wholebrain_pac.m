% Script to simulate whole-brain signals with univariate and bivariate PAC interactions. 
% Based on Franziska Pellegrini's PAC repository https://github.com/fpellegrini/PAC/blob/master/fp_pac_sim.m.
%
% Inputs:
%   sim_case - [1, 2, 3] univariate, bivariate or univariate + bivariate
%   n_univ   - [integer] number of univariate interactions, default is 1
%   n_biv    - [integer] number of bivariate interactions, default is 2
%   isnr     - [float] total signal-to-noise ratio
%
% Output:
%   signal_roi - (n_roi x epleng x n_epochs) simulated ROI data
%   fs         - [integer] sampling frequency
%   source     - ([epleng * n_epochs] x [n_univ + n_biv]) source data containing source PAC
%   filt       - filter settings

function [signal_sensor, fs, sources, filt, L] = sim_wholebrain_pac(sim_case, n_univ, n_biv, isnr)

    %% Signal generation

    % set parameters (see fp_pac_signal.m for documentation)
    params.case = sim_case; % univariate + bivariate case
    params.iReg = 1; 
    params.iss = 0.9; % brain noise to sensor noise ration (19 dB)
    params.isnr = isnr; % total SNR 
    params.t = 0; % [0, 1] true voxel pipeline or not
    if params.case == 3
        assert(~any([n_univ n_biv] == 0), 'Indicate the number of uni and bivariate interactions in the mixed case.')
        params.iInt = [n_univ n_biv]; % one univariate, two bivariate interactions
    else
        params.iInt = sum([n_univ n_biv]); 
    end

    % get atlas, voxel and roi indices; active voxel of each region is aleady selected here
    iReg = 1; 
    fprintf('Getting atlas positions... \n')
    D = fp_get_Desikan(iReg); 

    % signal generation
    fprintf('Signal generation... \n')
    [sig, brain_noise, sensor_noise, L, iroi_phase, iroi_amplt, D, fs, n_trials, filt, sources] = fp_pac_signal(params, D);

    % combine noise sources
    noise = params.iss*brain_noise + (1-params.iss)*sensor_noise;
    noise = noise ./ norm(noise(:),'fro');

    % combine signal and noise
    signal_sensor1 = params.isnr*sig + (1-params.isnr)*noise;
    signal_sensor1 = signal_sensor1 ./ norm(signal_sensor1(:), 'fro');

    % high-pass filter signal
    signal_sensor = (filtfilt(filt.bhigh, filt.ahigh, signal_sensor1'))';
    signal_sensor = signal_sensor / norm(signal_sensor, 'fro');
    
    % reshape (create epochs)
    signal_sensor = reshape(signal_sensor,[],size(signal_sensor,2)/n_trials,n_trials);
    [n_sensors, l_epoch, n_trials] = size(signal_sensor);

    %% Source reconstruction

    % select only voxels that belong to any roi
    L_backward = L(:, D.ind_cortex, :);
    A = fp_get_lcmv(signal_sensor, L_backward);

    % dimensionality reduction
    signal_roi = fp_dimred(signal_sensor,D,A,params.t);
end
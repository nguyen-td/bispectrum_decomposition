% Pipeline to run the decomposition method on simualted PAC data with a
% specified number of univariate and bivariate interactions.
%
% Input:
%   n_shuf - [integer] number of shuffles
%
% Optional inputs:
%   n      - [integer] model order/number of fitted sources, default is 5. Can also be an array of n's, e.g., [3 4 5]
%   n_univ - [integer] number of univariate interactions, default is 1
%   n_biv  - [integer] number of bivariate interactions, default is 2

function main_sim_pac(n_shuf, varargin)

    % setup
    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'n_univ'         'integer'       { }                1;
        'n_biv'          'integer'       { }                2;
        });
    if ischar(g), error(g); end
    
    % check case
    if ~g.n_univ == 0 && g.n_biv == 0
        sim_case = 1;
    elseif g.n_univ == 0 && ~g.n_biv == 0
        sim_case = 2;
    else
        sim_case = 3;
    end
    
    % generate simulated data
    [signal_sensor, fs, source] = sim_wholebrain_pac(sim_case, g.n_univ, g.n_biv);
    psd = pwelch(signal_sensor(:,:)', 100, 50, 2*fs, fs);

    % plot sources
end
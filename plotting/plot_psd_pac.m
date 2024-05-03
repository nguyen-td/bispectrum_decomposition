% Plot and save PSD of simulated PAC.
%
% Inputs:
%   psd    - (n_freqs x n_chans) power spectrum
%   frqs   - (n_freqs x 1) vector of frequency bins
%   DIROUT - output directory to save images
%
% Optional inputs:
%   name - [string] file name, default is ''

function plot_psd_pac(psd, frqs, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'name'    'string'     { }     '';
        });
    if ischar(g), error(g); end

    plot(frqs, 10 * log10(psd))
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    set(gca, 'FontSize', 15)
    grid on;
    
    % save figure
    save_err = [DIROUT 'PSD_sim_' g.name '.png'];
    exportgraphics(gcf, save_err)
end
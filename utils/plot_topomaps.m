% Create topomaps of patterns. Based on code from plot_sources.m
%
% Inputs:
%   A        - (n_chan x n) mixing matrix
%   n        - number of estimated sources
%   chanlocs - EEG channel locations from the EEG struct
%   name     - [string] name for plot, for example, "mixed", "demixed" etc.
%   DIROUT   - output directory to save images

function plot_topomaps(A, n, chanlocs, name, DIROUT)

    % topoplots of mixing matrix
    figure('Position', [600 100 1500 300]);
    tl2 = tiledlayout(1, n);
    for i = 1:n
        nexttile
        t = title(i);
        t.FontSize = 20;
        if i == n
            topoplot(A(:, i), chanlocs, 'electrodes', 'on'); colorbar; clim([min(A, [], 'all') max(A, [], 'all')])
        else
            topoplot(A(:, i), chanlocs, 'electrodes', 'on'); clim([min(A, [], 'all') max(A, [], 'all')])
        end
    end
    save_A = [DIROUT '/A_' name '.png'];
    exportgraphics(gcf, save_A)
    
end
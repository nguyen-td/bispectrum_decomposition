% Create topomaps of patterns. 
%
% Inputs:
%   A        - (n_chan x n) mixing matrix
%   n        - number of estimated sources
%   chanlocs - EEG channel locations from the EEG struct
%   cmap     - (n, 3) colormap 
%   isub     - subject ID, used for title
%   name     - [string] name for plot, for example, "mixed", "demixed" etc.
%   DIROUT   - output directory to save images
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is '_cross' (empty string)

function plot_topomaps_patterns(A, n, chanlocs, cmap, isub, name, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '';
        });
    if ischar(g), error(g); end

    % topoplots of mixing matrix
    figure('Position', [600 100 1500 300]);
    tl2 = tiledlayout(1, n);
    for i = 1:n
        nexttile
        t = title(i);
        t.FontSize = 20;
        if i == n
            topoplot(A(:, i), chanlocs, 'electrodes', 'on', 'colormap', cmap); colorbar; clim([min(A, [], 'all') max(A, [], 'all')])
        else
            topoplot(A(:, i), chanlocs, 'electrodes', 'on', 'colormap', cmap); clim([min(A, [], 'all') max(A, [], 'all')])
        end
    end
    save_A = [DIROUT 'A' g.bispec_type '_' name '_' int2str(isub) '.png'];
    exportgraphics(gcf, save_A)
    
end
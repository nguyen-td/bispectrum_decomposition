% Plot and save absolute values of the cross-bispectra.
%
% Inputs:
%   D        - (n x n x n) source cross-bispectrum
%   f1, f2   - frequencies in Hz, used for title
%   isub     - subject ID, used for title
%   name     - [string] name for plot, for example, "mixed", "demixed" etc.
%   DIROUT   - output directory to save images
%   cmap     - (n, 3) colormap 
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is '' 
%   istitle       - [boolean] whether to show a title, default is true
%   f_ext         - [string] file extension, default is .png.

function plot_bispectra(D, f1, f2, isub, name, DIROUT, cmap, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '';
        'istitle'        'boolean'    { }     true;
        'f_ext'          'string'     { }     '.png';
        });
    if ischar(g), error(g); end
    
    D_abs = abs(D);
    n = size(D, 1);
    if n == 3
        figure('Position', [600 100 1000 300]);
    else
        figure('Position', [600 100 1500 300]);
    end
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(squeeze(D_abs(i, :, :)));
        clim([0 max(D_abs, [], 'all')])
        if i == n
            c = colorbar;
%             c.Label.String = '-log10(p)';
        end
        set(gca, 'YDir','normal')
        colormap(cmap)
    end
    
    if g.istitle
        title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    end

    % saving
    save_D = [DIROUT 'D' g.bispec_type '_source_' name '_' int2str(isub) g.f_ext];
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_D)
    else
        exportgraphics(gcf, save_D)
    end
end

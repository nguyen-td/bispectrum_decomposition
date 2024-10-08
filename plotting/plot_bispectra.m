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
%   dim_chan      - [integer] channel over which the bispectrum will becollapsed. Default is 2.

function plot_bispectra(D, f1, f2, isub, name, DIROUT, cmap, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }          '';
        'istitle'        'boolean'    { }        true;
        'f_ext'          'string'     { }      '.png';
        'dim_chan'       'integer'    {1 2 3}       2; 
        });
    if ischar(g), error(g); end
    
    n = size(D, 1);
    D_abs = abs(D);

    if n == 3
        figure('Position', [600 100 1000 300]);
    else
        figure('Position', [600 100 1500 300]);
    end

    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;

        switch g.dim_chan
            case 1
                imagesc(squeeze(D_abs(i, :, :)));
            case 2
                imagesc(squeeze(D_abs(:, i, :)));
            otherwise
                imagesc(squeeze(D_abs(:, :, i)));
        end

        clim([0 max(D_abs, [], 'all')])
        if i == n
            c = colorbar;
%             c.Label.String = '-log10(p)';
        end
        set(gca, 'YDir','normal')
        colormap(cmap)
        xticks(1:n);
        yticks(1:n);

        t = title(i);
        t.FontSize = 20;
    end
    
    if g.istitle
        % title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    end

    % saving
    save_D = [DIROUT 'D' g.bispec_type '_source_' name '_' int2str(isub) g.f_ext];
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_D)
    else
        exportgraphics(gcf, save_D)
    end
end

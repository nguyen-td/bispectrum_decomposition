% Plot and save p-values computed from (source) cross-bispectra.
%
% Inputs:
%   f1, f2        - frequencies in Hz, used for title
%   isub          - subject ID, used for title
%   DIROUT        - output directory to save images
%   cmap          - (n, 3) colormap 
%   P_source_fdr  - (n x n x n) matrix of FDR-corrected p-values of the source cross-bispectrum
%   P_source      - (x x n x n) matrix of p-values of the source cross-bispectrum
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is '' (empty string)
%   istitle       - [boolean] whether to show a title, default is true
%   f_ext         - [string] file extension, default is .png.

function plot_pvalues_bispec_source(f1, f2, isub, DIROUT, cmap, P_source_fdr, P_source, varargin)
    
    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '';
        'istitle'        'boolean'    { }     true;
        'f_ext'          'string'     { }     '.png';
        });
    if ischar(g), error(g); end
    n = size(P_source_fdr, 1);

    % p-values of source cross-bispectrum
    if size(P_source_fdr) == 3
        figure('Position', [600 100 1000 300]);
    else
        figure('Position', [600 100 1500 300]);
    end
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(-log10(squeeze(P_source_fdr(i, :, :))));
%         hold on;
%         imagesc(-log10(squeeze(P_source(i, :, :))), 'AlphaData', 0.7);
        try
            clim([min(-log10(P_source_fdr), [], 'all') max(-log10(P_source_fdr), [], 'all')])
        end
        if i == n
            c = colorbar;
            c.Label.String = '-log10(p)';
        end
        colormap(cmap)
        set(gca, 'YDir','normal')
        xticks(1:n);
        yticks(1:n);
    end
    if g.istitle
        title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    end
    set(gca, 'YDir','normal', 'FontSize', 15)

    % saving figure
    save_P_source = [DIROUT 'P' g.bispec_type '_source_' int2str(isub) g.f_ext];
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_P_source)
    else
        exportgraphics(gcf, save_P_source)
    end
end

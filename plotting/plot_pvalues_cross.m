% Plot and save p-values computed from (source) cross-bispectra.
%
% Inputs:
%   f1, f2        - frequencies in Hz, used for title
%   isub          - subject ID, used for title
%   DIROUT        - output directory to save images
%   P_source_fdr  - (n x n x n) matrix of FDR-corrected p-values of the source cross-bispectrum
%   P_source      - (x x n x n) matrix of p-values of the source cross-bispectrum
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is '_cross' (empty string)


function plot_pvalues_cross(f1, f2, isub, DIROUT, P_source_fdr, P_source, varargin)
    
    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '';
        });
    if ischar(g), error(g); end
    n = size(P_source_fdr, 1);

    % p-values of source cross-bispectrum
    figure('Position', [600 100 1500 300]);
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(-log10(squeeze(P_source_fdr(i, :, :))));
%         hold on;
%         imagesc(-log10(squeeze(P_source(i, :, :))), 'AlphaData', 0.7);
        try
            clim([min(-log10(P_source_fdr), [], 'all') max(-log10(P_source_fdr), [], 'all')])
        end
        if i == 5
            c = colorbar;
            c.Label.String = '-log10(p)';
        end
        set(gca, 'YDir','normal')
    end
    title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    save_P_source = [DIROUT 'P' g.bispec_type '_source_' int2str(isub) '.png'];
    exportgraphics(gcf, save_P_source)
end

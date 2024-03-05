% Plot and save absolute values of the cross-bispectra.
%
% Inputs:
%   D        - (n x n x n) source cross-bispectrum
%   f1, f2   - frequencies in Hz, used for title
%   isub     - subject ID, used for title
%   name     - [string] name for plot, for example, "mixed", "demixed" etc.
%   DIROUT   - output directory to save images

function plot_bispectra(D, f1, f2, isub, name, DIROUT)
    
    D_abs = abs(D);
    n = size(D, 1);
    figure('Position', [600 100 1500 300]);
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(squeeze(D_abs(i, :, :)));
        clim([min(D_abs, [], 'all') max(D_abs, [], 'all')])
        if i == 5
            c = colorbar;
%             c.Label.String = '-log10(p)';
        end
        set(gca, 'YDir','normal')
    end
    title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    save_D = [DIROUT 'D_' name '_' int2str(isub) '.png'];
    exportgraphics(gcf, save_D)
end

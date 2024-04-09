% Create a topomap of a matrix given a selected seed. The first dimension
% of the matrix is then collapsed accordingly. 
%
% Inputs:
%   B         - (n_chan x n_chan) net bispectral matrix
%   seed_idx  - seed index
%   chanlocs  - EEG channel locations from the EEG struct
%   max_val   - maximum value for the topoplot
%   cmap      - (n, 3) colormap 
%   name      - [string] name for plot, for example, "mixed", "demixed" etc.
%   title_str - [string] full title
%   DIROUT    - output directory to save images

function plot_topomaps_seed(B, chanlocs, max_val, cmap, name, title_str, DIROUT)

    % get seed
    [~, max_idx] = max(B, [], 'all'); % get max. value in the row (1st dimension)
    [row, col] = ind2sub(size(B), max_idx);


    % topoplots of mixing matrix
    figure;
    t = title(title_str);
    t.FontSize = 20;
    topoplot(squeeze(B(row, :)), chanlocs, 'electrodes', 'on', 'emarker2', {row, 'o', 'black'}, 'colormap', cmap); colorbar; clim([0 max_val])
    save_B = [DIROUT 'B' name '_seed.png'];
    exportgraphics(gcf, save_B)
    
end
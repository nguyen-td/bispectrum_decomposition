% Create a topomap of a matrix given a selected seed. The first dimension
% of the matrix is then collapsed accordingly. 
%
% Inputs:
%   B         - (n_chan x n_chan) net bispectral matrix
%   seed_idx  - seed index
%   chanlocs  - EEG channel locations from the EEG struct
%   isub      - subject ID, used for title
%   name      - [string] name for plot, for example, "mixed", "demixed" etc.
%   title_str - [string] full title
%   DIROUT    - output directory to save images

function plot_topomaps_seed(B, seed_idx, chanlocs, name, title_str, DIROUT)

    % topoplots of mixing matrix
    figure;
    t = title(title_str);
    t.FontSize = 20;
    topoplot(squeeze(B(seed_idx, :)), chanlocs, 'electrodes', 'on', 'emarker2', {seed_idx, 'o', 'black'}); colorbar; clim([min(B, [], 'all') max(B, [], 'all')])
    save_B = [DIROUT 'B' name '_seed.png'];
    exportgraphics(gcf, save_B)
    
end
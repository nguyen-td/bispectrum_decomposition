% Plot and save mixed and demixed sources as brain plots.
%
% Inputs:
%   F_moca      - (n_voxels x n_dum x n) demixed sources
%   F           - (n_voxels x n_dum x n) mixed sources
%   n           - number of estimated sources
%   cortex75k   - 75k cortex structure from the NYhead for later plotting
%   cortex2k    - 2k subset of the cortex structure
%   source_inds - (n x 1) array of source indicess
%   cm          - colormap
%   isub        - subject ID, used for cortex plot
%   DIROUT      - output directory to save images
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is '' (empty string)

function plot_sources(F_moca, F, n, cortex75k, cortex2k, source_inds, cm, isub, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '';
        });
    if ischar(g), error(g); end

    % get max. absolute value
    max_val = max(sum(F_moca.^2, 2), [], 'all');

    for i = 1:n
        % plot and save demixed source
        source_moca = sum(F_moca(:, :, i).^2, 2); 
        f_name_moca = [DIROUT 'F' g.bispec_type '_demixed' int2str(i) '_' int2str(isub) '_'];
        allplots_cortex_nyhead_v2(cortex75k, source_moca(cortex2k.in_to_cortex75K_geod), [0 max_val], cm, 'demixed sources', 1, f_name_moca, ...
            {cortex75k.vc_smooth(cortex2k.in_from_cortex75K(source_inds), :)})
        
%         % plot and save mixed source
%         source_mixed = sum(F(:, :, i).^2, 2); 
%         f_name_mixed = [DIROUT 'F' g.bispec_type '_mixed' int2str(i) '_' int2str(isub) '_'];
%         allplots_cortex_nyhead_v2(cortex75k, source_mixed(cortex2k.in_to_cortex75K_geod), [0 max_val], cm, 'mixed sources', 1, f_name_mixed, ...
%             {cortex75k.vc_smooth(cortex2k.in_from_cortex75K(source_inds), :)})
    end

end
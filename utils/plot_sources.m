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

function plot_sources(F_moca, F, n, cortex75k, cortex2k, source_inds, cm, isub, DIROUT)

    for i = 1:n
        % plot and save demixed source
        source_moca = sum(F_moca(:, :, i).^2, 2); 
        f_name_moca = [DIROUT '/F_demixed' int2str(i) '_' int2str(isub) '_'];
        allplots_cortex_nyhead(cortex75k, source_moca(cortex2k.in_to_cortex75K_geod), [min(source_moca, [], 'all') max(source_moca, [], 'all')], cm, 'demixed sources', 1, f_name_moca, ...
            {cortex75k.vc_smooth(cortex2k.in_from_cortex75K(source_inds), :)})
        
        % plot and save mixed source
        source_mixed = sum(F(:, :, i).^2, 2); 
        f_name_mixed = [DIROUT '/F_mixed' int2str(i) '_' int2str(isub) '_'];
        allplots_cortex_nyhead(cortex75k, source_mixed(cortex2k.in_to_cortex75K_geod), [min(source_mixed, [], 'all') max(source_mixed, [], 'all')], cm, 'mixed sources', 1, f_name_mixed, ...
            {cortex75k.vc_smooth(cortex2k.in_from_cortex75K(source_inds), :)})
    end

end
% Function to compute true topomaps from the leadfield. 
%
% Inputs:
%   L        - (n_chan x n_voxels x n_dum) leadfield matrix
%   D        - atlas structure, in this case of the Desikan- Killiany atlas, see fp_get_Desikan.m for the full documentation
%   roi_inds - array of ROI indices, e.g., [27, 59]
%   DIROUT   - output directory to save images
%
% Outputs:
%   A_true   - (n_chan x n_roi_inds) true sensor topographies

function A_true = compute_topomap(L, D, roi_inds, DIROUT)
    B = L(:, D.sub_ind_cortex(roi_inds), :);

    % multiply with normal direction to get from three to one dipole dimension 
    normals = D.normals(D.sub_ind_cortex(roi_inds),:)'; 
    for b = 1:numel(D.sub_ind_cortex(roi_inds))
        A_true(:, b) = squeeze(B(:, b, :)) * squeeze(normals(:, b));
    end

    % plot topomap
    load cm17
    chanlocs = readlocs('channel_BrainProducts_ActiCap_97.mat');
    plot_topomaps_patterns(A_true, numel(roi_inds), chanlocs, cm17, '', 'true', DIROUT, 'f_ext', '.fig') 
end
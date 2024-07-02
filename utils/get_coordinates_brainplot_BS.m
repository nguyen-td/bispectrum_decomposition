% Extracts 3D voxels (Vertices) for plotting seed sources in the Brainstorm-based brain plots.
%
% Inputs:
%   cortex  - [struct] Brainstorm-based cortex structure
%   roi_vox - [cell] (1 x n_roi) array of voxels corresponding to the ROIs
%
% Output:
%   voxels  - [array] (n_roi x 3) voxels/coordinates

function voxels = get_coordinates_brainplot_BS(cortex, roi_vox)
    for iroi = 1:numel(roi_vox)
        voxels(iroi, :) = cortex.Vertices(roi_vox{iroi},:);
    end
end
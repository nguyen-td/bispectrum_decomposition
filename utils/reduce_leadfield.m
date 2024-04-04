% Reduce the leadfield to match the EEG channels of the data.
% Please make sure to download the leadfield beforehand: https://www.parralab.org/nyhead/
%
% Input:
%   EEG  - EEG struct
%
% Output:
%   L_3D            - (n_chans x n_voxels x dip_dir) leadfield tensor, dipole directions are typically 3
%   cortex75k       - cortex structure from the NYhead for later plotting
%   cortex2k        - 2k subset of the cortex structure

function [L_3D, cortex75k, cortex2k] = reduce_leadfield(EEG)

    % load the leadfield
    try
        load sa_nyhead
    catch
        warning("Please download the leadfield first: https://www.parralab.org/nyhead/")
    end
    
    % get data channels
    data_struct2cell = struct2cell(EEG.chanlocs);
    data_chans = data_struct2cell(1, 1, :);
    
    % compute intersection
    [intersect_chans, ~, idx_chans] = intersect(data_chans, sa.clab_electrodes);
    if ~(length(intersect_chans) == length(data_chans))
        warning('Not enough channels, use a larger leadfield matrix.')
    end
    
    % user lower-resolution cortex
    L_3D = sa.cortex75K.V_fem(idx_chans, sa.cortex2K.in_from_cortex75K, :);

    % save cortex structures for later plotting
    cortex75k = sa.cortex75K;
    cortex2k = sa.cortex2K;
end
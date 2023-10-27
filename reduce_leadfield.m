% Reduce the leadfield to match the EEG channels of the data.
%
% Input:
%   EEG  - EEG struct
%
% Output:
%   L_3D - (n_chans x n_voxels x dip_dir) leadfield tensor, dipole directions are typically 3

function L_3D = reduce_leadfield(EEG)

    % load the leadfield
    load sa_nyhead
    
    % get data channels
    data_struct2cell = struct2cell(EEG.chanlocs);
    data_chans = data_struct2cell(1, 1, :);
    
    % compute intersection
    [intersect_chans, ~, idx_chans] = intersect(data_chans, sa.clab_electrodes);
    if ~(length(intersect_chans) == length(data_chans))
        warning('Not enough channels, use a larger leadfield matrix.')
    end
    
    % 3D leadfield for the selected subset of channels
    L_3D = sa.cortex75K.V_fem(idx_chans, :, :);
end
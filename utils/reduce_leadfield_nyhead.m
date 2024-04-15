% Reduce the leadfield to match the EEG channels of the data.
% Please make sure to download the leadfield beforehand: https://www.parralab.org/nyhead/
%
% Input:
%   EEG  - EEG struct
%
% Optional input:
%   resolution - ['2k', '5k', '10k', '75k'] resolution of the cortex (it's a lower-case "k"), default is 2k
%
% Output:
%   L_3D            - (n_chans x n_voxels x dip_dir) leadfield tensor, dipole directions are typically 3
%   cortex75k       - cortex structure from the NYhead for later plotting
%   cortex2k        - 2k subset of the cortex structure

function [L_3D, cortex75k, cortex_reduced] = reduce_leadfield_nyhead(EEG, varargin)

    % setup
    g = finputcheck(varargin, { ...
        'resolution'   'string'     { '1k' '2k' '5k' '10k' '75k'}    '2k';
        });
    if ischar(g), error(g); end

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
    switch g.resolution
        case '1k'
            L_3D = sa.cortex75K.V_fem(idx_chans, sa.cortex1K.in_from_cortex75K, :);
            cortex_reduced = sa.cortex1K;
        case '5k'
            L_3D = sa.cortex75K.V_fem(idx_chans, sa.cortex5K.in_from_cortex75K, :);
            cortex_reduced = sa.cortex5K;
        case '10k'
            L_3D = sa.cortex75K.V_fem(idx_chans, sa.cortex10K.in_from_cortex75K, :);
            cortex_reduced = sa.cortex10K;
        case '75k'
            L_3D = sa.cortex75K.V_fem(idx_chans, :, :);
            cortex_reduced = sa.cortex75K; % here, cortex75K == cortex_reduced
        otherwise % '2k'
            L_3D = sa.cortex75K.V_fem(idx_chans, sa.cortex2K.in_from_cortex75K, :);
            cortex_reduced = sa.cortex2K;
    end


    % save cortex structures for later plotting
    cortex75k = sa.cortex75K;
end
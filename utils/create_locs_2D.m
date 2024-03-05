% Create matrix containing the 2D positions of sensors in a plane for later isualization.
% Please make sure to download the leadfield beforehand: https://www.parralab.org/nyhead/
%
% Input:
%   EEG  - EEG struct
%
% Output:
%   locs_2D - positions of sensors in a plan (for visualization of topographies using Guido's routines)

function locs_2D = create_locs_2D(EEG)

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
    
    % locs_2D structure for the selected subset of channels
    locs_2D = sa.locs_2D(idx_chans,:);
    locs_2D(:, 1) = 1:length(idx_chans); 
end
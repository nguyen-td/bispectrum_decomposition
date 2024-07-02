% Compute the Euclidian distance between the true and estimated sources (point sources). 
% For the estimated source, the maximum over the source distribution vector (which is used to create brain plots) 
% is computed and the corresponding sources (in three dimensions) are estimated.  
%
% Inputs:
%   F_true - (n_voxels x n_dum x n_sources) array of true source activations
%   F_est  - (n_voxels x n_dum x n_sources) array of estimated source activations
%
% Outputs:
%   eucl_dist  - Euclidian distance (MNI space) in mm

function eucl_dist = calc_source_eucl_dist(F_true, F_est)
    
    assert(all(size(F_true) == size(F_est)), 'Sizes of F_true and F_est must be equal.')
    n = size(F_true, 3);
    for i = 1:n
        % compute maximum of source magnitues (peaks on brain plots)
        [~, max_true] = max(sum(F_true(:, :, i).^2, 2)); 
        [~, max_est] = max(sum(F_est(:, :, i).^2, 2)); 

        % compute subspace angle
        eucl_dist = norm(F_true(max_true, :, i) - F_est(max_est, :, i));
    end
end
% Permute sources after using the permutation matrix from the Hungarian
% algorithm.
%
% Inputs:
%   F - (n_voxels x n_dum x n_sources) sources
%   P - (n_sources x n_sources) permutation matrix
%
% Output:
%   F_matched - (n_voxels x n_dum x n_source) matched sources

function F_matched = match_sources(F, P)
    
    [n_voxels, n_dum, n_sources] = size(F);
    resh = reshape(F, [], n_sources);
    F_matched = reshape((P * resh')', n_voxels, n_dum, n_sources);
end


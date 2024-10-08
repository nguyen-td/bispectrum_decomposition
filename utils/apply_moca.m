% Unmix the sources and get a new mixing matrix.
%
% Inputs:
%   L_3D  - (n_chans x n_voxels x n_dum) leadfield tensor, dipole directions are typically 3 
%   A     - (n_chan x n) mixing matrix, where n is the model order (number of estimated sources)
%   n     - number of estimated sources
%   regu  - regularization term for eLORETA, default is 0.05.
%
% Output:
%   A_moca - (n x n) demixing matrix, where n is the model order (number of estimated sources)
%   F_moca - (n_voxels x n_dum x n) demixed sources
%   F      - (n_voxels x n_dum x n) mixed sources

function [A_moca, F_moca, F] = apply_moca(L_3D, A, n, regu)

    % calculate 3D eLORETA inverse filter 
    if nargin < 4; regu = 0.05; end
    disp('Computing eLORETA...');
    P_eloreta = mkfilt_eloreta_v2(L_3D, regu); 

    % find mixed sources
    [nchan, nvoxel, ndum] = size(P_eloreta);
    
    % intitialze distributions for n sources 
    F = zeros(nvoxel, ndum, n); 

    for i = 1:n 
        for k=1:ndum
            F(:, k, i) = P_eloreta(:,:,k)' * A(:,i); 
        end
    end

    % now unmix the sources 
    if n > 1
        [F_moca, A_moca] = moca_ncomp(F);
    else
        warning('Since there is only a single source, it will not be demixed.')
        F_moca = F;
        A_moca = ones(n, n);
    end
end
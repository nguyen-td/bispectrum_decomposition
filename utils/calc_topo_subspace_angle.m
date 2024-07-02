% Calculate the subspace angle between two topomaps A and B. 
% 
% If one matrix has more columns than the other, the one with more columns will be
% truncated. For example, if A is a 68x3 matrix and B is a 68x2 matrix, the
% third column of A will be removed. A and B must have the same same number
% of rows.

function theta = calc_topo_subspace_angle(A, B)

    if ~isequal(size(A, 1), size(B, 1))
        error 'A and B must have the same number of rows.'
    end

    % ensure that both matrices have the same number of colums
    n_cols = min(size(A, 2), size(B, 2));
    A = A(:, 1:n_cols); 
    B = B(:, 1:n_cols);

    % compute subspace angle
    for i = 1:size(A, 2)
        theta(i) = subspace(A(:, i), B(:, i));
    end
end
% Function to set diagonals to a selected value.
%
% Input:
%   mat        - (n x n x n) real array
%   diag_value - [integer] value the diagonals should take
%
% Output:
%   mat_new - (n x n x n) real array where its diagonals are zero

function mat = set_diagonals(mat, diag_value)

    n = size(mat, 1);
    for i = 1:n
        mat(i, i, i) = diag_value;
    end
end
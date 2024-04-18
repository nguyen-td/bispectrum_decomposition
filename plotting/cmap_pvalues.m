% Determine which colormap to use for p-value plotting. If all entries are 1, use cm17, else use
% cm17a.
%
% Input:
%   P           - matrix of p-values
%   cm17, cm17a - colormaps

function cmap = cmap_pvalues(P, cm17, cm17a)
    if all(P(:) == P(1))
        cmap = cm17;
    else
        cmap = cm17a;
    end
end
% Extract frequency indices (bins) from inputs in Hz.
%
% Inputs:
%   f1, f2 - frequency in Hz, can either be single integers or an array, i.e., either f1 = 11 or f1 = [9 11];
%   frqs   - (n_freqs x 1) vector of frequencies, spaced according to the frequency resolution 
%
% Outputs:
%   freqpairs - (n_freqcombinations x 2) matrix containing all frequency combinations in bins

function freqpairs = get_freqindices(f1, f2, frqs)   

    % extract all individual frequencies in the selected bands
    size_low = size(f1, 2);
    size_high = size(f2, 2);
    mask_inds_low = frqs >= f1(1) & frqs <= f1(size_low);
    mask_inds_high = frqs >= f2(1) & frqs <= f2(size_high);
    frqs_low = frqs(mask_inds_low); 
    frqs_high = frqs(mask_inds_high);
    
    % determine all frequency combinations
    [k, j] = ndgrid(frqs_low, frqs_high);
    frqs_combs = [k(:), j(:)]; 
    n_combs = size(frqs_combs, 1);
    freqpairs = zeros(n_combs, 2);
    warning('The surrogate bispectra are going to be estimated on %d frequency pair(s).', n_combs);
    for i = 1:n_combs
        low = frqs_combs(i, 1);
        high = frqs_combs(i, 2);
        freqpairs(i, :) = [find(frqs == low), find(frqs == high)];
    end
end
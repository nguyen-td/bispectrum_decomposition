% Compute p-values from true value and shuffles. 
%
% Inputs:
%   true_val  - array containing the true value
%   shuf_vals - array containing the shuffled values, has to have the same dimension as true_val with 
%                an additional shuffle dimension, e.g., if true_val is (n_chan x n_chan x n_chan), then 
%                shuff_vals is (n_chan x n_chan x n_chan, n_shuf)
%   nshuf     - number of shuffles
%   alpha     - significance level, default is 0.05.
%
% Outputs:
%   P    - array containing p-values, has the same dimension as true_val
%   P_fdr - array containing FDR-corrected values, has the same dimension as P

function [P, P_fdr] = compute_pvalues(true_val, shuf_vals, n_shuf, alpha)

    % compute p-values, take mean over regions
    P = squeeze(sum(true_val < shuf_vals, 4) ./ n_shuf);
    P(P==0) = 1 / n_shuf;
    
    % correct for multiple comparisons
    [p_fdr, ~] = fdr(P, alpha);
    P_fdr = P;
    P_fdr(P > p_fdr) = 1;
end
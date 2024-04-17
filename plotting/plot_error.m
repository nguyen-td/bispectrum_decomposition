% Function to plot the error over iterations of the bsfit_freqbands.m
% function.
%
% Inputs:
%   errors       - [cell array] (n_freqcombs x n) errors over iterations
%   freqcomb_idx - [integer] index determining for which frequency combination the error will be plotted: errors{freqcomb_idx, :}
%   n            - [integer] model order
%   colors       - [array] list of line colors, e.g. ['r', 'b', 'k', 'g']
%   title_str    - [string] title
%   islog        - [boolean] whether to plot on log scale or not, default is 0

function plot_error(errors, freqcomb_idx, n, colors, title_str, islog)
    
    if nargin < 6
        islog = false;
    end

    max_err = max(cellfun(@max, errors), [], 'all');
    min_err = min(cellfun(@min, errors), [], 'all');
    figure; 
    hold on;
    for ind = 1:length(n)
        if islog
            plot(log(errors{freqcomb_idx, ind}), colors(ind), 'DisplayName', int2str(n(ind)))
        else
            plot(errors{freqcomb_idx, ind}, colors(ind), 'DisplayName', int2str(n(ind)))
        end
    end
    legend
    hold off;
    grid on;
    xlabel('Iteration')
    if islog
        ylabel('log(Error)')
        ylim([log(min_err)-1 log(max_err)+1])
    else
        ylabel('Error')
        ylim([0 max_err])
    end
    legend('Location', 'SouthEast')
    set(gca, 'FontSize', 15)
    t = title(title_str);
    t.FontSize = 15;
end
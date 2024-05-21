% Function to plot the error over iterations of the bsfit_freqbands.m
% function.
%
% Inputs:
%   errors       - [cell array] (n_freqcombs x n) errors over iterations 
%   freqcomb_idx - [integer] index determining for which frequency combination the error will be plotted: errors{freqcomb_idx, :}
%   n            - [integer] model order
%   colors       - [array] list of line colors, e.g. ['r', 'b', 'k', 'g']
%   title_str    - [string] title
%   isub         - subject ID, used for file name
%   DIROUT       - output directory to save images
%
% Optional input:
%   f_name  - [string] part of file name after 'err', default is '' (empty string)
%   islog   - [boolean] whether to plot on log scale or not, default is false
%   f_ext   - [string] file extension, default is .png.

function plot_error(errors, freqcomb_idx, n, colors, title_str, isub, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'f_name'    'string'     { }     '';
        'islog'     'boolean'    { }     false;
        'f_ext'    'string'     { }     '.png';
         });
    if ischar(g), error(g); end

    max_err = max(cellfun(@max, errors), [], 'all');
    min_err = min(cellfun(@min, errors), [], 'all');
    figure; 
    hold on;
    for ind = 1:length(n)
        if ~colors == 0
            if g.islog
                plot(log(errors{freqcomb_idx, ind}), colors(ind), 'DisplayName', int2str(n(ind)))
            else
                plot(errors{freqcomb_idx, ind}, colors(ind), 'DisplayName', int2str(n(ind)))
            end
        else
            if g.islog
                plot(log(errors{freqcomb_idx, ind}), 'DisplayName', int2str(n(ind)))
            else
                plot(errors{freqcomb_idx, ind}, 'DisplayName', int2str(n(ind)))
            end
        end
    end
    legend
    hold off;
    grid on;
    xlabel('Iteration')
    if g.islog
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
    
    % save figure
    save_err = [DIROUT 'err' g.f_name '_' int2str(isub) g.f_ext];
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_err)
    else
        exportgraphics(gcf, save_err)
    end
end
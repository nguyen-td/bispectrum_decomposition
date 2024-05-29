% Plot and save 2D bispectral indices according to an index scheme. For
% example, if the *B_ijj* slice should be plotted, the function will create
% and plot the following matrix:
%
% | 1x1x1 | 2x1x1 | 3x1x1 |
% | 1x2x2 | 2x2x2 | 3x2x2 |
% | 1x3x3 | 2x3x3 | 3x3x3 |
%
% If the *B_jij* slice should be plotted, the function will create the
% following matrix:
%
% | 1x1x1 | 2x1x2 | 3x1x3 |
% | 1x2x1 | 2x2x2 | 3x2x3 |
% | 1x3x1 | 2x3x2 | 3x3x3 |
%
% Note that i is the inner index (row index) and j is is the outer index (column index). 
% B(ichan, jchan) fills the lower triangle, B(jchan, ichan) the upper triangle  of the B_sliced 
% matrix.
%
% Inputs:
%   B          - (n_chan x n_chan x n_chan) cross-bispectrum, must be real-valued
%   slice_idx  - [idx idx idx] array that indicates the bispectral slice (read docstring above), e.g., [1 2 2] 
%                would extract and plot B_ijj slices.
%   isub       - subject ID, used for file name
%
% Optional inputs:
%   cbar_label  - [string] label of the colorbar, default is '|Bispectrum|'
%   isplot      - [boolean] whether to plot or not, default is true
%   f_name      - [string] part of file name after 'P_slice', default is '' (empty string)
%   f_ext       - [string] file extension, default is .png.
%   isplot_pval - [boolean] whether to plot p-values or bispectral slices
%
% Output:
%   B_sliced   - (n_chan x n_chan) bispectral slice

function B_sliced = plot_bispec_slices(B, slice_idx, cmap, isub, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'cbar_label'    'string'     { }     '|Bispectrum|';
        'isplot'        'boolean'    { }     true;    
        'f_name'        'string'     { }     '';
        'f_ext'         'string'     { }     '.png';
        'isplot_pval'   'boolean'    { }     true;
        });
    if ischar(g), error(g); end
    
    n_chan = size(B, 1);
    B_sliced = zeros(n_chan, n_chan);
    
    % create bispectral slice
    for jchan = 1:n_chan
        for ichan = jchan:n_chan
            if isequal(slice_idx, [1 2 2]) % B_ijj
                B_sliced(ichan, jchan) = B(jchan, ichan, ichan);
                B_sliced(jchan, ichan) = B(ichan, jchan, jchan);
            elseif isequal(slice_idx, [2 1 2]) % B_jij
                B_sliced(ichan, jchan) = B(ichan, jchan, ichan);
                B_sliced(jchan, ichan) = B(jchan, ichan, jchan);
            elseif isequal(slice_idx, [2 2 1]) % B_jji
                B_sliced(ichan, jchan) = B(ichan, ichan, jchan);
                B_sliced(jchan, ichan) = B(jchan, jchan, ichan);
            elseif isequal(slice_idx, [1 1 1]) % B_iii
                if jchan == ichan
                    B_sliced(ichan, jchan) = B(jchan, jchan, jchan);
                    B_sliced(jchan, ichan) = B(ichan, ichan, ichan);
                end
            else
                warning('Only B_iii, B_ijj, B_jij and B_jji slices can be extracted. If you want to extract any other slice, you will have to define it yourself first.')
            end
        end
    end
    
    % plotting
    if g.isplot
        figure; imagesc(B_sliced)
        colormap(cmap)
        c = colorbar();
        c.Label.String = g.cbar_label;
        set(gca, 'YDir','normal', 'FontSize', 15)
        xticks(1:n_chan);
        yticks(1:n_chan);
    end

    % save figure
    if g.isplot_pval
        save_B = [DIROUT 'P_slice' g.f_name '_' int2str(isub) g.f_ext]; 
    else
        save_B = [DIROUT 'B_slice' g.f_name '_' int2str(isub) g.f_ext]; 
    end
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_B)
    else
        exportgraphics(gcf, save_B)
    end
end
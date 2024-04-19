% Plot and save absolute values of the univariate (sensor) bispectrum. The mean is being computed over the first channel.
% Can also be used for cross-bispectra as long as the dimensions coincide. Check optional inputs for adjusting axis labels.
%
% Inputs:
%   B        - (n_freq x n_freq) univariate bispectrum, n_freq can be any unit, e.g., n_chan would also work
%   frqs     - (n_frqs x 1) array of frequencies, used for labels
%   isub     - subject ID, used for title
%   cmap     - (n, 3) colormap 
%   DIROUT   - output directory to save images
%
% Optional inputs:
%   bispec_type   - [string] type of bispectrum (for file name), default is "_univ"
%   label_x       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   label_y       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   custom_label  - [boolean] whether to show custom labels, only activate when plotting univariate bispectra where 
%                   both axes correspond to frequencies
%   title_str     - [string] title, default is "Univariate sensor bispectrum)"
%   mean_chan     - [1, 2, 3] channel along which the mean will be computed (for the title), default is 1 (channel 1)

function plot_bispectra_univ(B, frqs, isub, cmap, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '_univ';
        'label_x'        'string'     { }     '$f_1$ (Hz)';
        'label_y'        'string'     { }     '$f_2$ (Hz)';
        'custom_label'   'boolean'    { }     1;
        'title_str'      'string'     { }     'p-values (univariate sensor bispectrum)';
        'mean_chan'      'integer'    { 1, 2, 3 }  1;
        'isnorm'         'boolean'    { }      false;
        });
    if ischar(g), error(g); end

    figure;
    if g.isnorm
        imagesc(B, [0 1])
    else
        imagesc(B)
    end
    colormap(cmap)
    c = colorbar();
    c.Label.String = '|Bispectrum|';
    title(g.title_str)
    set(gca, 'YDir','normal')
    
    if g.custom_label
        label_idx = 1:20:size(B, 2); % show only every 20th label
        custom_label = 1:size(B, 2);
        set(gca, 'XTick', custom_label(label_idx));
        set(gca, 'XTickLabel', frqs(label_idx));
        set(gca, 'YTick', custom_label(label_idx));
        set(gca, 'YTickLabel', frqs(label_idx));
    end
    xlabel(g.label_x, 'Interpreter', 'latex');
    ylabel(g.label_y, 'Interpreter', 'latex');
    save_B_sensor = [DIROUT 'B' g.bispec_type '_chan' int2str(g.mean_chan) '_sensor_' int2str(isub) '.png'];
    exportgraphics(gcf, save_B_sensor)

end

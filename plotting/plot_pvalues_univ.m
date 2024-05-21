% Plot and save p-values computed from univariate surrogate (sensor) bispectra. Can also 
% be used for cross-bispectra as long as the dimensions coincide. Check optional inputs for adjusting axis labels.

%
% Inputs:
%   P_sens_fdr    - (n_freq x n_freq) matrix of FDR-corrected p-values of the sensor bispectrum
%   frqs          - (n_frqs x 1) array of frequencies, used for labels
%   isub          - subject ID, used for title
%   cmap          - (n, 3) colormap 
%   DIROUT        - output directory to save images
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is "_univ"
%   label_x       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   label_y       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   custom_label  - [boolean] whether to show custom labels, only activate when plotting univariate bispectra where 
%                   both axes correspond to frequencies
%   title_str     - [string] title, default is "p-values (univariate sensor bispectrum)"
%   f_ext         - [string] file extension, default is .png.
%   label_latex   - [boolean] whether x and y labels should be printed using the Latex interpreter

function plot_pvalues_univ(P_sens_fdr, frqs, isub, cmap, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '_univ';
        'label_x'        'string'     { }     '$f_1$ (Hz)';
        'label_y'        'string'     { }     '$f_2$ (Hz)';
        'custom_label'   'boolean'    { }     1;
        'title_str'      'string'     { }     'p-values (univariate sensor bispectrum)';
        'f_ext'          'string'     { }     '.png';
        'label_latex'    'boolean'    { }     true;
        });
    if ischar(g), error(g); end

    figure;
    imagesc(-log10(P_sens_fdr));
    colormap(cmap)
    c = colorbar();
    c.Label.String = '-log10(p)';
    title(g.title_str)
    set(gca, 'YDir','normal')
    
    if g.custom_label
        label_idx = 1:20:size(P_sens_fdr, 2); % show only every 20th label
        custom_label = 1:size(P_sens_fdr, 2);
        set(gca, 'XTick', custom_label(label_idx));
        set(gca, 'XTickLabel', frqs(label_idx));
        set(gca, 'YTick', custom_label(label_idx));
        set(gca, 'YTickLabel', frqs(label_idx));
    end
    if g.label_latex
        xlabel(g.label_x, 'Interpreter', 'latex');
        ylabel(g.label_y, 'Interpreter', 'latex');
    else
        xlabel(g.label_x);
        ylabel(g.label_y);
    end
    set(gca, 'YDir','normal', 'FontSize', 15)
    save_P_sensor = [DIROUT 'P' g.bispec_type '_sensor_' int2str(isub) g.f_ext];
    if strcmpi(g.f_ext, '.fig')
        saveas(gcf, save_P_sensor)
    else
        exportgraphics(gcf, save_P_sensor)
    end

end
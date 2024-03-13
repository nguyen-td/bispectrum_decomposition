% Plot and save p-values computed from univariate surrogate (sensor) bispectra. Can also 
% be used for cross-bispectra as long as the dimensions coincide. Check optional inputs for adjusting axis labels.

%
% Inputs:
%   P_sens_fdr    - (n_freq x n_freq) matrix of FDR-corrected p-values of the sensor bispectrum
%   frqs          - (n_frqs x 1) array of frequencies, used for labels
%   isub          - subject ID, used for title
%   DIROUT        - output directory to save images
%
% Optional input:
%   bispec_type   - [string] type of bispectrum (for file name), default is "_univ"
%   label_x       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   label_y       - [string] label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   custom_label  - [boolean] whether to show custom labels, only activate when plotting univariate bispectra where 
%                   both axes correspond to frequencies
%   title_str     - [string] title, default is "p-values (univariate sensor bispectrum)"

function plot_pvalues_univ(P_sens_fdr, frqs, isub, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     '_univ';
        'label_x'        'string'     { }     '$f_1$ (Hz)';
        'label_y'        'string'     { }     '$f_2$ (Hz)';
        'custom_label'   'boolean'    { }     1;
        'title_str'      'string'     { }     'p-values (univariate sensor bispectrum)';
        });
    if ischar(g), error(g); end

    figure;
    imagesc(-log10(P_sens_fdr)); 
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
    xlabel(g.label_x, 'Interpreter', 'latex');
    ylabel(g.label_y, 'Interpreter', 'latex');
    save_P_sensor = [DIROUT 'P' g.bispec_type '_sensor_' int2str(isub) '.png'];
    exportgraphics(gcf, save_P_sensor)

end
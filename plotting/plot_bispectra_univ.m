% Plot and save absolute values of the univariate (sensor) bispectrum. The mean is being computed over the first channel.
%
% Inputs:
%   frqs     - (n_frqs x 1) array of frequencies, used for labels
%   B        - (n_chan x n_freq x n_freq) univariate bispectrum
%   isub     - subject ID, used for title
%   name     - name, indicating if it is the normalized or unnormalized bispectrum (use upper case)
%   DIROUT   - output directory to save images
%
% Optional inputs:
%   bispec_type - type of bispectrum, default is "univariate sensor bispectrum"
%   label_x     - label x-axis, default is "f1 (Hz)" with LaTeX formatting
%   label_y     - label x-axis, default is "f1 (Hz)" with LaTeX formatting

function plot_bispectra_univ(frqs, B, isub, name, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'    'string'     { }     'univariate sensor bispectrum';
        'label_x'        'string'     { }     '$f_1$ (Hz)';
        'label_y'        'string'     { }     '$f_2$ (Hz)';
        });
    if ischar(g), error(g); end

    figure;
    imagesc(squeeze(mean(abs(B), 1)))
    c = colorbar();
    c.Label.String = '|Bispectrum|';
    title([name ' ' g.bispec_type ', subject ' int2str(isub)])
    set(gca, 'YDir','normal')

    label_idx = 1:20:size(B, 2); % show only every 20th label
    custom_label = 1:size(B, 2);
    set(gca, 'XTick', custom_label(label_idx));
    set(gca, 'XTickLabel', frqs(label_idx));
    set(gca, 'YTick', custom_label(label_idx));
    set(gca, 'YTickLabel', frqs(label_idx));
    xlabel(g.label_x, 'Interpreter', 'latex');
    ylabel(g.label_y, 'Interpreter', 'latex');
    save_B_sensor = [DIROUT 'B_sensor_' lower(name) '_' int2str(isub) '.png'];
    exportgraphics(gcf, save_B_sensor)

end

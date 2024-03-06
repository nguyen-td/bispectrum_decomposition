% Plot and save absolute values of the univariate (sensor) bispectrum.
%
% Inputs:
%   frqs     - (n_frqs x 1) array of frequencies, used for labels
%   B        - (n_chan x n_freq x n_freq) univariate bispectrum
%   isub     - subject ID, used for title
%   name     - name, indicating if it is the normalized or unnormalized bispectrum (use upper case)
%   DIROUT   - output directory to save images

function plot_bispectra_univ(frqs, B, isub, name, DIROUT)

    figure;
    imagesc(squeeze(mean(abs(B), 1)))
    c = colorbar();
    c.Label.String = '|Bispectrum|';
    title([name ' univariate sensor bispectrum, subject ' int2str(isub)])
    set(gca, 'YDir','normal')

    label_idx = 1:20:size(B, 2); % show only every 20th label
    custom_label = 1:size(B, 2);
    set(gca, 'XTick', custom_label(label_idx));
    set(gca, 'XTickLabel', frqs(label_idx));
    set(gca, 'YTick', custom_label(label_idx));
    set(gca, 'YTickLabel', frqs(label_idx));
    xlabel('$f_1$ (Hz)', 'Interpreter', 'latex');
    ylabel('$f_2$ (Hz)', 'Interpreter', 'latex');
    save_B_sensor = [DIROUT 'B_sensor_' lower(name) '_' int2str(isub) '.png'];
    exportgraphics(gcf, save_B_sensor)

end

% Plot and save p-values only for the univariate (sensor) bispectrum.
%
% Inputs:
%   frqs          - (n_frqs x 1) array of frequencies, used for labels
%   isub          - subject ID, used for title
%   DIROUT        - output directory to save images
%   P_sens_fdr    - (n_freq x n_freq) matrix of FDR-corrected p-values of the sensor bispectrum
%   P_sens        - (n_freq x n_freq) matrix of p-values of the sensor bispectrum

function plot_pvalues_univ(frqs, isub, DIROUT, P_sens_fdr)

    figure;
    imagesc(-log10(P_sens_fdr)); 
    c = colorbar();
    c.Label.String = '-log10(p)';
    title(sprintf('p-value (univariate sensor bispectrum), subject %d', isub))
    axis([frqs(2) frqs(end-1) frqs(2) frqs(end-1)]) % set axis limit
    set(gca, 'YDir','normal')
    xlabel('$f_1$ (Hz)', 'Interpreter', 'latex');
    ylabel('$f_2$ (Hz)', 'Interpreter', 'latex');
    save_P_sensor = [DIROUT 'P_sensor_' int2str(isub) '.png'];
    exportgraphics(gcf, save_P_sensor)

end
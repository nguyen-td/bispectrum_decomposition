% Plot and save p-values.
%
% Inputs:
%   A_sens        - (n_chan x n) mixing matrix
%   f1, f2        - frequencies in Hz, used for title
%   frqs          - (n_frqs x 1) array of frequencies, used for labels
%   isub          - subject ID, used for title
%   chanlocs      - EEG channel locations from the EEG struct
%   DIROUT        - output directory to save images
%   P_source_fdr  - (n x n x n) matrix of FDR-corrected p-values of the source cross-bispectrum
%   P_source      - (x x n x n) matrix of p-values of the source cross-bispectrum
%   P_sens_fdr    - (n_freq x n_freq) matrix of FDR-corrected p-values of the sensor bispectrum
%   P_sens        - (n_freq x n_freq) matrix of p-values of the sensor bispectrum

function plot_pvalues(A_sens, f1, f2, frqs, isub, chanlocs, DIROUT, P_source_fdr, P_source, P_sens_fdr, P_sens)
        
    if nargin < 10
        freq_manual = 1;
    else
        freq_manual = 0;
    end
    n = size(P_source_fdr, 1);

    % p-values of source cross-bispectrum
    figure('Position', [600 100 1500 300]);
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(-log10(squeeze(P_source_fdr(i, :, :))));
%         hold on;
%         imagesc(-log10(squeeze(P_source(i, :, :))), 'AlphaData', 0.7);
        try
            clim([min(-log10(P_source_fdr), [], 'all') max(-log10(P_source_fdr), [], 'all')])
        end
        if i == 5
            c = colorbar;
            c.Label.String = '-log10(p)';
        end
        set(gca, 'YDir','normal')
    end
    title(tl1, sprintf('Subject %d, f1 = %d Hz, f2 = %d Hz', isub, f1, f2))
    save_P_source = [DIROUT '/P_source_' int2str(isub) '.png'];
    exportgraphics(gcf, save_P_source)

    % topoplots of mixing matrix
    figure('Position', [600 100 1500 300]);
    tl2 = tiledlayout(1, n);
    for i = 1:n
        nexttile
        t = title(i);
        t.FontSize = 20;
        if i == n
            topoplot(A_sens(:, i), chanlocs, 'electrodes', 'on'); colorbar; clim([min(A_sens, [], 'all') max(A_sens, [], 'all')])
        else
            topoplot(A_sens(:, i), chanlocs, 'electrodes', 'on'); clim([min(A_sens, [], 'all') max(A_sens, [], 'all')])
        end
    end
    save_A = [DIROUT '/A_' int2str(isub) '.png'];
    exportgraphics(gcf, save_A)

    % p-values of sensor bispectrum
    if ~freq_manual
        figure;
        imagesc(-log10(P_sens_fdr)); 
        c = colorbar();
        c.Label.String = '-log10(p)';
        title('Univariate sensor bispectrum')
    %     hold on;
    %     imagesc(-log10(P_sens), 'AlphaData', 0.7)
        axis([frqs(2) frqs(end-1) frqs(2) frqs(end-1)]) % set axis limit
        set(gca, 'YDir','normal')
        xlabel('Frequency (Hz)')
        ylabel('Frequency (Hz)')
        save_P_sensor = [DIROUT '/P_sensor_' int2str(isub) '.png'];
        exportgraphics(gcf, save_P_sensor)
    end
end

% Plot and save p-values from the frequency selection step and for the
% final p-value tensor.
%
% Inputs:
%   P_sens_fdr    - (n_freq x n_freq) matrix of FDR-corrected p-values of the sensor bispectrum
%   P_sens        - (n_freq x n_freq) matrix of p-values of the sensor bispectrum
%   P_source_fdr  - (n x n x n) matrix of FDR-corrected p-values of the source cross-bispectrum
%   P_source      - (x x n x n) matrix of p-values of the source cross-bispectrum
%   A_hat         - (n_chan x n) mixing matrix
%   freq_manual   - manual frequency selection. Default is 'off', i.e., frequencies will be selected automatically.
%   f1, f2        - frequencies in Hz, used for title
%   isub          - subject ID, used for title

function plot_pvalues(P_sens_fdr, P_sens, P_source_fdr, P_source, A_hat, f1, f2, isub)

    n = size(P_source, 1);
    load('MotorImag/data/chanlocs.mat')

    % p-values of source cross-bispectrum
    tl1 = tiledlayout(1, n);
    for i = 1:n 
        nexttile;
        imagesc(-log10(squeeze(P_source_fdr(i, :, :))));
        hold on;
        imagesc(-log10(squeeze(P_source(i, :, :))), 'AlphaData', 0.7);
        clim([min(-log10(P_source), [], 'all') max(-log10(P_source), [], 'all')])
        if i == 5
            c = colorbar;
            c.Label.String = '-log10(p)';
        end
    end

    % topoplots of mixing matrix
    tl2 = tiledlayout(1, n);
    for i = 1:n
        nexttile
        t = title(i);
        t.FontSize = 20;
        if i == n
            topoplot(A_hat(:, i), chanlocs, 'electrodes', 'on'); colorbar; clim([min(A_hat, [], 'all') max(A_hat, [], 'all')])
        else
            topoplot(A_hat(:, i), chanlocs, 'electrodes', 'on'); clim([min(A_hat, [], 'all') max(A_hat, [], 'all')])
        end
    end

    % p-values of sensor bispectrum
    imagesc(-log10(P_sens_fdr))
    hold on;
    imagesc(-log10(P_sens), 'AlphaData', 0.7)
    xlabel('Frequency (Hz)')
    ylabel('Frequency (Hz)')
end

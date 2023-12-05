% Run ICA decomposition on the EEG data and plot the first n components.
%
% Inputs:
%   EEG        - EEGLAB struct
%   n          - model order/number of fitted sources
%   f_chanlocs - directory to load channel locations "chanlocs.mat"

function run_ica(EEG, n, f_chanlocs)

    load(f_chanlocs) 
    OUT_EEG = pop_runica(EEG, 'icatype', 'runica', 'pca', n);
    W = OUT_EEG.icawinv;
    
    % topoplots of mixing matrix
    figure('Position', [600 100 1500 300]);
    tl2 = tiledlayout(1, n);
    for i = 1:n
        nexttile
        t = title(i);
        t.FontSize = 20;
        if i == n
            topoplot(W(:, i), chanlocs, 'electrodes', 'on'); colorbar; clim([min(W, [], 'all') max(W, [], 'all')])
        else
            topoplot(W(:, i), chanlocs, 'electrodes', 'on'); clim([min(W, [], 'all') max(W, [], 'all')])
        end
    end
    save_WICA = [DIROUT '/W_ICA_' int2str(isub) '.png'];
    exportgraphics(gcf, save_WICA)
end
function plt_fulleeg
    % plot raw EEG signal of vp2
    addpath('home/tdnguyen/data')
    load('home/tdnguyen/data/vp2.mat');
    figure1 = figure;
    plot(cnt.x');
    saveas(figure1, 'vp2.fig')
end
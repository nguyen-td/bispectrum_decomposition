function plt_fulleeg
    % addpath
    addpath('/home/tdnguyen/data')

    % plot raw EEG signal of vp2
    load('/home/tdnguyen/data/vp2.mat');
    figure1 = figure;
    plot(cnt.x');
    saveas(figure1, 'vp2.fig')
end
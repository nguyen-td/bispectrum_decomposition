function plt_fulleeg
    figure1 = figure;
    x = 0:1:1000;
    y = x.^2;
    plot(y);
    saveas(figure1, 'vp2.jpg')
end
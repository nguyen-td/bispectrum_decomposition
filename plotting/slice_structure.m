load cm17.mat

A = zeros(3,3);
labels = {'$B_{111}$', '0', '0'; '0', '$B_{222}$', '0'; '0', '0', '$B_{333}$'}; % Matrix of strings
% labels = {'$B_{111}$', '$B_{122}$', '$B_{133}$'; '$B_{211}$', '$B_{222}$', '$B_{233}$'; '$B_{311}$', '$B_{322}$', '$B_{333}$'}; % Matrix of strings

figure;
imagesc(A)
colormap(cm17)
xticks([0 1.5 2.5 3.5]);
yticks([0 1.5 2.5 3.5]);
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
grid on
set(gca, 'YDir','normal', 'FontSize', 15)
title('$B_{iii}(f_1, f_2)$', 'Interpreter', 'latex', 'FontSize', 20)
% title('$B_{ijj}(f_1, f_2)$', 'Interpreter', 'latex', 'FontSize', 20)

[nRows, nCols] = size(labels);
for row = 1:nRows
    for col = 1:nCols
        text(col, row, labels{row, col}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex')
    end
end

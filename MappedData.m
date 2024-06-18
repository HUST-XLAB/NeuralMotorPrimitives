load('figuredata.mat');

figure();
imagesc(M5MappedData, [0, 1])
title('DSM 5 fingers');
set(gca, 'YTick', 1:5, 'YTickLabel', labels_5);
set(gca, 'XTick', 1:5, 'XTickLabel', labels_5);
%colormap(aa);
colorbar

figure();
imagesc(M4MappedData, [0, 1])
title('DSM 4 motion primitives');
set(gca, 'YTick', 1:5, 'YTickLabel', labels_5);
set(gca, 'XTick', 1:5, 'XTickLabel', labels_5);
%colormap(aa);
colorbar

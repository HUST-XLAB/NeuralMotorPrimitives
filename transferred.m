% Fig5 c
clc; clear;
c1 = [240, 170, 63] / 255;
c2 = [178, 178, 178] / 255;
num = [2 11; 13 6; 19 1];
figure('Color', 'w')
b = bar(num, 0.9);
set(b, 'FaceAlpha', 0.9);
set(b, 'edgecolor', 'none');
b(1).FaceColor = c1;
b(2).FaceColor = c2;
legend('unbroken', 'broken', 'Location', 'northwest', 'Orientation', 'vertical');
legend('boxoff');
ylim([0, 20]);
xlim([0.5, 3.5]);
ylabel('number of transferred objects in 2 min');

grid off; box off;
set(gca, 'XGrid', 'off');
set(gca, 'GridColor', 'k');
set(gca, 'GridLineStyle', '--');
set(gca, 'GridAlpha', 0.3);
set(gca, 'FontSize', 13);
set(gca, 'YTick', 0:10:20);
% set(gcf, 'position', [400 400 400 350]);

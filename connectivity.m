clc; clear;
load('Sub_EMG_fMRI_connectivity.mat')
%%
TypeName = 'Spearman';
RC_temp = zeros(length(Roilabel));

for i = eff_sub
    adjust_BOLDdata = [];

    for j = 1:4
        temp_data = (Roi_BOLD_data(i).(TaskName{j})(eff_timepoint, :)) .* T(j, :);
        adjust_BOLDdata = cat(1, adjust_BOLDdata, temp_data);
    end

    RC_temp = RC_temp + (round(atanh(corr(adjust_BOLDdata, 'type', TypeName)), 4));
end

RC_temp = RC_temp ./ numel(eff_sub);

RC = RC_temp;
RC_density = 0;
thershold = 0.99;
den_thershold = 0.3;

while (RC_density < den_thershold)
    RC = RC_temp;
    thershold = thershold - 0.01;
    RC(RC_temp < thershold) = 0;
    G = graph(RC, Roilabel, 'omitselfloops');
    L = numedges(G);
    RC_density = L / (length(Roilabel) * (length(Roilabel) - 1) / 2);
end

pos = [];
lag = 0;
Range = [80 80 80 80];

for k = 1:4
    distanceMatrix = 1 - RC((lag + 1):(lag + countTask(k)), (lag + 1):(lag + countTask(k)));
    distanceMatrix(logical(eye(size(distanceMatrix)))) = 0;
    coordinates = tsne(distanceMatrix, 'Algorithm', 'exact', 'Distance', 'euclidean', 'CacheSize', 'maximal'); %~~~~
    coordinates = rescale(coordinates, Range(k), -Range(k));
    lag = lag + countTask(k);
    pos = cat(1, pos, coordinates);
end

pos = [pos; 25, 25; -25, -25];

lag = 0;
e = 150;
exc = [0 0; e e; -e e; -e -e; e -e];

for kk = 1:5
    pos((lag + 1):(lag + countTask(kk)), :) = pos((lag + 1):(lag + countTask(kk)), :) + exc(kk, :);
    lag = lag + countTask(kk);
end

figure('Color', 'w')
h = plot(G, 'XData', pos(:, 1), 'YData', pos(:, 2));
h.MarkerSize = degree(G) ./ 2 + 1;
h.EdgeColor = "none";
h.EdgeAlpha = .1;
h.EdgeColor = [0.5, 0.5, 0.5];
h.LineWidth = 0.1;
h.NodeCData = [repmat(1, 1, countTask(1)), repmat(2, 1, countTask(2)), repmat(3, 1, countTask(3)), repmat(4, 1, countTask(4)), repmat(5, 1, countTask(5))];
CM = [211, 211, 211
      220 133 38 % F
      27 118 180 % A
      100 57 26 % L
      187 33 34] ./ 255; % M
colormap(CM)
axis off
axis equal
%% ++++
Fs = 1; Fc = 0.1; N = 3; [b, a] = butter(N, Fc / (Fs / 2), 'low');
TD = cell(numel(eff_sub), numel(countTask) - 1);
TD_task = cell(numel(eff_sub), numel(countTask) - 1);
mean_bold = cell(numel(eff_sub), numel(countTask) - 1);

for i = eff_sub

    for j = 1:4
        regoin_bold = [];
        lag = 0;

        for kk = 1:numel(countTask)
            data_temp = Roi_BOLD_data(i).(TaskName{j})(eff_timepoint, (lag + 1):(lag + countTask(kk)));
            data_filtered = filtfilt(b, a, data_temp);
            lag = lag + countTask(kk);
            regoin_bold = cat(2, regoin_bold, mean(data_filtered, 2));
        end

        lags_session = zeros(1, 4);

        for k = 1:4
            [c, lags] = xcov(regoin_bold(:, k + 1), regoin_bold(:, 1));
            maxPoint = find(c == max(c));
            x = [lags(maxPoint - 1), lags(maxPoint), lags(maxPoint + 1)];
            y = [c(maxPoint - 1), c(maxPoint), c(maxPoint + 1)];
            p = polyfit(x, y, 2);
            p_derivative = polyder(p);
            extrema = roots(p_derivative);
            extrema = extrema(imag(extrema) == 0);
            lags_session(k) = extrema;
        end

        mean_bold{i, j} = regoin_bold;
        TD{i, j} = lags_session;

        lags_session_task = zeros(1, 3);
        flag_num = 1;

        for k = setdiff(1:4, j)
            [c, lags] = xcov(regoin_bold(:, k + 1), regoin_bold(:, j + 1));
            maxPoint = find(c == max(c));
            x = [lags(maxPoint - 1), lags(maxPoint), lags(maxPoint + 1)];
            y = [c(maxPoint - 1), c(maxPoint), c(maxPoint + 1)];
            p = polyfit(x, y, 2);
            p_derivative = polyder(p);
            extrema = roots(p_derivative);
            extrema = extrema(imag(extrema) == 0);
            lags_session_task(flag_num) = extrema;
            flag_num = flag_num + 1;
        end

        TD_task{i, j} = lags_session_task;
    end

end

for i = eff_sub
    TD{i, 5} = (TD{i, 1} + TD{i, 2} + TD{i, 3} + TD{i, 4}) ./ 4;
    TD_task{i, 5} = [mean(TD_task{i, 1}), mean(TD_task{i, 2}), mean(TD_task{i, 3}), mean(TD_task{i, 4})];
end

sub_TD = cell2mat(TD(:, 5))';
mean_TD = mean(sub_TD, 2);
[h_TD, p_TD] = ttest(sub_TD');
figure('Units', 'normalized', 'Position', [.2, .3, .36, .45], 'Color', 'w');
CM = [211, 211, 211
      220 133 38 % F
      27 118 180 % A
      100 57 26 % L
      187 33 34] ./ 255; % M
NameList = {'F'; 'A'; 'L'; 'M'};
hold on
barHdl = bar(mean_TD([5, 4, 2, 3] - 1), 'EdgeColor', 'none', 'FaceAlpha', .8, 'BarWidth', .7);
barHdl.FaceColor = 'flat';
barHdl.CData(1:4, :) = CM([5, 4, 2, 3], :);
XA = barHdl.XEndPoints.' * ones(1, size(sub_TD, 2));
YA = sub_TD([5, 4, 2, 3] - 1, :);
scatter(XA(:), YA(:), 55, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', .5, ...
    'CData', repmat(CM([5, 4, 2, 3], :), 9, 1), 'XJitter', 'rand', 'XJitterWidth', 0.2, 'SizeData', 40)

fig = gcf;
ax1 = fig.Children;

ax1.XColor = 'none';
ax1.Children(end).BaseLine.LineWidth = 0.8;
ax1.LineWidth = 0.8;
axis([0, 5, -1, 1]);
ax1.YLabel.String = 'time delays(s)';
ax1.YLabel.FontSize = 15;

sub_TD_task = cell2mat(TD_task(:, 5))';
mean_TD_task = mean(sub_TD_task, 2);
[h_TD, p_TD] = ttest(sub_TD_task');

for i = eff_sub
    TD{i, 6} = [TD{i, 1}(1), TD{i, 2}(2), TD{i, 3}(3), TD{i, 4}(4)];
end

CMPs_TD = cell2mat(TD(:, 6))';
mean_CMPs_TD = mean(CMPs_TD, 2);
figure('Units', 'normalized', 'Position', [.2, .3, .36, .45], 'Color', 'w');
hold on
barHdl = bar(mean_CMPs_TD([5, 4, 2, 3] - 1), 'EdgeColor', 'none', 'FaceAlpha', .8, 'BarWidth', .7);
barHdl.FaceColor = 'flat';
barHdl.CData(1:4, :) = CM([5, 4, 2, 3], :);
XA = barHdl.XEndPoints.' * ones(1, size(CMPs_TD, 2));
YA = CMPs_TD([5, 4, 2, 3] - 1, :);
scatter(XA(:), YA(:), 55, 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', .5, ...
    'CData', repmat(CM([5, 4, 2, 3], :), 9, 1), 'XJitter', 'rand', 'XJitterWidth', 0.2, 'SizeData', 40)

fig = gcf;
ax1 = fig.Children;

ax1.XColor = 'none';
ax1.Children(end).BaseLine.LineWidth = 0.8;
ax1.LineWidth = 0.8;
axis([0, 5, -1, 1]);
ax1.YLabel.String = 'time delays(s)';
ax1.YLabel.FontSize = 15;

[h_TD, p_TD] = ttest(CMPs_TD');
%%
fmri_point = 45:164;
emg_point = 1:300;
ROI_connect = struct('sub', {}, 'fMRI_data_series', {}, 'EMG_data_series', {}, 'normalized_xcov_matrix', {});
mean_nor_matrix = zeros(5, 4);
effnum = 0;
snr_sub = [];

for i = eff_sub
    ROI_connect(i).sub = Roi_BOLD_data(i).sub;
    adjust_BOLDdata = [];

    for j = 1:numel(TaskName)
        temp_data = (Roi_BOLD_data(i).(TaskName{j})(fmri_point, :)) .* T(j, :);
        adjust_BOLDdata = cat(1, adjust_BOLDdata, temp_data);
        data_series = [];
        lag = 0;

        for k = 1:numel(countTask)
            data_series = cat(2, data_series, mean(temp_data(:, (lag + 1):(lag + countTask(k))), 2));
            lag = lag + countTask(k);
        end

        ROI_connect(i).fMRI_data_series = cat(1, ROI_connect(i).fMRI_data_series, data_series);
    end

    ROI_connect(i).fMRI_data_series = normalize(ROI_connect(i).fMRI_data_series); %normalize

    snr_temp = zeros(4);
    windowWidth = 100;

    for j = 1:numel(TaskName)
        roi_emg_data = zeros(numel(fmri_point), 4);
        data_temp = [];

        for k = 1:4
            adjust_EMGdata = smooth(nanmean(Roi_EMG_data(i).(TaskName{j}){k}, 2), windowWidth, 'moving');
            EMGsnr = calculate_series_snr(adjust_EMGdata, windowWidth);
            data_temp = cat(2, data_temp, 10 ^ (EMGsnr) * adjust_EMGdata);
            snr_temp(k, j) = 10 ^ (EMGsnr);
        end

        data_temp = data_temp(emg_point, :);

        for k = 1:4
            roi_emg_data(:, k) = interp1(data_temp(:, k), linspace(1, size(data_temp(:, k), 1), length(temp_data)), 'spline');
        end

        normEMG = mapminmax(roi_emg_data(:)', 0, 1);
        roi_emg_data = reshape(normEMG, size(roi_emg_data));
        ROI_connect(i).EMG_data_series = cat(1, ROI_connect(i).EMG_data_series, roi_emg_data);
    end

    ROI_connect(i).EMG_data_series = normalize(ROI_connect(i).EMG_data_series);
    snr_sub = cat(3, snr_sub, snr_temp);

    ROI_connect(i).normalized_xcov_matrix = zeros(size(ROI_connect(i).fMRI_data_series, 2), size(ROI_connect(i).EMG_data_series, 2));

    for m = 1:size(ROI_connect(i).fMRI_data_series, 2)

        for n = 1:size(ROI_connect(i).EMG_data_series, 2)
            [c, lags] = xcorr(ROI_connect(i).fMRI_data_series(:, m), ROI_connect(i).EMG_data_series(:, n), 60, 'normalized');
            ROI_connect(i).normalized_xcov_matrix(m, n) = max(abs(c));
        end

    end

    mean_nor_matrix = mean_nor_matrix + ROI_connect(i).normalized_xcov_matrix;
    effnum = effnum + 1;
end

mean_nor_matrix = mean_nor_matrix ./ effnum

% cross-modal connectivity
TypeName = 'Spearman';
RC_temp = zeros(length(Roilabel));

for i = eff_sub
    adjust_BOLDdata = [];

    for j = 1:4
        temp_data = (Roi_BOLD_data(i).(TaskName{j})(fmri_point, :)) .* T(j, :);
        adjust_BOLDdata = cat(1, adjust_BOLDdata, temp_data);
    end

    RC_temp = RC_temp + (round((corr(adjust_BOLDdata, 'type', TypeName)), 4));
end

RC_temp = RC_temp ./ numel(eff_sub);

RC = RC_temp;
RC_density = 0;
thershold = 0.99;
den_thershold = 0.4;

while (RC_density < den_thershold)
    RC = RC_temp;
    thershold = thershold - 0.01;
    RC(RC_temp < thershold) = 0;
    G = graph(RC, Roilabel, 'omitselfloops');
    L = numedges(G);
    RC_density = L / (length(Roilabel) * (length(Roilabel) - 1) / 2);
end

tempRC = [];
plotorder = [1 5 4 2 3];
col = [1, 13; 14, 19; 20, 25; 26, 32; 33, 34];

for k = plotorder
    tempRC = cat(1, tempRC, RC(col(k, 1):col(k, 2), :));
end

plotRC = [];

for k = plotorder
    plotRC = cat(2, plotRC, tempRC(:, col(k, 1):col(k, 2)));
end

figure; heatmap(Roilabel, Roilabel, plotRC);
colormap('sky')

pos = [];
lag = 0;
Range = [200 200 200 200];

for k = 1:4
    distanceMatrix = 1 - (RC((lag + 1):(lag + countTask(k)), (lag + 1):(lag + countTask(k))));
    distanceMatrix(logical(eye(size(distanceMatrix)))) = 0;
    coordinates = tsne(distanceMatrix, 'Algorithm', 'exact', 'Distance', 'euclidean', 'CacheSize', 'maximal'); %~~~~
    coordinates = rescale(coordinates, Range(k), -Range(k));
    lag = lag + countTask(k);
    pos = cat(1, pos, coordinates);
end

pos = [pos, zeros(size(pos, 1), 1)];
pos = [pos; 25, 25, 0; -25, -25, 0];

lag = 0;
e = 400;
palmE = 2 * e;

RC_center = zeros(length(countTask));

for i = eff_sub
    RC_center = RC_center + abs(round(atanh(corr(ROI_connect(i).fMRI_data_series, 'type', TypeName)), 4));
end

RC_center = RC_center ./ numel(eff_sub);
RC_center(logical(eye(size(RC_center)))) = 0;
center = tsne(RC_center, 'Algorithm', 'exact', 'Distance', 'chebychev', 'CacheSize', 'maximal', 'Exaggeration', 5, 'Standardize', true);
center = rescale(center, -palmE, palmE);
exc = [center, [0 0 0 0 0]'];

for kk = 1:5
    pos((lag + 1):(lag + countTask(kk)), :) = pos((lag + 1):(lag + countTask(kk)), :) + exc(kk, :);
    lag = lag + countTask(kk);
end

figure('Color', 'w')
h = plot(G, 'XData', pos(:, 1), 'YData', pos(:, 2), 'ZData', pos(:, 3), 'NodeLabel', []);
h.MarkerSize = degree(G) ./ 2 + 1;
h.EdgeAlpha = .1;
h.EdgeColor = [0.2, 0.2, 0.2];
h.LineWidth = 0.25;
h.NodeCData = [repmat(1, 1, countTask(1)), repmat(2, 1, countTask(2)), repmat(3, 1, countTask(3)), repmat(4, 1, countTask(4)), repmat(5, 1, countTask(5))];
CM = [143, 59, 127
      220 133 38 % F
      27 118 180 % A
      100 57 26 % L
      187 33 34] ./ 255; % M
colormap(CM)
axis off
axis equal
hold on

s = scatter3(center(2:end, 1)', center(2:end, 2)', -palmE * ones(size(center, 1) - 1, 1)', 80, CM(2:end, :), "filled", 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75); %EMG ROI center point
s.MarkerFaceAlpha = 0.9;

% connect fMRI ROI center with EMG ROI center using [mean_nor_matrix]
connect_thershold = 0.6;
connect_matrix = mean_nor_matrix;
connect_matrix(connect_matrix < connect_thershold) = 0;

for m = 2:size(connect_matrix, 1)

    for n = 1:size(connect_matrix, 2)

        if connect_matrix(m, n) > 0
            plot3([center(m, 1), center(n + 1, 1)], [center(m, 2), center(n + 1, 2)], [0, -palmE], 'LineWidth', exp(connect_matrix(m, n)), 'Color', [CM(m, :), exp(connect_matrix(m, n)) - 1.5])
        end

    end

end

hold off

hold on
C = mean(center);
EdgeLong = max(pos, [], 'all');
EdgeLong = max(abs(pos), [], 'all') + 0.5 * 10 ^ (floor(log10(EdgeLong)));
[x, y] = meshgrid([C(1) - EdgeLong, C(1) + EdgeLong], [C(2) - EdgeLong, C(2) + EdgeLong]);

z1 = zeros(size(x));
co = 0.9 * ones(size(z1, 1), size(z1, 2), 3);
a = 0.2;
h = surf(x, y, z1, co, 'EdgeColor', 'k');
alpha(h, a);

z2 = -palmE * ones(size(x));
h = surf(x, y, z2, co, 'EdgeColor', 'k');
alpha(h, a);
view(125, 15);

figure('Color', 'w')
h = plot(G, 'XData', pos(:, 1), 'YData', pos(:, 2), 'NodeLabel', []);
h.MarkerSize = degree(G) ./ 2 + 1;
h.EdgeAlpha = .1;
h.EdgeColor = [0.2, 0.2, 0.2];
h.LineWidth = 0.25;
h.NodeCData = [repmat(1, 1, countTask(1)), repmat(2, 1, countTask(2)), repmat(3, 1, countTask(3)), repmat(4, 1, countTask(4)), repmat(5, 1, countTask(5))];
CM = [143, 59, 127
      220 133 38 % F
      27 118 180 % A
      100 57 26 % L
      187 33 34] ./ 255; % M
colormap(CM)
axis off
axis equal

%  snr encode
mean_snrEMGdata = zeros(size(ROI_connect(1).EMG_data_series));

for i = eff_sub
    mean_snrEMGdata = mean_snrEMGdata + ROI_connect(i).EMG_data_series;
end

mean_snrEMGdata = mean_snrEMGdata ./ effnum;
figure('Color', 'w')
adjustorder = [4 3 1 2];
datas_seg = [1:120; 121:240; 241:360; 361:480];
index = reshape(1:16, 4, 4)';
datacell = cell(4);

for m = 1:4

    for n = 1:4
        subplot(4, 4, index(m, n))
        plot(mean_snrEMGdata(datas_seg(adjustorder(n), :), adjustorder(m)));
        ylim([-1, 2])
        box off
        axis off
        datacell{m, n} = mean_snrEMGdata(datas_seg(adjustorder(n), :), adjustorder(m));
    end

end

snrMatrix = zeros(4);

for m = 1:4

    for n = 1:4
        snrMatrix(m, n) = snr_contrast(datacell{m, m}, datacell{m, n});
    end

end

figure('Color', 'w')
imagesc(snrMatrix);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
axis equal
box on
axis off
colormap((slanCM('Greens')))
colorbar
%%
%  normalized cross-correlation coefficient between EMG and fMRI data
fmri_point = 45:164;
emg_point = 1:300;
ROI_connect = struct('sub', {}, 'fMRI_data_series', {}, 'EMG_data_series', {}, 'normalized_xcov_matrix', {});
mean_nor_matrix = zeros(5, 4);

for i = eff_sub
    ROI_connect(i).sub = Roi_BOLD_data(i).sub;
    adjust_BOLDdata = [];

    for j = 1:numel(TaskName)
        temp_data = (Roi_BOLD_data(i).(TaskName{j})(fmri_point, :));
        adjust_BOLDdata = cat(1, adjust_BOLDdata, temp_data);
        data_series = [];
        lag = 0;

        for k = 1:numel(countTask)
            data_series = cat(2, data_series, mean(temp_data(:, (lag + 1):(lag + countTask(k))), 2));
            lag = lag + countTask(k);
        end

        ROI_connect(i).fMRI_data_series = cat(1, ROI_connect(i).fMRI_data_series, data_series);
    end

    ROI_connect(i).fMRI_data_series = normalize(ROI_connect(i).fMRI_data_series); %normalize

    for j = 1:numel(TaskName)
        roi_emg_data = zeros(numel(fmri_point), 4);
        data_temp = [];

        for k = 1:4
            adjust_EMGdata = nanmean(Roi_EMG_data(i).(TaskName{j}){k}, 2);
            data_temp = cat(2, data_temp, adjust_EMGdata);
        end

        data_temp = data_temp(emg_point, :);

        for k = 1:4
            roi_emg_data(:, k) = interp1(data_temp(:, k), linspace(1, size(data_temp(:, k), 1), length(temp_data)), 'spline');
        end

        normEMG = mapminmax(roi_emg_data(:)', 0, 1);
        roi_emg_data = reshape(normEMG, size(roi_emg_data));
        ROI_connect(i).EMG_data_series = cat(1, ROI_connect(i).EMG_data_series, roi_emg_data);
    end

    ROI_connect(i).EMG_data_series = normalize(ROI_connect(i).EMG_data_series);

    ROI_connect(i).normalized_xcov_matrix = zeros(size(ROI_connect(i).fMRI_data_series, 2), size(ROI_connect(i).EMG_data_series, 2));

    for m = 1:size(ROI_connect(i).fMRI_data_series, 2)

        for n = 1:size(ROI_connect(i).EMG_data_series, 2)
            [c, lags] = xcorr(ROI_connect(i).fMRI_data_series(:, m), ROI_connect(i).EMG_data_series(:, n), 60, 'normalized');
            ROI_connect(i).normalized_xcov_matrix(m, n) = max(abs(c));
        end

    end

    mean_nor_matrix = mean_nor_matrix + ROI_connect(i).normalized_xcov_matrix;
end

mean_nor_matrix = mean_nor_matrix ./ numel(eff_sub)

TypeName = 'Spearman';
RC_temp = zeros(length(Roilabel));

for i = eff_sub
    adjust_BOLDdata = [];

    for j = 1:4
        temp_data = (Roi_BOLD_data(i).(TaskName{j})(fmri_point, :));
        adjust_BOLDdata = cat(1, adjust_BOLDdata, temp_data);
    end

    RC_temp = RC_temp + (round((corr(adjust_BOLDdata, 'type', TypeName)), 4));
end

RC_temp = RC_temp ./ numel(eff_sub);

RC = RC_temp;
RC_density = 0;
thershold = 0.99;
den_thershold = 0.4;

while (RC_density < den_thershold)
    RC = RC_temp;
    thershold = thershold - 0.01;
    RC(RC_temp < thershold) = 0;
    G = graph(RC, Roilabel, 'omitselfloops');
    L = numedges(G);
    RC_density = L / (length(Roilabel) * (length(Roilabel) - 1) / 2);
end

pos = [];
numDimensions = 3;
lag = 0;
Range = [200 200 200 200];

for k = 1:4
    distanceMatrix = 1 - (RC((lag + 1):(lag + countTask(k)), (lag + 1):(lag + countTask(k))));
    distanceMatrix(logical(eye(size(distanceMatrix)))) = 0;
    coordinates = tsne(distanceMatrix, 'Algorithm', 'exact', 'Distance', 'euclidean', 'CacheSize', 'maximal'); %~~~~
    coordinates = rescale(coordinates, Range(k), -Range(k));
    lag = lag + countTask(k);
    pos = cat(1, pos, coordinates);
end

pos = [pos, zeros(size(pos, 1), 1)];
pos = [pos; 25, 25, 0; -25, -25, 0];

lag = 0;
e = 400;
palmE = 2 * e;

RC_center = zeros(length(countTask));

for i = eff_sub
    RC_center = RC_center + abs(round((corr(ROI_connect(i).fMRI_data_series, 'type', TypeName)), 4));
end

RC_center = RC_center ./ numel(eff_sub);
RC_center(logical(eye(size(RC_center)))) = 0;
center = tsne(RC_center, 'Algorithm', 'exact', 'Distance', 'chebychev', 'CacheSize', 'maximal', 'Exaggeration', 5, 'Standardize', true);
center = rescale(center, -palmE, palmE);
exc = [center, [0 0 0 0 0]'];

for kk = 1:5
    pos((lag + 1):(lag + countTask(kk)), :) = pos((lag + 1):(lag + countTask(kk)), :) + exc(kk, :);
    lag = lag + countTask(kk);
end

figure('Color', 'w')
h = plot(G, 'XData', pos(:, 1), 'YData', pos(:, 2), 'ZData', pos(:, 3), 'NodeLabel', []);
h.MarkerSize = degree(G) ./ 2 + 1;
h.EdgeAlpha = .1;
h.EdgeColor = [0.2, 0.2, 0.2];
h.LineWidth = 0.25;
h.NodeCData = [repmat(1, 1, countTask(1)), repmat(2, 1, countTask(2)), repmat(3, 1, countTask(3)), repmat(4, 1, countTask(4)), repmat(5, 1, countTask(5))];
CM = [143, 59, 127
      220 133 38 % F
      27 118 180 % A
      100 57 26 % L
      187 33 34] ./ 255; % M
colormap(CM)
axis off
axis equal
hold on

s = scatter3(center(2:end, 1)', center(2:end, 2)', -palmE * ones(size(center, 1) - 1, 1)', 80, CM(2:end, :), "filled", 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75); %EMG ROI center point
s.MarkerFaceAlpha = 0.9;

% connect fMRI ROI center with EMG ROI center using [mean_nor_matrix]
connect_thershold = 0;
connect_matrix = mean_nor_matrix;
connect_matrix(connect_matrix < connect_thershold) = 0;

for m = 2:size(connect_matrix, 1)

    for n = 1:size(connect_matrix, 2)

        if connect_matrix(m, n) > 0
            plot3([center(m, 1), center(n + 1, 1)], [center(m, 2), center(n + 1, 2)], [0, -palmE], 'LineWidth', exp(5 * connect_matrix(m, n)), 'Color', [CM(m, :), exp(connect_matrix(m, n)) - 1])
        end

    end

end

hold off

hold on
C = mean(center);
EdgeLong = max(pos, [], 'all');
EdgeLong = max(pos, [], 'all') + 0.5 * 10 ^ (floor(log10(EdgeLong)));
[x, y] = meshgrid([C(1) - EdgeLong, C(1) + EdgeLong], [C(2) - EdgeLong, C(2) + EdgeLong]);

z1 = zeros(size(x));
co = 0.9 * ones(size(z1, 1), size(z1, 2), 3);
a = 0.2;
h = surf(x, y, z1, co, 'EdgeColor', 'k');
alpha(h, a);
z2 = -palmE * ones(size(x));
h = surf(x, y, z2, co, 'EdgeColor', 'k');
alpha(h, a);
view(125, 15);
hold off

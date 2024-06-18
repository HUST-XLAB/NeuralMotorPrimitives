clear
clc
config = cosmo_config();
study_path = fullfile(config.tutorial_data_path, 'ak6');
output_path = config.output_data_path;

% reset citation list
cosmo_check_external('-tic');

%% Preprocessing for DISTATIS: RSM analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubNum = {
          'HUST03' %1
          'HUST04' %2-
          'HUST05' %3-
          'HUST06' %4-
          'HUST07' %5
          'HUST08' %6
          'HUST09' %7-
          'HUST10' %8-
          'HUST12' %9-
          'HUST13' %10-
          'HUST14' %11
          'HUST15' %12-
          'HUST16' %13-
          'HUST17' %14-
          };

mask_path = 'U:\Data Base\fMRI\DataSet\data_processing\mask\';
mask_fn = fullfile(mask_path, 'rM1_mask.nii');
DirHead = 'U:\Data Base\fMRI\DataSet\data_processing\';

ds_rsms = cell(8, 1);

task = [1, 6, 7, 8];

for i = 1:8

    for j = 1:4

        TaskName = strcat('task', num2str(task(j), '%02d'));

        SpmPath = strcat(DirHead, SubNum{n(i), 1}, '\task\nii\', TaskName, '\All\classical-C6-LR-M1\');
        data_fn = fullfile(SpmPath, 'SPM.mat');
        ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
        ds.samples(j, :) = [1, 0] * ds_temp.samples;
        ds.sa.labels{j, 1} = cellstr(TaskName);
        ds.sa.targets(j, 1) = j;
        ds.sa.chunks(j, 1) = 1;
        ds.a = ds_temp.a;
        ds.fa = ds_temp.fa;
    end

    ds_rsm = cosmo_dissimilarity_matrix_measure(ds);
    ds_rsm.samples = pdist(ds.samples, 'correlation')';
    ds_rsm.sa.chunks = i * ones(size(ds_rsm.samples, 1), 1);
    ds_rsms{i} = ds_rsm;
end

% combine data from all subjects
all_ds = cosmo_stack(ds_rsms);

distatis = cosmo_distatis(all_ds);

%% show comprimise distance matrix
[compromise_matrix, dim_labels, values] = cosmo_unflatten(distatis, 1);

labels = {'F', 'A', 'L', 'M'};
n_labels = numel(labels);
figure();

load('Normal_M4.mat');
imagesc(M4MappedData, [0, 1])

title('DSM');
set(gca, 'YTick', 1:n_labels, 'YTickLabel', labels);
set(gca, 'XTick', 1:n_labels, 'XTickLabel', labels);
ylabel(dim_labels{1});
xlabel(dim_labels{2});

load('colorbar_imola.mat');
colormap(aa);

colorbar

if cosmo_check_external('@stats', false)
    figure();
    hclus = linkage(compromise_matrix);
    dendrogram(hclus, 'labels', labels, 'orientation', 'left');
    title('dendrogram');

    figure();
    F = cmdscale(squareform(compromise_matrix));
    text(F(:, 1), F(:, 2), labels);
    title('2D MDS plot');
    mx = max(abs(F(:)));
    xlim([-mx mx]); ylim([-mx mx]);
end

cosmo_check_external('-cite');

% [compromise_matrix,F]=dist_all([1 6 7 8],'rM1_mask.nii')
% [compromise_matrix,F]=dist_all(1:5,'rM1_mask.nii')
% [compromise_matrix,F]=dist_all(1:8,'rM1_mask.nii')
%% Set data paths
% The function cosmo_config() returns a struct containing paths to tutorial
% data. (Alternatively the paths can be set manually without using
% cosmo_config.)
function [compromise_matrix, F] = dist_all(task_index, mask)
    config = cosmo_config();
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

    output_path = 'U:\Data Base\fMRI\DataSet\data_processing\mvpa\dist';
    DirHead = 'U:\Data Base\fMRI\DataSet\data_processing\mvpa\svm\nii\';
    mask_path = 'U:\Data Base\fMRI\DataSet\data_processing\mask\';
    mask_fn = fullfile(mask_path, mask);

    % reset citation list
    cosmo_check_external('-tic');
    %% Preprocessing for DISTATIS: RSM analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    task_class = length(task_index);
    ds_rsms = cell(10, 1); % allocate space for output

    for i = 1:10

        for j = 1:task_class

            if (task_class == 8) && (j < 6)
                TaskName = strcat('task', num2str(task_index(j)));

                data_fn = strcat(DirHead, SubNum{n(i), 1}, '_class', num2str(task_class), '_', TaskName, '_rM1.nii');
                ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
                ds.samples(j, :) = ds_temp.samples;
                ds.sa.labels{j, 1} = cellstr(TaskName);
                ds.sa.targets(j, 1) = j;
                ds.sa.chunks(j, 1) = 1;
                ds.a = ds_temp.a;
                ds.fa = ds_temp.fa;
            else
                TaskName = strcat('task', num2str(task_index(j)));
                data_fn = strcat(DirHead, SubNum{n(i), 1}, '_class', num2str(task_class), '_', TaskName, '_rM1.nii');
                ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
                ds.samples(j, :) = ds_temp.samples;
                ds.sa.labels{j, 1} = cellstr(TaskName);
                ds.sa.targets(j, 1) = j;
                ds.sa.chunks(j, 1) = 1;
                ds.a = ds_temp.a;
                ds.fa = ds_temp.fa;
            end

        end

        ds_rsm = cosmo_dissimilarity_matrix_measure(ds);
        ds_rsm.samples = pdist(ds.samples, 'cosine')';
        ds_rsm.sa.chunks = i * ones(size(ds_rsm.samples, 1), 1);
        ds_rsms{i} = ds_rsm;
    end

    % combine data from all subjects
    all_ds = cosmo_stack(ds_rsms);
    %all_ds.samples=(ds_rsms{1,1}.samples+ds_rsms{2,1}.samples+ds_rsms{3,1}.samples+ds_rsms{4,1}.samples)/4;
    temp = 0;

    for n = 1:10
        temp = temp + ds_rsms{n, 1}.samples;
    end

    all_ds.samples = temp ./ 10;

    %% Run DISTATIS
    distatis = cosmo_distatis(all_ds);

    %% show comprimise distance matrix
    [compromise_matrix, dim_labels, values] = cosmo_unflatten(distatis, 1);
    compromise_matrix = squareform(all_ds.samples);

    switch task_class
        case 4
            labels = {'F', 'A', 'L', 'M'};
        case 5
            labels = {'Thumb', 'Index finger', 'Middle finger', 'Ring finger', 'Little finger'};
        case 8
            labels = {'Little finger', 'Ring finger', 'Middle finger', 'Index finger', 'Thumb/F', 'A', 'L', 'M'};
    end

    n_labels = numel(labels);
    maskname = extractBefore(mask, '_mask.nii');
    figure();
    imagesc(compromise_matrix);
    title('DSM');
    set(gca, 'YTick', 1:n_labels, 'YTickLabel', labels);
    set(gca, 'XTick', 1:n_labels, 'XTickLabel', labels);
    ylabel(dim_labels{1});
    xlabel(dim_labels{2});
    colorbar
    %savefig([output_path,'\All-DSM_','_class',num2str(task_class),'_',maskname,'.fig']);
    %close
    % skip if stats toolbox is not present
    if cosmo_check_external('@stats', false)
        figure();
        hclus = linkage(compromise_matrix); %
        dendrogram(hclus, 'labels', labels, 'orientation', 'left');
        title('dendrogram');
        xlim([0 inf]);
        %savefig([output_path,'\All-Dendrogram_','_class',num2str(task_class),'_',maskname,'.fig']);
        %close
        figure();
        F = cmdscale(squareform(compromise_matrix));
        text(F(:, 1), F(:, 2), labels);
        title('2D MDS plot');
        mx = max(abs(F(:)));
        xlim([-mx mx]); ylim([-mx mx]);
        %savefig([output_path,'\All-MDS_','_class',num2str(task_class),'_',maskname,'.fig']);
        %close
    end

    return

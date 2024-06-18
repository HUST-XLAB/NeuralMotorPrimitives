function [rsm, dis, compromise, Allcompromise] = dist_SLmask1678(task_index)
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
    DirHead = 'U:\Data Base\fMRI\DataSet\data_processing\';
    mask_path = 'U:\Data Base\fMRI\DataSet\data_processing\mask\';
    M1mask_fn = fullfile(mask_path, 'rM1_mask.nii');
    mask_fn = fullfile(mask_path, 'T1678_mask.nii');
    % reset citation list
    cosmo_check_external('-tic');

    task_class = length(task_index);
    [Allcompromise, fm] = dist_all(task_index, 'rM1_mask.nii');
    figure();

    for task_num = 1:task_class

        tmask_fn = strcat(mask_path, 'T', num2str(task_index(task_num)), '_mask.nii');
        test = cosmo_fmri_dataset(tmask_fn, 'mask', mask_fn);

        for i = 1:4

            for j = 1:10
                TaskName = strcat('task', num2str(task_index(task_num), '%02d'));

                SpmPath = strcat(DirHead, SubNum{n(j), 1}, '\task\nii\', TaskName, '\All\split4-1\', num2str(i));
                data_fn = fullfile(SpmPath, 'SPM.mat');

                ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
                ds.samples(j, :) = ds_temp.samples .* test.samples;
                ds.sa.labels{j, 1} = cellstr(TaskName);
                ds.sa.targets(j, 1) = j;
                ds.sa.chunks(j, 1) = 1;
                ds.a = ds_temp.a;
                ds.fa = ds_temp.fa;
            end

            ds_rsm = cosmo_dissimilarity_matrix_measure(ds);
            ds_rsm.samples = pdist(ds.samples, 'cosine')';
            ds_rsm.sa.chunks = ones(size(ds_rsm.samples, 1), 1);
            ds_rsms{i} = ds_rsm;
        end

        all_ds = cosmo_stack(ds_rsms);
        all_ds = (ds_rsms{1, 1}.samples + ds_rsms{1, 2}.samples + ds_rsms{1, 3}.samples + ds_rsms{1, 4}.samples) / 4;
        all_ds = ds_rsm;
        all_ds.samples = (ds_rsms{1, 1}.samples + ds_rsms{1, 2}.samples + ds_rsms{1, 3}.samples + ds_rsms{1, 4}.samples) / 4
        %% Run DISTATIS
        rsm(task_num).data = ds_rsms;

        distatis = cosmo_distatis(all_ds);
        dis(task_num).data = distatis;
        %% show comprimise distance matrix
        [compromise_matrix, dim_labels, values] = cosmo_unflatten(distatis, 1);
        compromise(task_num).data = compromise_matrix;

        labels = {'F', 'A', 'L', 'M'};

        if cosmo_check_external('@stats', false)

            F = cmdscale(squareform(compromise_matrix));
            Y = [F(:, 1), F(:, 2)];
            C = [[18 125 194] ./ 255; [230 137 33] ./ 255; [110 58 27] ./ 255; [198 29 32] ./ 255];
            ConfidenceRegion(Y, 10 * fm(task_num, 1), 10 * fm(task_num, 2), 0.05, 'exp')
            text(10 * fm(task_num, 1), 10 * fm(task_num, 2), labels{task_num});

            hold on

        end

    end

    title('2D MDS plot');
    return

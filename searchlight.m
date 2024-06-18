function ouput = searchlight(subname, task_index, mask)
    DirHead = 'U:\Data Base\fMRI\DataSet\data_processing\';
    mask_path = 'U:\Data Base\fMRI\DataSet\data_processing\mask\';
    mask_fn = fullfile(mask_path, mask);
    %% Set data paths
    cosmo_check_external('-tic');

    %% LDA classifier searchlight analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    task_class = length(task_index);

    for i = 1:3

        for j = 1:task_class
            TaskName = strcat('task', num2str(task_index(j), '%02d'));
            SpmPath = strcat(DirHead, subname, '\task\nii\', TaskName, '\All\split4-1\', num2str(i));
            data_fn = fullfile(SpmPath, 'SPM.mat');
            ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
            ds.samples(j + (i - 1) * task_class, :) = ds_temp.samples;
            ds.sa.labels{j + (i - 1) * task_class, 1} = cellstr(strcat(TaskName, ' for test ', num2str(i)));
            ds.sa.targets(j + (i - 1) * task_class, 1) = j;
            ds.sa.chunks(j + (i - 1) * task_class, 1) = j + (i - 1) * task_class;
        end

    end

    ds.a = ds_temp.a;
    ds.fa = ds_temp.fa;
    train_ds = ds;

    measure = @cosmo_crossvalidation_measure;
    measure_args = struct();

    measure_args.classifier = @cosmo_classify_svm;

    %measure_args.partitions=cosmo_nchoosek_partitioner(ds,1);
    measure_args.partitions = cosmo_independent_samples_partitioner(ds, 'test_count', 1, 'fold_count', 20 * task_class);
    test_index = ((3 * task_class + 1):4 * task_class)'; %

    for k = 1:length(measure_args.partitions.train_indices)
        measure_args.partitions.test_indices{1, k} = test_index;
    end

    for t = 1:task_class
        CV_task = task_index(t);
        test_ds = train_ds;

        for j = 1:task_class
            TaskName = strcat('task', num2str(CV_task, '%02d'));
            SpmPath = strcat(DirHead, subname, '\task\nii\', TaskName, '\All\split4-1\', num2str(4));
            data_fn = fullfile(SpmPath, 'SPM.mat');
            ds_temp = cosmo_fmri_dataset(data_fn, 'mask', mask_fn);
            test_ds.samples(j + 3 * task_class, :) = ds_temp.samples;
            test_ds.sa.labels{j + 3 * task_class, 1} = cellstr(strcat(TaskName, ' for CV'));
            test_ds.sa.targets(j + 3 * task_class, 1) = t;
            test_ds.sa.chunks(j + 3 * task_class, 1) = j + 3 * task_class;
        end

        nvoxels_per_searchlight = 6;
        nbrhood = cosmo_spherical_neighborhood(test_ds, ...

            lda_results = cosmo_searchlight(test_ds, nbrhood, measure, measure_args);
        ouput(t).info = cellstr(strcat('Class ', num2str(task_class), ', ', TaskName, ' for CV'));
        ouput(t).results = lda_results;

    end

    return

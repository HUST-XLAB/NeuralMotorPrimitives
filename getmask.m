function getmask(filename)
    mask_path = 'U:\Data Base\fMRI\DataSet\data_processing\mask\';
    file_fn = fullfile(mask_path, '1_5\', filename);
    mask_fn = fullfile(mask_path, 'rM1_mask.nii');
    test = cosmo_fmri_dataset(file_fn, 'mask', mask_fn);
    temp = test.samples;
    test.samples(temp < 6.180659) = 0;
    test.samples(temp > 6.18) = 1;
    maskname = extractBefore(filename, '_out.nii')
    cosmo_map2fmri(test, strcat(mask_path, '1_5\', maskname, '_mask.nii'))
    return

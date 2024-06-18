%-----------------------------------------------------------------------
% Job saved on 10-Sep-2023 16:18:52 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.write.subj.def = '<UNDEFINED>';
% {'U:\Data Base\fMRI\DataSet\ZIC_data_processing\HUST04 CAI JUNJIE\task\nii\task01\ALL\structure\y_sHUST04-02.nii'};
%%
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = '<UNDEFINED>';

matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.', 'val', '{}', {1}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), substruct('()', {1}, '.', 'files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{3}.spm.stats.fmri_spec.dir = '<UNDEFINED>';
% {'U:\Data Base\fMRI\DataSet\ZIC_data_processing\HUST03 SUN BAIYANG\task\nii\task01\ALL\classical-All'};
matlabbatch{3}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{3}.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch{3}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{3}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{3}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.', 'val', '{}', {2}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), substruct('.', 'files'));
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).name = 'RAction';
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).onset = [5
                                                         45
                                                         85
                                                         125];
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).duration = 20;
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).name = 'LAction';
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).onset = [169
                                                         209
                                                         249
                                                         289];
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).duration = 20;
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{3}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
matlabbatch{3}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{3}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{3}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{3}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{3}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{3}.spm.stats.fmri_spec.volt = 1;
matlabbatch{3}.spm.stats.fmri_spec.global = 'None';
matlabbatch{3}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{3}.spm.stats.fmri_spec.mask = {''};
matlabbatch{3}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{4}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.', 'val', '{}', {3}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), substruct('.', 'spmmat'));
matlabbatch{4}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{4}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{5}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.', 'val', '{}', {4}, '.', 'val', '{}', {1}, '.', 'val', '{}', {1}), substruct('.', 'spmmat'));
matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'RAction> Rest';
matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = [1 0];
matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{2}.tcon.name = 'Rest >RAction';
matlabbatch{5}.spm.stats.con.consess{2}.tcon.weights = [-1 0];
matlabbatch{5}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{3}.tcon.name = 'LAction> Rest';
matlabbatch{5}.spm.stats.con.consess{3}.tcon.weights = [0 1];
matlabbatch{5}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{4}.tcon.name = 'Rest >LAction';
matlabbatch{5}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
matlabbatch{5}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.delete = 0;

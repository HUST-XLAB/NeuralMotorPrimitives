
% 实施例 compromise_matrix=dist_SLmask1_5(1:5)
%% Set data paths
% The function cosmo_config() returns a struct containing paths to tutorial
% data. (Alternatively the paths can be set manually without using
% cosmo_config.)
function compromise_matrix=dist_SLmask1_5(task_index)
config=cosmo_config();
 SubNum={
      'HUST03'%1
      'HUST04'%2-
      'HUST05'%3-
      'HUST06'%4-
      'HUST07'%5
      'HUST08'%6
      'HUST09'%7-
      'HUST10'%8-
      'HUST12'%9-
      'HUST13'%10-
      'HUST14'%11
      'HUST15'%12-
      'HUST16'%13-
      'HUST17'%14-
      };
 % 只需要4,5,6,9,10,12,13,15,16,17号受试者
  n=[2,3,4,7,8,9,10,12,13,14];%
output_path='U:\Data Base\fMRI\DataSet\data_processing\mvpa\dist';
DirHead = 'U:\Data Base\fMRI\DataSet\data_processing\';  
mask_path='U:\Data Base\fMRI\DataSet\data_processing\mask\';
M1mask_fn=fullfile(mask_path,'rM1_mask.nii');  %重采样后的mask文件地址
mask_fn=fullfile(mask_path,'T1_5_mask.nii');  %%这里用了1-5公共mask
% reset citation list
cosmo_check_external('-tic');
%% Preprocessing for DISTATIS: RSM analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 准备数据
%% 读10个人的同一动作，根据人分类，度量不同人之间的距离
task_class=length(task_index);
[~,fm]=dist_all(task_index,'rM1_mask.nii');
figure();



for task_num=[1,2,3]  % 4画不了，5也画不了
%for task_num=1:task_class %原始程序
    
     tmask_fn=strcat(mask_path,'1_5\T',num2str(task_index(task_num)),'_mask.nii');
     test=cosmo_fmri_dataset(tmask_fn,'mask', mask_fn);
     for i=1:4  
         for j=1:10 %10个人
            TaskName=strcat('task',num2str(task_index(task_num),'%02d'));
            %data_fn=strcat(DirHead,SubNum{n(j),1},'_class',num2str(task_class),'_',TaskName,'_rM1.nii');
            SpmPath=strcat(DirHead,SubNum{n(j),1},'\task\nii\',TaskName,'\All\split4-1\',num2str(i));
            data_fn=fullfile(SpmPath,'SPM.mat');
            
            ds_temp=cosmo_fmri_dataset(data_fn,'mask', mask_fn);
            ds.samples(j,:)=ds_temp.samples.*test.samples;
            ds.sa.labels{j, 1}=cellstr(TaskName);
            ds.sa.targets(j,1)= j;
            ds.sa.chunks(j,1)= 1;
            ds.a=ds_temp.a;
            ds.fa=ds_temp.fa;
        end
        ds_rsm=cosmo_dissimilarity_matrix_measure(ds); % 形成距离度量下三角阵，分别以label和值存储，6*6矩阵对应6*5/2=15个比较对
        ds_rsm.samples= pdist(ds.samples,'cosine')';
        ds_rsm.sa.chunks=ones(size(ds_rsm.samples,1),1);
        ds_rsms{i}=ds_rsm;
     end
    all_ds=cosmo_stack(ds_rsms);
    all_ds=(ds_rsms{1,1}.samples+ds_rsms{1,2}.samples+ds_rsms{1,3}.samples+ds_rsms{1,4}.samples)/4;
    all_ds=ds_rsm;
    all_ds.samples=(ds_rsms{1,1}.samples+ds_rsms{1,2}.samples+ds_rsms{1,3}.samples+ds_rsms{1,4}.samples)/4
        %% Run DISTATIS
        distatis=cosmo_distatis(all_ds);

        %% show comprimise distance matrix
        [compromise_matrix,dim_labels,values]=cosmo_unflatten(distatis,1); %将得到的距离阵写成矩阵的形式，更多功能有待探讨！！
        % 10个人
       %labels={'F','A', 'L', 'M'};
       labels={'Thumb', 'Index finger', 'Middle finger', 'Ring finger','Little finger'};

% % %maskname=extractBefore(mask,'_mask.nii');
% % figure();
% % imagesc(compromise_matrix);
% % title('DSM');
% % set(gca,'YTick',1:n_labels,'YTickLabel',labels);
% % set(gca,'XTick',1:n_labels,'XTickLabel',labels);
% % ylabel(dim_labels{1});
% % xlabel(dim_labels{2});
% % colorbar
% % %savefig([output_path,'\DSM_',sub_id,'_class',num2str(task_class),'_',maskname,'.fig']);
% % %close 
% skip if stats toolbox is not present
    if cosmo_check_external('@stats',false)
% %         figure();
% %         hclus = linkage(compromise_matrix); % 层次聚类树
% %         dendrogram(hclus,'labels',labels,'orientation','left');
% %         title('dendrogram');
% %         xlim([0 inf]);
% %         %savefig([output_path,'\Dendrogram_',sub_id,'_class',num2str(task_class),'_',maskname,'.fig']);
% %         %close
% %         figure();
% %         F = cmdscale(squareform(compromise_matrix)); %多维缩放降维
% %         text(F(:,1), F(:,2), labels);
% %         title('2D MDS plot');
% %         mx = max(abs(F(:)));
% %         xlim([-mx mx]); ylim([-mx mx]);


        F = cmdscale(squareform(compromise_matrix));
        Y=[F(:,1), F(:,2)];
        % for t=1:5
        % text(10*fm(t,1),10*fm(t,2), labels{t});
        % end
        % hold on
        ConfidenceRegion(Y,10*fm(task_num,1),10*fm(task_num,2),0.05,'exp')
        %text(10*fm(task_num,1),10*fm(task_num,2), labels{task_num});
%          ConfidenceRegion(Y,fm(task_num,1),fm(task_num,2),0.05,'exp')
%         text(fm(task_num,1),fm(task_num,2), labels{task_num});
        hold on

    %savefig([output_path,'\MDS_',sub_id,'_class',num2str(task_class),'_',maskname,'.fig']);
    %close
    end

end
 hold on
       for t=1:5
        text(10*fm(t,1),10*fm(t,2), labels{t});
        end
xlim([-8 10]) 
ylim([-7 10])
    title('2D MDS plot');
return
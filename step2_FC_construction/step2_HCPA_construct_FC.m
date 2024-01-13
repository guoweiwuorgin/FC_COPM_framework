clear;clc;
addpath('/home/wuguowei/code/task_predict_code');
addpath(genpath('/root_dir/wuguowei/python_code/cifti-matlab/'));
rest_conn_dir = '/root_dir/wuguowei/code/';
%%
template_dir = '/root_dir/wuguowei/python_code/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/';
scheafer_400 = cifti_read([template_dir 'Schaefer2018_200Parcels_17Networks_order.dscalar.nii']);% 200 and 400, 200 for replication

rest_data = '/root_dir/wuguowei/code/HCPA/';
all_subject = dir([rest_data 'sub-*']);
all_subject = {all_subject(:).name}';
subject_list = importdata('/root_dir/wuguowei/code/FC_ontology/FC_ontology_subjects.csv');
subject_name = subject_list.data(:,1);
out_dir = '/root_dir/wuguowei/code/FC_ontology/FC_net_SP_200';
if ~exist(out_dir,'dir')
	mkdir(out_dir);
end
n_iter_2 = 1;
n_iter_1 = 1;
scheafer_net_parcel = scheafer_400.cdata;
for sub = 1:numel(subject_name)
    subname = ['sub-' num2str(subject_name(sub))];
    sub_dir = [rest_data subname filesep];
    rest1 = dir([sub_dir '*-REST1*space-fsLR_den-91k_desc-residual_smooth_den-91k_bold.dtseries.nii']);
    rest2 = dir([sub_dir '*-REST2*space-fsLR_den-91k_desc-residual_smooth_den-91k_bold.dtseries.nii']);
    all_label=scheafer_net_parcel;
    if numel(rest1)==2 && ~exist([out_dir filesep subname filesep  'schaefer_400_run1.mat'],'file') && numel(rest2)==2 && ~exist([out_dir filesep subname filesep  'schaefer_400_run2.mat'],'file')
       rest1_1 = cifti_read([rest1(1).folder filesep rest1(1).name]);
       rest1_1_L = cifti_struct_dense_extract_surface_data(rest1_1,'CORTEX_LEFT');
       rest1_1_R = cifti_struct_dense_extract_surface_data(rest1_1,'CORTEX_RIGHT');
       rest1_1 =[rest1_1_L;rest1_1_R];
      
       rest1_1 = compute_mat_base_label(rest1_1,unique(all_label),all_label);
       %
       rest1_2 = cifti_read([rest1(2).folder filesep rest1(2).name]); 
       rest1_2_L = cifti_struct_dense_extract_surface_data(rest1_2,'CORTEX_LEFT');
       rest1_2_R = cifti_struct_dense_extract_surface_data(rest1_2,'CORTEX_RIGHT');
       rest1_2 =[rest1_2_L;rest1_2_R];
       rest1_2 = compute_mat_base_label(rest1_2,unique(all_label),all_label);
       rest1 = [rest1_1,rest1_2];
       rest1_conn = convet_matrix_to_vector(atanh(corr(rest1(2:end,:)')));
       out_rest1 = [out_dir filesep subname filesep  'schaefer_400_run1.mat'];
        mkdir([out_dir filesep subname filesep]);
        save(out_rest1,'rest1_conn');

       rest2_1 = cifti_read([rest2(1).folder filesep rest2(1).name]);
       rest2_1_L = cifti_struct_dense_extract_surface_data(rest2_1,'CORTEX_LEFT');
       rest2_1_R = cifti_struct_dense_extract_surface_data(rest2_1,'CORTEX_RIGHT');
       rest2_1 =[rest2_1_L;rest2_1_R];
       rest2_1 = compute_mat_base_label(rest2_1,unique(all_label),all_label);
       %
       rest2_2 = cifti_read([rest2(2).folder filesep rest2(2).name]); 
       rest2_2_L = cifti_struct_dense_extract_surface_data(rest2_2,'CORTEX_LEFT');
       rest2_2_R = cifti_struct_dense_extract_surface_data(rest2_2,'CORTEX_RIGHT');
       rest2_2 =[rest2_2_L;rest2_2_R];
       rest2_2 = compute_mat_base_label(rest2_2,unique(all_label),all_label);
       rest2 = [rest2_1,rest2_2];
       rest2_conn = convet_matrix_to_vector(atanh(corr(rest2(2:end,:)')));
       if ~exist([out_dir filesep subname],'dir')
            mkdir([out_dir filesep subname]);
       end
        all_rest = [rest1 rest2];
        rest_conn = convet_matrix_to_vector(atanh(corr(all_rest(2:end,:)')));
       out_rest2 = [out_dir filesep subname filesep  'schaefer_400_run2.mat'];
       save(out_rest2,'rest2_conn');
       out_rest= [out_dir filesep subname filesep  'schaefer_400.mat'];
       save(out_rest,'rest_conn');
    else
        no_file_sub{n_iter_2,2} = subname;
        n_iter_2 = n_iter_2+1;
    end  
    disp(subname);
end




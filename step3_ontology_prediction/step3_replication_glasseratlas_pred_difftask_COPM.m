clear;clc;
path1  = '/root_dir/wuguowei/code/FC_ontology/';
load([path1 filesep 'Datasets/Contrast_task.mat']);
sublist = Contrast_task.data;
func_data = load([path1 filesep 'HCPA_glasser_360_all_FC.mat']);%
repeat_lag = randperm(100);
behave_items = [5:30];%
behave_items_name = Contrast_task.name(5:end);
load('/root_dir/wuguowei/code/FC_ontology/FC_net_glasser/glasser_Cog_top_10.mat');%load('Cognition_pre_coef_200.mat');top_10_edge_cognition
mkdir('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_glasser');
cd('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_glasser');
save('repeat_lag.mat','repeat_lag');
%%
all_sub_func_rest = func_data.func_all(:,top_10_edge);
for n =1:numel(behave_items)
    behave_n = behave_items(n);
    if ~exist(['Real_top_10_R_noage_' behave_items_name{behave_n-4}],'dir')
        mkdir(['Real_top_10_R_noage_' behave_items_name{behave_n-4}]);
        cd(['Real_top_10_R_noage_' behave_items_name{behave_n-4}]);
        system('cp ../../codes/run_single_parcel.sh ./');
        system('cp ../../codes/run_ability_predict_100.m ./');
        save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
        cmd{1} = ['cd(''/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_glasser/Real_top_10_R_noage_'  behave_items_name{behave_n-4}  ''')'];
        cmd{2} = 'load(''use.mat'')';
        for p=1:100
            cmd{3} = ['[Results,rand_loction] = run_ability_predict_100(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),2:4)'];
            cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
            fid = fopen(['real_test_' num2str(p) '.m'],'wt');
            fprintf(fid,'%s\n',cmd{:});
            fclose(fid);
            %pause(10)
            %system(['sbatch -p node1 -c 1 '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
        end
            %pause(50)
    end
    disp(['Real_top_10_R_noage_' behave_items_name{behave_n-4}]);
    cd('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_glasser');
end


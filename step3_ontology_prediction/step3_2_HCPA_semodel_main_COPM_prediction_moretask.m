clear;clc;
path1  = '/ibmgpfs/cuizaixu_lab/wuguowei/code/FC_ontology/';
sublist = csvread([path1 filesep 'more_control_task.csv'],1);
func_data = csvread([path1 filesep 'HCPA_FC_all_data.csv'],1);%
sublist_cov = csvread([path1 filesep 'motor_combined.csv'],1);
sublist = [sublist_cov(:,1:4),sublist(:,1:18)];
load('more_control_task_name.mat');
load('top_10_edge_cognition.mat');%load('Cognition_pre_coef_200.mat');top_10_edge_cognition
%[sorted_value,index] = sort(coef_200,'descend');
%top_10_edge = index(1:1990);
repeat_lag = randperm(100);
behave_items = [5:21];%
behave_items_name = name;
%load('repeat_lag.mat');
cd('/ibmgpfs/cuizaixu_lab/wuguowei/code/FC_ontology/FC_net_Ontology_400');
save('repeat_lag.mat','repeat_lag');
%%
all_sub_func_rest = func_data(:,top_10_edge);
for n =1:numel(behave_items)
    behave_n = behave_items(n);
    if ~exist(['Real_top_10_R_noage_' behave_items_name{behave_n-4}],'dir')
        mkdir(['Real_top_10_R_noage_' behave_items_name{behave_n-4}]);
        cd(['Real_top_10_R_noage_' behave_items_name{behave_n-4}]);
        system('cp ../../codes/run_single_parcel.sh ./');
        system('cp ../../codes/run_ability_predict_100.m ./');
        save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
        cmd{1} = ['cd(''/ibmgpfs/cuizaixu_lab/wuguowei/code/FC_ontology/FC_net_Ontology_400/Real_top_10_R_noage_'  behave_items_name{behave_n-4}  ''')'];
        cmd{2} = 'load(''use.mat'')';
        for p=1:100
            cmd{3} = ['[Results,rand_loction] = run_ability_predict_100(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),2:4)'];
            cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
            fid = fopen(['real_test_' num2str(p) '.m'],'wt');
            fprintf(fid,'%s\n',cmd{:});
            fclose(fid);
            pause(5)
            system(['sbatch -p node1 -c 2 '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
        end
            pause(30)
    end

    cd('/ibmgpfs/cuizaixu_lab/wuguowei/code/FC_ontology/FC_net_Ontology_400');
end


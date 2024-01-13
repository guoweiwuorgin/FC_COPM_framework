clear;clc;
path1  = '/root_dir/wuguowei/code/FC_ontology/';
load([path1 filesep 'Datasets/Contrast_task.mat']);
load([path1 filesep 'Datasets/HCPA_Schaefer_400_FC.mat']);
load([path1 filesep 'Datasets/secondfactor_coef.mat']);
sublist = Contrast_task.data;
repeat_lag = randperm(100);
behave_items = [5:30];%
behave_items_name = Contrast_task.name(5:end);
%load('repeat_lag.mat');
mkdir('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_400_difftreshold');
cd('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_400_difftreshold');
save('repeat_lag.mat','repeat_lag');
%%
[sorted_value,index] = sort(abs(mean(ceof_400,1)),'descend');
top_10_edge = index(1:7980);
thresholds=[0.2,0.3,0.4,0.5];
for fold=1:4
    top_edge = index(1:(79800*thresholds(fold)));
    all_sub_func_rest = func_data(:,top_edge);
    for n =[3 5 8 9 12 13 20 23]
        behave_n = behave_items(n);
        data_dir_tar=['Real_top_10_R_noage_' behave_items_name{behave_n-4} '_thr_' num2str(thresholds(fold)*10)];
        if ~exist(data_dir_tar,'dir')
            mkdir(data_dir_tar);
            cd(data_dir_tar);
            system('cp ../../codes/run_single_parcel.sh ./');
            system('cp ../../codes/run_ability_predict_100.m ./');
            save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
            cmd{1} = ['cd(''/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_400_difftreshold/'  data_dir_tar ''')'];
            cmd{2} = 'load(''use.mat'')';
            for p=1:100
                cmd{3} = ['[Results,rand_loction] = run_ability_predict_100(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),2:4)'];
                cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
                fid = fopen(['real_test_' num2str(p) '.m'],'wt');
                fprintf(fid,'%s\n',cmd{:});
                fclose(fid);
                pause(5)
                system(['sbatch -p node1 -c 1 '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
            end
                pause(30)
        end
        cd('/root_dir/wuguowei/code/FC_ontology/FC_net_Ontology_400_difftreshold');
    end
end


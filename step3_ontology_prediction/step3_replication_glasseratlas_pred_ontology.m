clear;clc;
path1  = '/root_dir/wuguowei/code/FC_ontology/';
sublist = csvread([path1 filesep 'motor_combined.csv'],1);
func_data = load([path1 filesep 'HCPA_glasser_360_all_FC.mat']);%
load([path1 filesep 'motor_combined_data.mat']);
repeat_lag = randperm(100);
behave_items = [5:15];%
behave_items_name = name;
%load('repeat_lag.mat');
cd('/root_dir/wuguowei/code/FC_ontology/FC_net_glasser');
save('repeat_lag.mat','repeat_lag');
%%
all_sub_func_rest = func_data.func_all;
for n =11%:numel(behave_items)
    behave_n = behave_items(n);
    if ~exist(['Real_R_noage_runall_F' behave_items_name{behave_n}],'dir')
        mkdir(['Real_R_noage_runall_F' behave_items_name{behave_n}]);
        cd(['Real_R_noage_runall_F' behave_items_name{behave_n}]);
        system('cp ../../codes/run_single_parcel.sh ./');
        system('cp ../../codes/run_ability_predict_100.m ./');
        save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
        cmd{1} = ['cd(''/root_dir/wuguowei/code/FC_ontology/FC_net_glasser/Real_R_noage_runall_F'  behave_items_name{behave_n}  ''')'];
        cmd{2} = 'load(''use.mat'')';
        for p=1:100
            cmd{3} = ['[Results,rand_loction] = run_ability_predict_100(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),2:4)'];
            cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
            fid = fopen(['real_test_' num2str(p) '.m'],'wt');
            fprintf(fid,'%s\n',cmd{:});
            fclose(fid);
            pause(10)
            system(['sbatch -p node1 '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
        end
    end
    pause(150)
    cd('/root_dir/wuguowei/code/FC_ontology/FC_net_glasser');
end


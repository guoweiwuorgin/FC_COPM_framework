clear;clc;
path1  = '/root_dir/wuguowei/code/FC_ontology/';
sublist = csvread([path1 filesep 'motor_combined.csv'],1);
func_data = csvread([path1 filesep 'HCPA_FC_all_data.csv'],1);%
load([path1 filesep 'motor_combined_data.mat']);
repeat_lag = randperm(200);
behave_items = [8:15];%
behave_items_name = name;
%load('repeat_lag.mat');
cd('/root_dir/wuguowei/code/FC_ontology/FC_net_400');
save('repeat_lag_permu.mat','repeat_lag');
%%
all_sub_func_rest = func_data;
for n =1:numel(behave_items)
    behave_n = behave_items(n);
    if ~exist(['Permu_R_noage_' behave_items_name{behave_n}],'dir')
        mkdir(['Permu_R_noage_' behave_items_name{behave_n}]);
        cd(['Permu_R_noage_' behave_items_name{behave_n}]);
        system('cp ../../codes/run_single_parcel.sh ./');
        system('cp ../../codes/run_ability_predict_permu.m ./');
        save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
        cmd{1} = ['cd(''/root_dir/wuguowei/code/FC_ontology/FC_net_400/Permu_R_noage_'  behave_items_name{behave_n}  ''')'];
        cmd{2} = 'load(''use.mat'')';
        for p=1:200
            cmd{3} = ['[Results,rand_loction] = run_ability_predict_permu(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),2:4)'];
            cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
            fid = fopen(['real_test_' num2str(p) '.m'],'wt');
            fprintf(fid,'%s\n',cmd{:});
            fclose(fid);
            pause(10)
            system(['sbatch '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
        end
    end
    pause(50)
    cd('/root_dir/wuguowei/code/FC_ontology/FC_net_400');
end


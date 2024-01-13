clear;clc;
path1  = '/root_dir/wuguowei/code/FC_ontology/';
load([path1 filesep 'HCPD_pred_test/HCPD_all_behav.mat']);
HCPD_func_dir = [path1 filesep 'FC_net_HCPD'];
repeat_lag = randperm(100);
cd([path1 filesep 'HCPD_pred_test_all/']);
save('repeat_lag.mat','repeat_lag');
%%
for n =[2:5,7:9]
    sublist = all_behav{n,3};
    behave_name = all_behav{n,1};
    behave_subin = all_behav{n,2};
    func_data = zeros(length(behave_subin),79800);
    for sub =1:length(behave_subin)
        func_dir = [HCPD_func_dir filesep 'sub-' behave_subin(sub,:) filesep 'schaefer_400.mat'];
        load(func_dir);
        func_data(sub,:) = rest_conn;
    end
    all_sub_func_rest = func_data;
    [a,b] = find(isnan(sublist));
    if length(a)>0
        sublist(a,:)=[];
        all_sub_func_rest(a,:)=[];
    end
    behave_n=4;
    if ~exist(['Real_R_noage_' behave_name],'dir')
        mkdir(['Real_R_noage_' behave_name]);
        cd(['Real_R_noage_' behave_name]);
        system('cp ../../codes/run_single_parcel.sh ./');
        system('cp ../../codes/run_ability_predict_100.m ./');
        save('use.mat','behave_n','sublist','all_sub_func_rest','repeat_lag');
        cmd{1} = ['cd(''/root_dir/wuguowei/code/FC_ontology/HCPD_pred_test_all/Real_R_noage_'  behave_name  ''')'];
        cmd{2} = 'load(''use.mat'')';
        for p=1:100
            cmd{3} = ['[Results,rand_loction] = run_ability_predict_100(behave_n,sublist,all_sub_func_rest,repeat_lag(' num2str(p) '),1:3)'];
            cmd{4} = ['save(''data_re_' num2str(p) '.mat'',''Results'',''rand_loction'')'];
            fid = fopen(['real_test_' num2str(p) '.m'],'wt');
            fprintf(fid,'%s\n',cmd{:});
            fclose(fid);
            pause(10)
           system(['sbatch '  'run_single_parcel.sh '  'real_test_' num2str(p)]);
        end
    end
    pause(10)
    cd('/root_dir/wuguowei/code/FC_ontology/HCPD_pred_test_all');
end


clear;clc;
all_items = dir('Real_*');
iter_n = 1;
for n = 1:numel(all_items)
    folder_name = all_items(n).name;
    cd(folder_name);
    real_R{1,iter_n} = strrep(folder_name,'Real_R_noage_','');
    ceof_all{1,iter_n} = strrep(folder_name,'Real_R_noage_','');
    all_per = dir('data_re*.mat');
    for p=1:numel(all_per)
        load(all_per(p).name);
        R1 =  ~isnan(Results.y_test_all(:,1));
        R2 =  ~isnan(Results.y_test_all(:,2));
        R(1)=corr(Results.y_test_all(R1,1),Results.yhat(R1,1));
        R(2)=corr(Results.y_test_all(R2,2),Results.yhat(R2,2));
%         ceof_all{p+1,iter_n} = mean(Results.coef,2);
%         ceof_400(p,:) = mean(Results.coef,2);
        real_R{p+1,iter_n} = mean(R);
%         real_contribution = mean(Results.coef,2);
%         all_fold_data(n,:,p) = real_contribution;
    end
    cd('..');
    iter_n = iter_n+1;
end
[sorted_value,index] = sort(ceof_400,'descend');
top_10_edge = index(1:7980);

save('secondfactor_coef.mat','ceof_400','index','sorted_value','top_10_edge');


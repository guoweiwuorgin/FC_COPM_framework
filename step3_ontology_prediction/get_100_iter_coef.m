function all_fold_data = get_100_iter_coef(dir_folder)
    cd(dir_folder);
    n=1
    iter_n=1
    all_per = dir('data_re*.mat');
    for p=1:numel(all_per)
        load(all_per(p).name);
        R1 =  ~isnan(Results.y_test_all(:,1));
        R2 =  ~isnan(Results.y_test_all(:,2));
        R(1)=corr(Results.y_test_all(R1,1),Results.yhat(R1,1));
        R(2)=corr(Results.y_test_all(R2,2),Results.yhat(R2,2));
        ceof_all{p+1,iter_n} = mean(Results.coef,2);
        real_R{p+1,iter_n} = mean(R);
        real_contribution = mean(Results.coef,2);
        all_fold_data(n,:,p) = real_contribution;
    end
end
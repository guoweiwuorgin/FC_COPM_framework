function data_b_label = compute_mat_base_label(data,label,label_vertex)
    data_b_label = zeros(length(label),size(data,2));
    for l = 1:length(label)
        label_all_ver = find(label_vertex==label(l));
        if label(l)==0
            data_b_label(l,:) = zeros(1,size(data,2));
        else
            data_b_label(l,:) = mean(data(label_all_ver,:));
        end
    end
end
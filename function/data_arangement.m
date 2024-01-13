clear;clc;
%% data arrangement
path1  = '/root_dir/wuguowei/code/FC_ontology/';
sublist = csvread([path1 filesep 'motor_combined.csv'],1);
func_data = csvread([path1 filesep 'HCPA_FC_all_data.csv'],1);%
load([path1 filesep 'motor_combined_data.mat']);

Ontology.data = sublist(:,[1:4,9:15]);
Ontology.name = name([1:4,9:15]);
save('./Datasets/Ontology_behaviors.mat','Ontology');
save('./Datasets/Schaefer_400_FC.mat','func_data');
% get motor behavoirs 
motordata = sublist(:,8);
motordatanames = name(8);
% emotion data
sublist_all = csvread([path1 filesep 'Ontology_use_data.csv'],1);
load('Test_use_ontology.mat')
emotiondata = sublist_all(:,10:12);
emotionname = name(10:12);

othercogtask = sublist_all(:,17:20);
othercogtaskname = name(17:20);
% all control tasks
sublist = csvread([path1 filesep 'more_control_task.csv'],1);
sublist_cov = csvread([path1 filesep 'motor_combined.csv'],1);
sublist = [sublist_cov(:,1:4),sublist(:,1:18),motordata,emotiondata,othercogtask];
load('more_control_task_name.mat');
load('top_10_edge_cognition.mat');
behave_items_name = [Ontology.name(1:4);name;motordatanames;emotionname;othercogtaskname];
Ontologytop10_schaefer_func_rest = func_data(:,top_10_edge);
Contrast_task.data = sublist;
Contrast_task.name = behave_items_name;
save('./Datasets/Contrast_task.mat','Contrast_task');
save('./Datasets/Ontology_egde_top10.mat','top_10_edge');




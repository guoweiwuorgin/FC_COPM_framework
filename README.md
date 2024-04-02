# FC_ontology_framework
This repo is using for the FC ontology project replication  
1. Data Avaliable
   The raw data of  T1-weight MRI, functional MRI, diffusion-weighted MRI, and behaviral data of the HCP-YA and HCP-D datasets are available at https://db.humanconnectome.org/. Gene expression data can be downloaded from the AHBA (http://human.brain-map.org) and processed using the abagen toolbox (https://github.com/rmarkello/abagen).
   The motion information of subject, the ontology score, the behaviral score of HCPA and HCPD dataset are aviliable in Datasets folder.
   
3. Code use instruction
   The ontology construction code for all three different model in the paper can be found in folder: step1_ontology model.
   
   The functional connectivity extraction code can be found in folder: step2_FC_construction.

   The ridge regression model for predicting the cognitive ontology score, the Cognitive Ontology based Prediction Model (COPM) model, full-FC feature and the Ridge-regression based Connectome Predictive Model all are avilable in the folder: step3_ontology_prediction

   step_4_0_secondorder_model_plot_and_networkpreporty_corr.R was used to plot  Figure 2 and Figure 3.

   step_4_1_bifactor_plot.R was used to plot Supplementary Figure 1.

   step_4_2_More_task_interpretation.R was used to plot figure 3 and figure 4.

   step_4_1_HCPD_plot.R was used to plot figure5.

   step_final_diffthreshold_and_Glasser_atlas.R was used to replicate our main analysis on different threhold value and parcel resolution.

   
   

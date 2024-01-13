setwd("D:/FC_ontology_project")
library(bruceR)
all_behave_data <- import("HCPA_all_behav_data.xlsx")
all_behave_data %>% filter(`3T_Full_MR_Compl`==TRUE,Handedness>0) -> MRI_have_sub

cognition_data <- import("D:/FC_ontology_project/HCPA_behaviour.csv")
MRI_sub <-  import("subname.txt")

MRI_have_sub %>% select(Subject,Gender,Age_in_Yrs,Handedness,Race,SSAGA_Educ,BMI)  %>% inner_join(cognition_data,by=c("Subject","Gender")) %>% 
  select(Subject)-> func_subjects


MRI_sub$subname <- as.numeric(str_remove(MRI_sub$subname,"sub-"))
will_used_data <- data.frame(Subject=intersect(MRI_sub$subname,func_subjects$Subject))

library(EFAtools)

more_cognition <- c("PicSeq_Unadj","CardSort_Unadj","Flanker_Unadj","PMAT24_A_CR",
                    "PMAT24_A_SI","PMAT24_A_RTCR","ReadEng_Unadj","PicVocab_Unadj",
                    "ProcSpeed_Unadj","VSPLOT_CRTE","SCPT_TPRT","SCPT_SPEC","SCPT_SEN","SCPT_LRNR",
                    "VSPLOT_TC","VSPLOT_CRTE","VSPLOT_OFF",
                    "IWRD_RTC","ListSort_Unadj","Language_Task_Story_Avg_Difficulty_Level",
                    "Language_Task_Story_Acc","Language_Task_Story_Median_RT",
                    "Relational_Task_Match_Median_RT","Relational_Task_Rel_Median_RT",
                    "WM_Task_2bk_Acc","WM_Task_2bk_Median_RT",
                    "Dexterity_Unadj","Strength_Unadj",'GaitSpeed_Comp','Endurance_Unadj',"DDisc_AUC_40K",
                    "PercStress_Unadj","Loneliness_Unadj","PercReject_Unadj","InstruSupp_Unadj","LifeSatisf_Unadj",
                    "IWRD_TOT"
                    )

emotion_items <- c("AngAffect_Unadj","FearAffect_Unadj","Sadness_Unadj")


MRI_have_sub %>% select(Subject,Gender,Age_in_Yrs,Handedness,Race,SSAGA_Educ,BMI)  %>% inner_join(cognition_data,by=c("Subject","Gender")) %>% 
  select(Subject,Gender,Age_in_Yrs,Race,SSAGA_Educ,all_of(more_cognition),emotion_items)-> all_func_data_items

mean_motion_data<- import("project_code/motion_info_func.xlsx")
mean_motion_data$Subject <- as.numeric(str_remove(mean_motion_data$Subject,"sub-"))
all_func_data_items <- all_func_data_items %>% inner_join(mean_motion_data)%>% filter(mean_motion<0.2)%>%na.omit()

all_func_data_items %>% select(Subject,Gender,Age_in_Yrs,SSAGA_Educ,Race) %>% group_by(Gender,Race)%>%
  summarise(mean_age=mean(Age_in_Yrs),n_gender=n(),min_age=min(Age_in_Yrs),max_age=max(Age_in_Yrs),mean_edu=mean(SSAGA_Educ))
i<-0
for( n in unique(all_func_data_items$Race)){
  all_func_data_items$Race[all_func_data_items$Race==n]<-i
  i<- i+1}
all_func_data_items$Race <- as.numeric(all_func_data_items$Race)
all_func_data_items %>% export("all_func_data_items.xlsx")
all_func_data_items %>% export("all_func_data_items.csv")

all_func_data_items %>% 
  select(more_cognition,emotion_items)-> all_func_data_items_new

all_func_data_items_new %>% select(Language_Task_Story_Avg_Difficulty_Level,VSPLOT_CRTE,IWRD_RTC,SCPT_TPRT)%>% export("Contrast_task_data.mat")
all_func_data_items_new %>% select(Language_Task_Story_Avg_Difficulty_Level,VSPLOT_CRTE,IWRD_RTC,SCPT_TPRT)-> contrast_task
all_func_data_items_new %>% na.omit() %>% apply(2,scale)%>% as.data.frame() -> z_score_data

correlations <- cor(z_score_data)
EFAtools::BARTLETT(correlations,N=nrow(all_func_data_items_new))
EFAtools::KMO(correlations)
PARALLEL(correlations, eigen_type = "EFA", n_datasets = 1000,N = nrow(all_func_data_items_new))

EFA_AV <- EFA_AVERAGE(as.matrix(correlations), n_factors = 9, N = nrow(all_func_data_items_new),
                      method = c("PAF", "ML", "ULS"), rotation = "promax",varimax_type = c("svd"),start_method = c("psych"),
                      show_progress = FALSE)
library(lavaan)
library(foreign) 
library(lavaan)
m2a  <- ' EF1 =~ CardSort_Unadj+Flanker_Unadj + ProcSpeed_Unadj
          EF2 =~ WM_Task_2bk_Acc + ListSort_Unadj+PicSeq_Unadj
          FI=~PMAT24_A_CR+PMAT24_A_SI+PMAT24_A_RTCR
          Relation =~ Relational_Task_Match_Median_RT+Relational_Task_Rel_Median_RT
          LAN =~ PicVocab_Unadj+ReadEng_Unadj
          Cog =~ EF1+EF2+Relation+LAN+FI '


m2a_bi  <- 'EF1 =~ CardSort_Unadj+Flanker_Unadj + ProcSpeed_Unadj
          EF2 =~ WM_Task_2bk_Acc + ListSort_Unadj+PicSeq_Unadj
          FI=~PMAT24_A_CR+PMAT24_A_SI+PMAT24_A_RTCR
          Relation =~ Relational_Task_Match_Median_RT+Relational_Task_Rel_Median_RT
          LAN =~ PicVocab_Unadj+ReadEng_Unadj
          Cog =~ CardSort_Unadj+Flanker_Unadj+ProcSpeed_Unadj+WM_Task_2bk_Acc+ListSort_Unadj+PicSeq_Unadj+PMAT24_A_CR+PMAT24_A_SI+PMAT24_A_RTCR+Relational_Task_Match_Median_RT+Relational_Task_Rel_Median_RT+PicVocab_Unadj+ReadEng_Unadj'

iter<-1
z_score_data_residual <- cbind(z_score_data,all_func_data_items[,2:3])
for(variable in names(z_score_data_residual[,1:32])){
  lm_reg_formula <- as.formula(paste0(variable,'~Gender+Age_in_Yrs'))
  lm_reg_age_sex_IQ <- lm(lm_reg_formula,z_score_data_residual)
  if(iter==1){
    reg_cov_time3 <- data.frame(lm_reg_age_sex_IQ$residuals)
  }else{
    reg_cov_time3 <- cbind(reg_cov_time3,data.frame(lm_reg_age_sex_IQ$residuals))
  }
  iter <- iter+1
}
names(reg_cov_time3)<- names(z_score_data_residual[,1:32])
z_score_data <- select(reg_cov_time3,CardSort_Unadj,Flanker_Unadj,ProcSpeed_Unadj,
                        WM_Task_2bk_Acc,ListSort_Unadj,PicSeq_Unadj,PMAT24_A_CR,
                        PMAT24_A_SI,PMAT24_A_RTCR,Relational_Task_Match_Median_RT,
                        Relational_Task_Rel_Median_RT,PicVocab_Unadj,ReadEng_Unadj)
fit_behav2 <- lavaan::cfa(m2a_bi, data = z_score_data, std.lv=TRUE,estimator='ML', orthogonal = TRUE , mimic =c("MPlus"),check.gradient = FALSE)
fitMeasures(fit_behav2, c('chisq', 'df', 'pvalue', 'cfi', 'rmsea', 'srmr', 'AIC'))
fit_behav3 <- lavaan::cfa(m2a, data = z_score_data, std.lv=TRUE, orthogonal = TRUE , mimic =c("MPlus"),check.gradient = FALSE)

bruceR::CFA(z_score_data,model =m2a)

predicted_latent_var1 <- lavPredict(fit_behav2)
predicted_latent_var2 <- lavPredict(fit_behav3)

motor_combined <- data.frame(all_func_data_items[,c(1:3,40,29:32)],predicted_latent_var1)
 
ontology_included_function <- c("CardSort_Unadj",'Flanker_Unadj','ProcSpeed_Unadj','WM_Task_2bk_Acc','ListSort_Unadj',
                                'PicSeq_Unadj','PMAT24_A_CR','PMAT24_A_SI','PMAT24_A_RTCR','AngAffect_Unadj','FearAffect_Unadj',
                                'Sadness_Unadj','Relational_Task_Match_Median_RT','Relational_Task_Rel_Median_RT','PicVocab_Unadj',
                                'ReadEng_Unadj')
contrast_task_function <- c('Language_Task_Story_Avg_Difficulty_Level','VSPLOT_CRTE','IWRD_RTC','SCPT_TPRT','Endurance_Unadj')
Ontology_use_data <- all_func_data_items %>% select(c(ontology_included_function,contrast_task_function))
export(Ontology_use_data,"Ontology_use_data.csv")
export(data.frame(name=names(Ontology_use_data)),"Test_use_ontology.mat")

motor_combined$Gender[motor_combined$Gender=="M"]<-0
motor_combined$Gender[motor_combined$Gender=="F"]<-1
motor_combined$Gender<-as.numeric(motor_combined$Gender)
export(motor_combined,"motor_combined.csv")


sum_fit <- summary(fit_behav2, fit.measures=TRUE, standardized=TRUE)
library(sem)

library(semPlot)
pdf("psycho_ontology.pdf", height=60, width=70)
semPaths(fit_behav3, what = "equality",
         whatLabels = "std",             # std = factor loading
         groups = "latents",             # 依latent 分組上色
         pastel = TRUE,                  # 柔和色
         mar = c(3, 3, 3, 3),            # 邊界
         intercepts = FALSE,              # 殘差顯示
         edge.label.cex = 0.6,             # label 的字體大小
         sizeMan = 6,                    # measurement variable 顯示大小
         sizeLat = 8,                   # latent variable 顯示大小
         exoCov = FALSE,                 # 相互影響的線不要
         optimizeLatRes = TRUE,          # 讓latent 的殘差長整齊一點
         style="ram",                 # 讓殘差格式跟LISREL一樣
         layout="tree2",
         rotation = 1,                   # 橫向
         structural=FALSE)
dev.off()


pdf("psycho_ontology_bimodel.pdf", height=60, width=70)
semPaths(fit_behav2, what = "equality",
         whatLabels = "std",             # std = factor loading
         groups = "latents",             # 依latent 分組上色
         pastel = TRUE,                  # 柔和色
         mar = c(3, 3, 3, 3),            # 邊界
         intercepts = FALSE,              # 殘差顯示
         edge.label.cex = 0.6,             # label 的字體大小
         sizeMan = 6,                    # measurement variable 顯示大小
         sizeLat = 8,                   # latent variable 顯示大小
         exoCov = FALSE,                 # 相互影響的線不要
         optimizeLatRes = TRUE,          # 讓latent 的殘差長整齊一點
         style="ram",                 # 讓殘差格式跟LISREL一樣
         layout="tree2",
         rotation = 1,                   # 橫向
         structural=FALSE,
         bifactor="Cog")
dev.off()

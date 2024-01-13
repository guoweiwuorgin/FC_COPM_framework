setwd("D:/FC_ontology_project/New_results")
library(bruceR)
HCP_full_model <- import("HCPA_all_FC_400.xlsx")
HCP_CPM_model <- import("HCPA_CPM.xlsx")
HCP_ontology_model <- import("bimodel_ontology_prediction.xlsx")
tasknames <- import("more_task.csv")
names(HCP_full_model) %>% intersect(names(HCP_CPM_model))%>% 
  intersect(names(HCP_ontology_model)) %>% 
  intersect(tasknames$`names(useful_data)`) -> task_use_more
task_use_more[task_use_more!="IWRD_TOT"] ->task_use_more
HCP_full_model <- HCP_full_model %>% select(all_of(task_use_more))
HCP_CPM_model <- HCP_CPM_model %>% select(all_of(task_use_more))
HCP_ontology_model <- HCP_ontology_model %>% select(all_of(task_use_more))
ontology_score_pred_permu <- import("Permutation_R_ontology_score.xlsx")%>% mutate(Type="Permuataion",pred_folder=c(1:200))



HCP_full_model <- HCP_full_model %>%mutate(pred_type = "Full_FC",pred_folder=c(1:100))
HCP_CPM_model <- HCP_CPM_model %>%mutate(pred_type = "rCPM",pred_folder=c(1:100))
HCP_ontology_model <- HCP_ontology_model %>%mutate(pred_type = "Ontology",pred_folder=c(1:100))

all_data <- rbind(HCP_full_model,HCP_CPM_model) %>% rbind(HCP_ontology_model)%>% 
  pivot_longer(cols = 1:25,names_to = "Behavior",values_to = "Predict_R")
all_data$Behavior <- str_remove(all_data$Behavior ,"_Unadj")
re_COPM_F_w <- c()
re_COPM_F_p <- c()
re_COPM_rCPM_w <- c()
re_COPM_rCPM_p <- c()
name_behave <- c()
i<-1
for (name_var in unique(all_data$Behavior)) {
  tmp_data <- filter(all_data,Behavior==name_var)
  tmp_data_COPM_F <- dplyr::filter(tmp_data,pred_type!="rCPM")
  re_COPM_F <- wilcox.test(Predict_R~pred_type,tmp_data_COPM_F,paired=T)
  re_COPM_F_w[i]<-re_COPM_F$statistic
  re_COPM_F_p[i]<-re_COPM_F$p.value
  tmp_data_COPM_rCPM <- filter(tmp_data,pred_type!="Full_FC")
  re_COPM_rCPM <- wilcox.test(Predict_R~pred_type,tmp_data_COPM_rCPM,paired=T)
  re_COPM_rCPM_w[i]<-re_COPM_rCPM$statistic
  re_COPM_rCPM_p[i]<-re_COPM_rCPM$p.value
  name_behave[i] <- name_var
  i <- i+1
}

paired_re <- data.frame(name_behave,re_COPM_F_w,re_COPM_F_p,re_COPM_rCPM_w,re_COPM_rCPM_p)

paired_re <- mutate(paired_re,adj_re_COPM_F_p=p.adjust(re_COPM_F_p,method = "fdr")<0.05,adj_re_COPM_F_p1=p.adjust(re_COPM_F_p,method = "fdr")<0.01,
                    adj_re_COPM_F_p2=p.adjust(re_COPM_F_p,method = "fdr")<0.001,adj_re_COPM_rCPM_p=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.05,
                    adj_re_COPM_rCPM_p1=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.01,adj_re_COPM_rCPM_p2=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.001)



all_data <- all_data %>% 
  mutate(Ontology_type = rep(c("EMO","EMO","VEM","LAN","SATT","EMO","SpO","Motor","SRL","Motor","Motor","SR",
                               "LAN","LAN","LifeSatisf","SR","SR","SS","SATT","SATT","SATT","Motor",
                               "SpO","SpO","EF"),300))
cognitive_ontology_relation <- import("More_task_cognitive_ontology_relation.xlsx")
names(cognitive_ontology_relation)[1]<-"Behavior"
all_data$Behavior <- factor(all_data$Behavior,levels = c("WM_Task_2bk_Median_RT","Language_Task_Story_Avg_Difficulty_Level",
                                                         "Language_Task_Story_Acc","Language_Task_Story_Median_RT","DDisc_AUC_40K",
                                                         "SCPT_LRNR","SCPT_SEN","SCPT_SPEC","SCPT_TPRT","IWRD_RTC","VSPLOT_CRTE","VSPLOT_OFF",
                                                         "VSPLOT_TC","Dexterity","Endurance","GaitSpeed_Comp","Strength","AngAffect",
                                                         "FearAffect","Sadness","Loneliness","PercReject","PercStress",
                                                         "InstruSupp","LifeSatisf"))

all_data$pred_type <- factor(all_data$pred_type,levels = c("Ontology","Full_FC","rCPM"),ordered = TRUE)
all_data$Ontology_type <-  factor(all_data$Ontology_type,levels = c("EF","LAN","SRL","SATT","VEM","SpO",
                                                                "Motor","EMO","SR","SS"),ordered = TRUE)

all_data <- all_data[order(all_data$Ontology_type),]

cognitive_ontology_relation$abs_R <- round(cognitive_ontology_relation$abs_R,2)

color_sets2 <- c("#9b1ed1","#ca7aec","#e8c4f7")
comparelist <- list(c("Full_FC","Ontology"),c("rCPM","Ontology"))
all_data %>% 
  ggplot(aes(x=pred_type,y=Predict_R,fill=pred_type)) + 
  geom_flat_violin(position=position_nudge(x=0, y=0), alpha=0.6, width=1.1,color='white') +
  stat_summary(
    alpha = 1,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_nudge(0.06),width = 0.05
  )+facet_wrap(~Behavior,nrow = 5)+scale_y_continuous(limits = c(-0.1,0.6))+
  theme_classic(base_size = 20)+geom_hline(yintercept = 0,linetype=2)+
  theme(strip.background = element_blank())+
  ggsignif::geom_signif(comparisons = list(c("Full_FC","Ontology")),test = wilcox.test,
                        map_signif_level = TRUE,
                        textsize = 4,y_position = 0.55 )+
  ggsignif::geom_signif(comparisons = list(c("rCPM","Ontology")),test = wilcox.test,
                        map_signif_level = TRUE,
                        textsize = 4,y_position = 0.5 )+
  scale_fill_manual(values = color_sets2)->pred_all



ggsave(filename = "moretask_pred_bimodel.pdf",pred_all,width = 500, 
       height = 600, dpi = 300, units = "mm", device='pdf')

all_data  %>% group_by(pred_type,Behavior,Ontology_type) %>% 
  summarise(mean_R=mean(Predict_R),sd=sd(Predict_R))%>% 
  pivot_wider(id_cols = c(Ontology_type,Behavior),names_from ="pred_type", 
              values_from = "mean_R") -> describe_R
export(describe_R,"bimodel_Mean_pred_R.xlsx")
all_data  %>% group_by(pred_type,Behavior,Ontology_type) %>% 
  summarise(mean_R=mean(Predict_R),sd=sd(Predict_R))%>% 
  pivot_wider(id_cols = c(Ontology_type,Behavior),names_from ="pred_type", 
              values_from = "sd") -> describe_SD_R
export(describe_SD_R,"bimodel_SD_pred_R.xlsx")




setwd("D:/FC_ontology_project/New_results")
library(bruceR)
HCPD_CPM <- import("Random_HCPD.xlsx")%>%mutate(pred_type = "rCPM",pred_folder=c(1:100))
HCPD_ALLFC <- import("HCPD_ALLFC.xlsx")%>%mutate(pred_type = "Full_FC",pred_folder=c(1:100))%>% select(names(HCPD_CPM))
HCPD_Cogpred <- import("HCPD_cognition_pred.xlsx")%>%mutate(pred_type = "COPM",pred_folder=c(1:100))%>% select(names(HCPD_CPM))


all_HCPD_data <- rbind(HCPD_CPM,HCPD_ALLFC,HCPD_Cogpred) %>% 
  pivot_longer(cols = 1:8,names_to = "Behavioral",values_to = "Prediction_R") %>% 
  mutate(Task_type = rep(c("EF","Emotion","EF","EF","EF","LAN","LAN","Motor"),300))

all_HCPD_data$Behavioral <- factor(all_HCPD_data$Behavioral,levels = c("dccs_score","flanker_score","picseq_score",
                                                                       "lswmt_score","pic_volcabulary_score","reading_score",
                                                                       "emo_score","motor"),ordered = TRUE)
task_id <- rep(c(1,0,1,1,1,1,1,0),100)
library(rstatix)
all_HCPD_data %>% filter(pred_type=='COPM')%>%mutate(task_id=task_id) %>% group_by(pred_folder,task_id)%>%
  summarise(mean_Prediction_R=mean(Prediction_R)) %>% ungroup( )->High_low_comp

wilcox.test(mean_Prediction_R~task_id,High_low_comp,paired=T)
wilcox_effsize(High_low_comp,mean_Prediction_R~task_id,paired=T)

re_COPM_F_w <- c()
re_COPM_F_p <- c()
re_COPM_rCPM_w <- c()
re_COPM_rCPM_p <- c()
name_behave <- c()
i<-1
i<-1
for (name_var in unique(all_HCPD_data$Behavioral)) {
  tmp_data <- filter(all_HCPD_data,Behavioral==name_var)
  tmp_data_COPM_F <- dplyr::filter(tmp_data,pred_type!="rCPM")
  re_COPM_F <- wilcox.test(Prediction_R~pred_type,tmp_data_COPM_F,paired=T)
  re_COPM_r <- wilcox_effsize(tmp_data_COPM_F,Prediction_R~pred_type,paired=T)
  re_COPM_F_w[i]<-re_COPM_r$effsize[1]
  re_COPM_F_p[i]<-re_COPM_F$p.value
  tmp_data_COPM_rCPM <- filter(tmp_data,pred_type!="Full_FC")
  re_COPM_rCPM <- wilcox.test(Prediction_R~pred_type,tmp_data_COPM_rCPM,paired=T)
  re_COPM_rCPMr <- wilcox_effsize(tmp_data_COPM_rCPM,Prediction_R~pred_type,paired=T)
  re_COPM_rCPM_w[i]<-re_COPM_rCPMr$effsize[1]
  re_COPM_rCPM_p[i]<-re_COPM_rCPM$p.value
  name_behave[i] <- name_var
  i <- i+1
}

paired_re <- data.frame(name_behave,re_COPM_F_w,re_COPM_F_p,re_COPM_rCPM_w,re_COPM_rCPM_p)

paired_re <- mutate(paired_re,adj_re_COPM_F_p=p.adjust(re_COPM_F_p,method = "fdr"),adj_re_COPM_F_p1=p.adjust(re_COPM_F_p,method = "fdr")<0.01,
                    adj_re_COPM_F_p2=p.adjust(re_COPM_F_p,method = "fdr")<0.001,adj_re_COPM_rCPM_p=p.adjust(re_COPM_rCPM_p,method = "fdr"),
                    adj_re_COPM_rCPM_p1=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.01,adj_re_COPM_rCPM_p2=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.001)
export(paired_re,"Ontology_paired_re_HCPD.xlsx")


re_COPM_F_w <- c()
re_COPM_F_p <- c()
re_COPM_rCPM_w <- c()
re_COPM_rCPM_p <- c()
name_behave <- c()
i<-1
for (name_var in unique(all_HCPD_data$Behavioral)) {
  tmp_data <- filter(all_HCPD_data,Behavioral==name_var)
  tmp_data_COPM_F <- dplyr::filter(tmp_data,pred_type!="rCPM")
  re_COPM_F <- wilcox.test(Prediction_R~pred_type,tmp_data_COPM_F,paired=T)
  re_COPM_F_w[i]<-re_COPM_F$statistic
  re_COPM_F_p[i]<-re_COPM_F$p.value
  tmp_data_COPM_rCPM <- filter(tmp_data,pred_type!="Full_FC")
  re_COPM_rCPM <- wilcox.test(Prediction_R~pred_type,tmp_data_COPM_rCPM,paired=T)
  re_COPM_rCPM_w[i]<-re_COPM_rCPM$statistic
  re_COPM_rCPM_p[i]<-re_COPM_rCPM$p.value
  name_behave[i] <- name_var
  i <- i+1
}

paired_re <- data.frame(name_behave,re_COPM_F_w,re_COPM_F_p,re_COPM_rCPM_w,re_COPM_rCPM_p)

paired_re <- mutate(paired_re,adj_re_COPM_F_p=p.adjust(re_COPM_F_p,method = "fdr")<0.05,adj_re_COPM_F_p1=p.adjust(re_COPM_F_p,method = "fdr")<0.01,
                    adj_re_COPM_F_p2=p.adjust(re_COPM_F_p,method = "fdr")<0.001,adj_re_COPM_rCPM_p=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.05,
                    adj_re_COPM_rCPM_p1=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.01,adj_re_COPM_rCPM_p2=p.adjust(re_COPM_rCPM_p,method = "fdr")<0.001)





library(ggsignif)
HCPD_stats <- data.frame(cpm_Ontology,ALLFC_Ontology,F_tabel,behav_name)
library(PupillometryR)
color_HCPD <- c("#e295a5","#54bd95","#98a9e1","#87E293")
library(ggsignif)
all_HCPD_data$pred_type <- factor(all_HCPD_data$pred_type,levels = c("COPM","Full_FC","rCPM"))
color_sets2 <- c("#9b1ed1","#ca7aec","#e8c4f7")
all_HCPD_data %>% 
  ggplot(aes(x=pred_type,y=Prediction_R,fill=pred_type)) + 
  geom_flat_violin(position=position_nudge(x=0, y=0), alpha=0.6, width=1.1,color='white') +
  stat_summary(
    alpha = 1,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_nudge(0.06),width = 0.05
  )+facet_wrap(~Behavioral,nrow = 2)+scale_y_continuous(limits = c(-0.2,0.45))+
  theme_classic(base_size = 40)+geom_hline(yintercept = 0,linetype=2)+
  theme(strip.background = element_blank())+
  ggsignif::geom_signif(comparisons = list(c("Full_FC","COPM")),test = wilcox.test,
                        map_signif_level = TRUE,
                        textsize = 4,y_position = 0.4 )+
  ggsignif::geom_signif(comparisons = list(c("rCPM","COPM")),test = wilcox.test,
                        map_signif_level = TRUE,
                        textsize = 4,y_position = 0.35 )+
  scale_fill_manual(values = color_sets2)->HCPD_ontology_pred_sig

ggsave(filename = "moretask_pred_semodel_HCPD.pdf",HCPD_ontology_pred_sig,width = 400, 
       height = 300, dpi = 300, units = "mm", device='pdf')


all_HCPD_data %>% group_by(Behavioral,pred_type )%>%summarise(mean_R=mean(Prediction_R),
                                                  sd_R =sd(Prediction_R))%>%
  export("HCPD_pred_RE.csv")









ggsave(filename = paste0("HCPD_ontology_pred_all_sig.png"),HCPD_ontology_pred_sig,width = 600, 
       height = 420, dpi = 300, units = "mm", device='png')
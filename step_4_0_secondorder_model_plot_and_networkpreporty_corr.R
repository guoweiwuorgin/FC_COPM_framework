setwd("E:/FC_ontology_project/New_results")
library(bruceR)
library(PupillometryR)
ontology_score_pred <- import("Ontology_Score_prediction.xlsx") %>% mutate(Type="Real",pred_folder=c(1:100))
names(ontology_score_pred)<-str_remove(names(ontology_score_pred),'runall_F')
ontology_score_pred_permu <- import("Permutation_R_ontology_score.xlsx")%>% mutate(Type="Permuataion",pred_folder=c(1:200))
ontology_score_pred %>% rbind(ontology_score_pred_permu) -> ontology_score_pred_all
ontology_score_pred_all %>% pivot_longer(cols = c(1:8),names_to = "Ontology_categary",values_to = "Predict_R")->long_ontology_score_pred_all
long_ontology_score_pred_all %>%mutate(Type_comb=ifelse(Type=="Real",paste0(Type,Ontology_categary),Type))->tmp_data

for (type in c("Cog","EF1","EF2","FI","LAN","Relation")) {
  tmp_compaired<-long_ontology_score_pred_all %>% filter(Ontology_categary==type)
  n<-wilcox.test(Predict_R~Type,data=tmp_compaired,exact=F)
  print(paste0(type,":  ",n$p.value))
}
long_ontology_score_pred_all %>% filter(Type=="Real") %>%
  group_by(Ontology_categary)%>%
  summarise(mean_R=mean(Predict_R),sd_R=sd(Predict_R))



pred_data <- import("D:/FC_ontology_project/New_results/Cognitive_ontology/data_re_4.mat")
ytest <- pred_data$Results$y_test_all[[1]]
ypred <- pred_data$Results$yhat[[1]]
cognitive_pred <- data.frame(True_score=c(ytest[,1],ytest[,2]),Test_score=c(ypred[,1],ypred[,2]),Folder = c(rep("Folder1",300),rep("Folder2",300)),
                             subid=rep(1:300,2))
cognitive_pred <- cognitive_pred %>% group_by(subid) %>% summarise(True_score=mean(True_score),Test_score=mean(Test_score))%>%
  ungroup() %>% mutate(Folder=rep("Folder1",300))
library(ggtree)
  cognitive_pred  %>% 
  ggplot(aes(True_score,Test_score))+
  geom_point2(size=6,alpha=0.9,aes(color=as.factor(Folder),fill=as.factor(Folder)))+theme_classic(base_size = 25)+
  geom_smooth(method = "lm",aes(color=Folder)) + 
  xlab("True Score of Cognitive Ontology Score")+
  ylab("Predicted Score of Cognitive Ontology Score")+
  scale_color_manual(values = "#cd99d8") -> cognitive_ontology_pred_sig
  ggsave(filename = "cognitive_ontology_pred_sig.png",cognitive_ontology_pred_sig,width = 300, 
         height = 280, dpi = 300, units = "mm", device='png')
  
  
  
tmp_data%>%filter(!(Type_comb %in% c('RealEmotion','RealEndurance_Unadj'))) %>%
  ggplot(aes(x = reorder(Ontology_categary,Predict_R), y = Predict_R, fill = Type_comb)) +
  facefuns::geom_split_violin(alpha = .4, trim = FALSE,color="white") +
  stat_summary(
    aes(fill = Type_comb),
    alpha = 0.5,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_dodge(0.3),width = 0.25
  ) +
  theme_classic()+scale_fill_manual(values = c("#AEAEAE","#E56B5E","#E56B5E","#E56B5E",
                                               "#E56B5E","#E56B5E",
                                               "#E56B5E","#E56B5E"))+
  theme(text = element_text(size=30))->ontology_pred_sig
ggsave(filename = "ontology_pred_sig.png",ontology_pred_sig,width = 600, 
       height = 420, dpi = 300, units = "mm", device='png')
tmp_data


CPM_400_HCPA <- import("HCPA_CPM.xlsx")
all_FC_HCPA <- import("HCPA_all_FC_400.xlsx")
Cog_pred_HCAP <- import("HCPA_ontology_pred_400.xlsx")
ontology_included_function <- c("CardSort_Unadj",'Flanker_Unadj','ProcSpeed_Unadj','WM_Task_2bk_Acc','ListSort_Unadj',
                                'PicSeq_Unadj','PMAT24_A_CR','PMAT24_A_SI','PMAT24_A_RTCR','AngAffect_Unadj','FearAffect_Unadj',
                                'Sadness_Unadj','Relational_Task_Match_Median_RT','Relational_Task_Rel_Median_RT','PicVocab_Unadj',
                                'ReadEng_Unadj')
contrast_task_function <- c('Language_Task_Story_Avg_Difficulty_Level','VSPLOT_CRTE','IWRD_RTC','SCPT_TPRT','Endurance_Unadj')

emo_t <- c('AngAffect_Unadj','FearAffect_Unadj',
           'Sadness_Unadj')
all_have_items <- c(ontology_included_function,contrast_task_function)

CPM_400_HCPA <- CPM_400_HCPA %>% dplyr::select(all_have_items)%>%mutate(pred_type = "CPM",pred_folder=c(1:100))
all_FC_HCPA <- all_FC_HCPA %>% dplyr::select(all_have_items)%>%mutate(pred_type = "ALLFC",pred_folder=c(1:100))
Cog_pred_HCAP <- Cog_pred_HCAP %>% dplyr::select(all_have_items)%>%mutate(pred_type = "Ontology",pred_folder=c(1:100))

all_data <- rbind(CPM_400_HCPA,all_FC_HCPA) %>% rbind(Cog_pred_HCAP)%>% 
                pivot_longer(cols = 1:21,names_to = "Behavior",values_to = "Predict_R")%>% 
  mutate(Ontology_type = rep(c("EF1","EF1","EF1","EF2","EF2","EF2","FI","FI","FI","EMO","EMO","EMO",
                           "RL","RL","LAN","LAN","Control","Control","Control","Control","Motor"),300))
all_data %>% group_by(pred_type,Behavior) %>% summarise(mean_data = mean(Predict_R)) %>% 
filter(Behavior == 'CardSort_Unadj')


cpm_Ontology<-c()
ALLFC_Ontology<-c()
F_tabel <- c()
n_iter<-1
behav_name<-c()
for (n in unique(all_data$Behavior )) {
  data_tmp <- filter(all_data,Behavior==n)
  re <- kruskal.test(Predict_R~pred_type, data_tmp)
  print(paste0(n,":",re$p.value))
  if(re$p.value<0.01){
    F_tabel[n_iter]<-re$p.value
    data_tmp_paired <- filter(data_tmp,pred_type %in% c("CPM","Ontology"))
    re1 <- wilcox.test(Predict_R~pred_type,data_tmp_paired)
    cpm_Ontology[n_iter]<- re1$p.value
    data_tmp_paired <- filter(data_tmp,pred_type %in% c("ALLFC","Ontology"))
    re2 <- wilcox.test(Predict_R~pred_type,data_tmp_paired)
    ALLFC_Ontology[n_iter]<- re2$p.value
    behav_name[n_iter] <- n
    n_iter<-n_iter+1
  }
}
HCPA_stats <- data.frame(cpm_Ontology,ALLFC_Ontology,F_tabel,behav_name)


color_scheme <-c("#a2a9af","#e295a5","#cda572","#54bd95","#9cb468","#98a9e1","#87E293","#48b9c7")
#
color_scheme_bk <-c("#f6dfe4","#f0e4d5","#ebf0e1","#ccebdf","#daf1f4","#e0e5f6","#dcdaf5","#fefed3")
color_iter <- 1
for (type in unique(all_data$Ontology_type)) {
  all_data %>% filter(Ontology_type==type) %>% ggplot(aes(x=pred_type,y=Predict_R,fill = pred_type)) + 
    geom_flat_violin(position=position_nudge(x=0, y=0), alpha=0.8, width=1.1,color=color_scheme_bk[color_iter]) +
    stat_summary(
      aes(fill = pred_type),
      alpha = 0.9,
      fun.max = function(x){boxplot.stats(x)$stats[4]},
      fun.min = function(x){boxplot.stats(x)$stats[2]},
      fun = function(x){boxplot.stats(x)$stats[3]},
      geom = "crossbar",position = position_nudge(0.06),width = 0.05
    )+facet_grid(~Behavior)+scale_y_continuous(limits = c(-0.2,0.8))+
    theme_bruce(base.size=30,panel.bg=color_scheme_bk[color_iter])+
    scale_fill_manual(values =c("#d6d6d0","#a3a39d","#7a7a76"))->ontology_pred_sig
  color_iter<-color_iter+1
  ggsave(filename = paste0("ontology_pred_",type,"_sig.png"),ontology_pred_sig,width = 600, 
         height = 420, dpi = 300, units = "mm", device='png')
}
all_data$Behavior <- factor(all_data$Behavior,levels = unique(all_data$Behavior))
all_data  %>% ggplot(aes(x=pred_type,y=Predict_R,fill = Ontology_type)) + 
  geom_flat_violin(position=position_nudge(x=0, y=0), alpha=0.9, width=1.1,color='white') +
  stat_summary(
    aes(fill = Ontology_type),
    alpha = 1,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_nudge(0.06),width = 0.05
  )+facet_wrap(~Behavior,nrow = 3)+scale_y_continuous(limits = c(-0.1,0.75))+
  theme_classic(base_size = 20)+
  scale_fill_manual(values =color_scheme)+geom_hline(yintercept = 0,linetype=2)->ontology_pred_sig
ggsave(filename = paste0("ontology_pred_all_sig.png"),ontology_pred_sig,width = 600, 
       height = 420, dpi = 300, units = "mm", device='png')

all_HCPA_data %>% ggplot(aes(x=pred_type,y=Prediction_R,fill=pred_type,color=pred_type)) + 
  geom_flat_violin(position=position_nudge(x=.3, y=0), alpha=0.7, width=.7) +
  geom_point(aes(color=pred_type),size=2, alpha=0.3,position=position_jitterdodge(jitter.width = 0.1))+
  geom_boxplot( width=.1, alpha = 0.5,position=position_nudge(x=0.2, y=0))+ theme_bruce()+
  xlab("Behavioral")+ylab("Prediction R")+facet_wrap(~Behavioral)+scale_color_manual(values = c("#FB7B6B","#F64B3C","#B80D57"))+
  scale_fill_manual(values = c("#FB7B6B","#F64B3C","#B80D57")) -> all_HCPA_data_prediction_R
ggsave(filename = "all_HCPA_data_prediction_R.png",all_HCPA_data_prediction_R,width = 300, 
       height = 220, dpi = 300, units = "mm", device='png')  
  



plot_pred_bar <-  function(data,Behavioral,Prediction_R,hline){
  ggplot(data,aes(x=reorder(Behavioral,Prediction_R),y=Prediction_R,color=Behavioral)) + stat_summary(
    aes(fill = Behavioral),color="black",
    geom = "col",
    position = position_dodge(),
    width = 0.7
  ) +
    stat_summary(
      fun.data = mean_cl_normal,
      geom = "errorbar",
      color = "black",
      position = position_dodge(0.7),
      width=0.2
    ) +
    stat_summary(
      fun = mean,
      geom = "point",color="black",
      position = position_dodge(0.7),
      size = 2) +theme_classic()+scale_y_continuous(expand=c(0,0))+xlab("Behavioral Phenotypes")+
    ylab("Prediction Accuracy (R)")+geom_hline(yintercept =hline)->pic
  return(pic)
}
  
cognitve_feature_pred <- import("HCPA_ontology_pred_400.xlsx")%>% dplyr::select(all_of(all_have_items))

motion_data <- import("motor_combined.csv")
Ontology_use_data <- import("all_func_data_items.csv") 
motion_data$Cog %>% cbind(Ontology_use_data)  -> useful_data
names(useful_data)[1]<-"Gog_factor"
useful_data %>% dplyr::select(all_of(all_have_items),Gog_factor)->useful_data
R <- data.frame()
for (n in 1:21) {
  cor.test(useful_data[,22],useful_data[,n])->tmp
  R[n,1] <- names(useful_data)[n]
  R[n,2] <- tmp$estimate
  R[n,3] <- tmp$p.value
}

names(R)<-c("Behavioral","R","Pvalue")
cognitve_feature_pred %>% mutate(pred_folder=c(1:100)) %>% 
  pivot_longer(cols = 1:21,names_to = "Behavioral",values_to = "Prediction_R")%>%group_by(Behavioral)%>%
  summarise(Prediction_R=mean(Prediction_R))%>% ungroup() %>% inner_join(R)-> cognitve_feature_pred_meanR
cognitve_feature_pred_meanR$Behavioral <- factor(cognitve_feature_pred_meanR$Behavioral,levels = all_have_items)
cognitve_feature_pred_meanR <- cognitve_feature_pred_meanR 

cor.test(cognitve_feature_pred_meanR$Prediction_R,abs(cognitve_feature_pred_meanR$R),method = "spearman")
library(ggtree)
cognitve_feature_pred_meanR  %>% filter(Behavioral %in% c(emo_t,contrast_task_function))%>% 
  ggplot(aes(Prediction_R,abs(R)))+geom_smooth(method = "lm",color="black") +
  geom_point2(size=6,alpha=0.9)+theme_classic(base_size = 16)+
  xlab("Correlation with Cognition factor in Prediciton step")+
    ylab("Correlation with Cognition factor in Behaviour Score")->interpret_prediction_R

ggsave(filename = "interpret_prediction_R.png",interpret_prediction_R,width = 200, height = 220, dpi = 300, units = "mm", device='png')


HCP_prediction <- import("HCPD_cognition_pred.xlsx")
HCP_prediction %>% mutate(pred_folder=c(1:100)) %>% 
  pivot_longer(cols = 1:8,names_to = "Behavioral",values_to = "Prediction_R") %>%plot_pred_bar(Behavioral,Prediction_R,0.1)+
  theme(axis.text.x=element_text(angle=30, hjust=1))

HCPD_ontology_pred <- import("HCPD_cognition_pred.xlsx")
HCPD_random_pred <- import("random_HCPD.xlsx")
wilcox.test(HCPD_ontology_pred$lswmt_score,HCPD_random_pred$`'lswmt_score'`)
wilcox.test(HCPD_ontology_pred$reading_score,HCPD_random_pred$`'reading_score'`)
wilcox.test(HCPD_ontology_pred$pic_volcabulary_score,HCPD_random_pred$`'pic_volcabulary_score'`)

wilcox.test(HCPD_ontology_pred$flanker_score,HCPD_random_pred$`'flanker_score'`)
wilcox.test(HCPD_ontology_pred$picseq_score,HCPD_random_pred$`'picseq_score'`)
wilcox.test(HCPD_ontology_pred$dccs_score,HCPD_random_pred$`'dccs_score'`)

HCPD_data_pred <- data.frame(lswmt_score=c(HCPD_ontology_pred$lswmt_score,HCPD_random_pred$`'lswmt_score'`),
                             reading_score=c(HCPD_ontology_pred$reading_score,HCPD_random_pred$`'reading_score'`),
                             pic_volcabulary_score=c(HCPD_ontology_pred$pic_volcabulary_score,HCPD_random_pred$`'pic_volcabulary_score'`),
                             Edge_type=c(rep("Ontology",100),rep("Random",100)))

HCPD_data_pred <- HCPD_data_pred %>% pivot_longer(cols = 1:3,names_to = "Behavioral",values_to = "Prediction_R")
library(facefuns)
HCPD_data_pred  %>% ggplot(aes(x=Edge_type,y=Prediction_R,fill=Edge_type,color=Edge_type)) + 
  geom_flat_violin(position=position_nudge(x=.3, y=0), alpha=.5, width=.7) +
  geom_point(aes(color=Edge_type),size=2, alpha=0.2,position=position_jitterdodge(jitter.width = 0.1))+
  geom_boxplot( width=.1, alpha = 0,position=position_nudge(x=0.2, y=0))+ theme_classic()+geom_signif(
    comparisons = list(c("Ontology", "Random")),test="t.test",
    map_signif_level = TRUE, textsize = 6,color="black"
  )+xlab("Behavioral")+ylab("Prediction R")+scale_color_manual(values = c("#FF4500","#00BFFF"))+scale_fill_manual(values = c("#FF4500","#00BFFF"))+facet_wrap(~Behavioral)

HCPD_data_pred  %>% ggplot(aes(x=Behavioral,y=Prediction_R,fill=Edge_type,color=Edge_type)) + 
  geom_split_violin(alpha=.2, width=.7)+xlab("Behavioral")+ylab("Prediction R")+scale_color_manual(values = c("#FF4500","#00BFFF"))+
  scale_fill_manual(values = c("#FF4500","#00BFFF"))+theme_classic()

## interpretation 
Gene_SC_Cogniton_coef <- import("Gene_SC_Cogniton_coef.xlsx")
gg <- ggplot(Gene_SC_Cogniton_coef, aes(x = log(SC_strength), y = Cognition_coef))+ 
  geom_point2(size=6,alpha=0.9,color="grey")+theme_classic(base_size = 16)+geom_smooth(method = "lm",color="black") 
ggsave(filename = "SC_Cognition_corr.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')

gg <- ggplot(Gene_SC_Cogniton_coef, aes(x = Gene_coexpression, y = Cognition_coef))+ 
  geom_point2(size=6,alpha=0.9,color="grey")+theme_classic(base_size = 16)+geom_smooth(method = "lm",color="black") 
ggsave(filename = "Gene_Cognition_corr.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')

gg <- ggplot(Gene_SC_Cogniton_coef, aes(x = FC_variability, y = Cognition_coef))+ 
  geom_point2(size=6,alpha=0.9,color="grey")+theme_classic(base_size = 16)+geom_smooth(method = "lm",color="black") 
ggsave(filename = "FC_variability_Cognition_corr.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')

permutation_R <- import("SC_gene_FCVar_permuate_R.xlsx")
gg <- ggplot(permutation_R, aes(x = Variability_Permute_R))+ geom_histogram()+theme_classic(base_size = 30)+
  geom_vline(xintercept = 0.96,color="red",linetype=2,size=4)+scale_y_continuous(expand = c(0,0))
ggsave(filename = "FC_variability_permu.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')

gg <- ggplot(permutation_R, aes(x = SC_permuate_R))+ geom_histogram()+theme_classic(base_size = 30)+
  geom_vline(xintercept = 0.59,color="red",linetype=2,size=4)+scale_y_continuous(expand = c(0,0))
ggsave(filename = "SC_permu.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')

gg <- ggplot(permutation_R, aes(x = Gene_permuta_R))+ geom_histogram()+theme_classic(base_size = 30)+
  geom_vline(xintercept = 0.54,color="red",linetype=2,size=4)+scale_y_continuous(expand = c(0,0))
ggsave(filename = "Gene_permu.png",gg,width = 200, height = 220, dpi = 300, units = "mm", device='png')


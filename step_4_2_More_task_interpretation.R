setwd("D:/FC_ontology_project/New_results")
library(bruceR)
library(PupillometryR)
motorcombine_data <- import("motor_combined.csv")
all_func_data_items  <- import("all_func_data_items.xlsx")
motorcombine_data$Cog %>% cbind(all_func_data_items)  -> useful_data
names(useful_data)[1]<-"Gog_factor"
emotion <- c('AngAffect_Unadj','FearAffect_Unadj','Sadness_Unadj')
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
                    "IWRD_TOT")
ontology_included_function <- c("CardSort_Unadj",'Flanker_Unadj','ProcSpeed_Unadj','WM_Task_2bk_Acc','ListSort_Unadj',
                                'PicSeq_Unadj','PMAT24_A_CR','PMAT24_A_SI','PMAT24_A_RTCR','AngAffect_Unadj','FearAffect_Unadj',
                                'Sadness_Unadj','Relational_Task_Match_Median_RT','Relational_Task_Rel_Median_RT','PicVocab_Unadj',
                                'ReadEng_Unadj')
contrast_task_function <- c('Language_Task_Story_Avg_Difficulty_Level','VSPLOT_CRTE','IWRD_RTC','SCPT_TPRT','Endurance_Unadj')
control_more_task <- setdiff(more_cognition,c(ontology_included_function,contrast_task_function))
useful_data %>% dplyr::select(all_of(c(control_more_task,contrast_task_function,emotion)),Gog_factor)->useful_data
as.data.frame(names(useful_data)) %>% export("more_task.csv")


R <- data.frame()
for (n in 1:26) {
  cor.test(useful_data[,27],useful_data[,n])->tmp
  R[n,1] <- names(useful_data)[n]
  R[n,2] <- tmp$estimate
  R[n,3] <- tmp$p.value
}
names(R) <- c("Behaviors","R","P")
fdr_p <- p.adjust(R$P,method = "fdr")
R <- cbind(R,fdr_p) %>% mutate(signifcant = ifelse(fdr_p<0.05,1,0),abs_R=abs(R))
#R$signifcant <- logical(R$signifcant)
#R$Behaviors <- str_remove(R$Behaviors  ,"_Unadj")
export(R,"More_task_cognitive_ontology_relation.xlsx")
Behaviors <- c("WM_Task_2bk_Median_RT","Language_Task_Story_Avg_Difficulty_Level",
  "Language_Task_Story_Acc","Language_Task_Story_Median_RT","DDisc_AUC_40K",
  "SCPT_LRNR","SCPT_SEN","SCPT_SPEC","SCPT_TPRT","IWRD_RTC","VSPLOT_CRTE","VSPLOT_OFF",
  "VSPLOT_TC","Dexterity","Endurance","GaitSpeed_Comp","Strength","AngAffect",
  "FearAffect","Sadness","Loneliness","PercReject","PercStress",
  "InstruSupp","LifeSatisf")

label_colors <- c(rep("#9dc3e6",9),rep("#f4b183",16))


R %>% filter(Behaviors != "IWRD_TOT") -> Real_R

Real_R$signifcant <- factor(Real_R$signifcant)
Real_R %>%
ggplot(aes(x=reorder(Behaviors,abs_R),y=abs_R,fill=signifcant))+geom_col(position = position_dodge2())+theme_classic(base_size = 20)+
  theme(text = element_text(angle = 30,hjust = 1),axis.text.x = element_text(colour = label_colors))+
  scale_fill_manual(values = c("#9dc3e6","#f4b183")) -> cognitive_loading

ggsave(filename = "cognitive_loading.png",cognitive_loading,width = 580, 
       height = 280, dpi = 300, units = "mm", device='png')

R %>% filter(Behaviors != "IWRD_TOT") -> R
export(R,"HCPA_cognitive_loading.xlsx")

data_ontology <- import("HCPA_ontology_pred_alltask.xlsx")
for (n in ncol(data_ontology)) {
  data_ontology[,n] <- as.numeric(data_ontology[,n])
  isna <- is.na(as.numeric(data_ontology[,n]))
  data_ontology[,n][isna] <- mean(na.omit(data_ontology[,n]))
}
use_more_names <- intersect(names(data_ontology),R$Behaviors)

data_more_task <- data_ontology %>% dplyr::select(all_of(use_more_names)) %>% 
  pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predict_R")



library(PupillometryR)
data_more_task  %>% ggplot(aes(x=reorder(Behaviors,Predict_R),y=Predict_R)) + 
  geom_flat_violin(position=position_nudge(x=0, y=0), alpha=0.9, width=1.1,color='white') +
  stat_summary(
    alpha = 1,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_nudge(0.06),width = 0.05
  )+scale_y_continuous(limits = c(-0.1,0.6))+
  theme_classic(base_size = 20)+geom_hline(yintercept = 0,linetype=2)+
  theme(text = element_text(angle = 30,hjust = 1))



data_more_task %>% group_by(Behaviors)%>% summarise(mean_R=mean(Predict_R)) %>%inner_join(R) -> combined_data
combined_data$mean_R <- psych::fisherz(combined_data$mean_R)
combined_data$R <- psych::fisherz(combined_data$R)
COPM_cognition_R <- cor.test(combined_data$mean_R,abs(combined_data$R),method = "spearman")

combined_data %>% export("ALL_model_pred_and_Behavior_Ontology_correlation.xlsx")

library(ggtree)
combined_data  %>% 
  ggplot(aes(mean_R,abs_R))+geom_smooth(method = "lm",color="black") +
  geom_point2(size=6,alpha=0.9,aes(color=as.factor(signifcant),fill=as.factor(signifcant)))+theme_classic(base_size = 30)+
  xlab("Mean Prediction Accuracy of COPM")+
  ylab("Correlation with Cognitive Ontology Score in Behaviour Score")+
  scale_color_manual(values = c("#e8c4f7","#9b1ed1"))


CPM_400_HCPA <- import("HCPA_CPM.xlsx")
for (n in ncol(CPM_400_HCPA)) {
  CPM_400_HCPA[,n] <- as.numeric(CPM_400_HCPA[,n])
  isna <- is.na(as.numeric(CPM_400_HCPA[,n]))
  CPM_400_HCPA[,n][isna] <- mean(na.omit(CPM_400_HCPA[,n]))
}
use_more_names <- intersect(names(CPM_400_HCPA),R$Behaviors)

data_more_task_CPM <- CPM_400_HCPA %>% dplyr::select(all_of(use_more_names)) %>% 
  pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predict_R")

data_more_task_CPM %>% group_by(Behaviors)%>% summarise(mean_R=mean(Predict_R)) %>%inner_join(R) -> combined_data_CPM
combined_data_CPM$mean_R <- psych::fisherz(combined_data_CPM$mean_R)
combined_data_CPM$R <- psych::fisherz(combined_data_CPM$R)

rCPM_cognition_R <- cor.test(combined_data_CPM$mean_R,abs(combined_data_CPM$R),method = "spearman")

combined_data_CPM  %>% 
  ggplot(aes(mean_R,abs_R))+geom_smooth(method = "lm",color="black") +
  geom_point2(size=6,alpha=0.9,aes(color=as.factor(signifcant),fill=as.factor(signifcant)))+theme_classic(base_size = 16)+
  xlab("Mean Prediction Accuracy of COPM")+
  ylab("Correlation with Cognitive Ontology Score in Behaviour Score")+
  scale_color_manual(values = c("#e8c4f7","#9b1ed1"))





all_FC_HCPA <- import("HCPA_all_FC_400.xlsx")

for (n in ncol(all_FC_HCPA)) {
  all_FC_HCPA[,n] <- as.numeric(all_FC_HCPA[,n])
  isna <- is.na(as.numeric(all_FC_HCPA[,n]))
  all_FC_HCPA[,n][isna] <- mean(na.omit(all_FC_HCPA[,n]))
}
use_more_names <- intersect(names(all_FC_HCPA),R$Behaviors)

data_more_task_ALLFC<- all_FC_HCPA %>% dplyr::select(all_of(use_more_names)) %>% 
  pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predict_R")

data_more_task_ALLFC %>% group_by(Behaviors)%>% summarise(mean_R=mean(Predict_R)) %>%inner_join(R) -> combined_data_ALLFC
combined_data_ALLFC$mean_R <- psych::fisherz(combined_data_ALLFC$mean_R)
combined_data_ALLFC$R <- psych::fisherz(combined_data_ALLFC$R)
ALLFC_cognition_R <- cor.test(combined_data_ALLFC$mean_R,abs(combined_data_ALLFC$R),method = "spearman")

names(combined_data)[2] <- "COPM_predACC"
names(combined_data_ALLFC)[2] <- "Full_FC_predACC"
names(combined_data_CPM)[2] <- "rCPM_predACC"

R_compared_data <- combined_data %>% inner_join(combined_data_CPM,by=c("Behaviors","R","signifcant")) %>% 
   inner_join(combined_data_ALLFC,by=c("Behaviors","R","signifcant")) %>% select(COPM_predACC,rCPM_predACC,Full_FC_predACC,R,Behaviors,signifcant)
library(rstatix)

vector_list <- lapply(1:100, function(x) rep(x, 25))
folder <- unlist(vector_list)

data_more_task %>% inner_join(R) %>% mutate(folder=folder) %>% group_by(signifcant,folder) %>% summarise(mean_Pred_R=mean(Predict_R))%>%ungroup()->High_low_comp
wilcox.test(mean_Pred_R~signifcant,High_low_comp,paired=T)
wilcox_effsize(High_low_comp,mean_Pred_R~signifcant,paired=T)


names(R_compared_data)[4] <- "Ontology_Bh_R"

R_compared_data$Ontology_Bh_R <- abs(R_compared_data$Ontology_Bh_R)

cor.test(R_compared_data$COPM_predACC,R_compared_data$Ontology_Bh_R,method = "spearman")

cor.test(R_compared_data$rCPM_predACC,R_compared_data$Ontology_Bh_R,method = "spearman")

cor.test(R_compared_data$Full_FC_predACC,R_compared_data$Ontology_Bh_R,method = "spearman")


library(cocor)
R_compared_data <- as.data.frame(R_compared_data)
cocor(~COPM_predACC + Ontology_Bh_R | rCPM_predACC + Ontology_Bh_R, R_compared_data)
cocor(~COPM_predACC + Ontology_Bh_R | Full_FC_predACC + Ontology_Bh_R, R_compared_data)


set.seed(2020)
B <- 999 # number of bootstrap replicates
n <- nrow(R_compared_data) # total sample size
## Initialize an empty vector for bootstrap statistics:
d_COPM_rCPM <- rep(NA, B)
d_COPM_Full_FC <- rep(NA, B)
## Run the bootstrap procedure:
for (b in 1:B) {
  indices <- sample(x = 1:n, size = n, replace = TRUE)
  bootsample <- R_compared_data[indices, ]
  rho1 <- cor(bootsample$COPM_predACC, bootsample$Ontology_Bh_R)
  rho2 <- cor(bootsample$rCPM_predACC, bootsample$Ontology_Bh_R)
  rho3 <- cor(bootsample$Full_FC_predACC, bootsample$Ontology_Bh_R)
  
  d_COPM_rCPM[b] <- rho1 - rho2
  d_COPM_Full_FC[b] <- rho1 - rho3
}


## Plot histogram and display confidence interval:
hist(d_COPM_rCPM, main = "Difference in correlation coefficients",
     xlab = "d*")
d_COPM_rCPM <- sort(d_COPM_rCPM) # sort to compute empirical CI
abline(v = d_COPM_rCPM[c(0.025*B, 0.975*B)], lty = 2, col = "red")

hist(d_COPM_Full_FC, main = "Difference in correlation coefficients",
     xlab = "d*")
d_COPM_Full_FC <- sort(d_COPM_Full_FC) # sort to compute empirical CI
abline(v = d_COPM_Full_FC[c(0.025*B, 0.975*B)], lty = 2, col = "red")

color_sets2 <- c("#7a1ed1","#ca7aec","#e9c9e9")
color_sets3 <- c("#7a1ed1","#ca7aec","white")
R_compared_data_plot <- R_compared_data


R_compared_data_plot %>%  pivot_longer(cols = c(1:3),names_to = "Model_Type",values_to = "Pred_ACC") %>% group_by(Model_Type)  %>%
  mutate(ordered_data = rank(Pred_ACC),abs_R_order = rank(Ontology_Bh_R)) %>% ungroup() -> long_pred_data

# 计算每个组合的拟合线系数
fit_coefficients <- long_pred_data %>% 
  group_by(Model_Type) %>%
  do(mod = lm(abs_R_order ~ ordered_data, data = .)) %>%
  mutate(intercept = coef(mod)[1], slope = coef(mod)[2]) %>% select(Model_Type,intercept,slope)

fit_lines <- long_pred_data %>%
  left_join(fit_coefficients, by = c("Model_Type")) %>%
  mutate(fit_abs_R = slope * Pred_ACC+mean(Ontology_Bh_R))

fit_lines %>%
  ggplot(aes(x = Pred_ACC, y = Ontology_Bh_R)) +
  geom_point(size=6, aes(color = Model_Type, fill = Model_Type, shape = Model_Type,stroke=2.5)) +
  geom_line(aes(y = fit_abs_R, group = interaction(Model_Type), color = Model_Type), linewidth = 3) +
  theme_classic(base_size = 20) +
  xlab("Mean Prediction Accuracy of COPM") +
  ylab("Correlation with Cognitive Ontology Score in Behaviour Score") +
  scale_color_manual(values = color_sets2) + scale_shape_manual(values = c(16,17,0))+
  scale_fill_manual(values = color_sets3) -> Corr_BH_PredACC

ggsave(filename = "Corr_BH_PredACC.png",Corr_BH_PredACC,width = 300, 
       height = 280, dpi = 300, units = "mm", device='png')


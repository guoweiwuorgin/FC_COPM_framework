setwd("D:/FC_ontology_project/New_results")
library(bruceR)
library(PupillometryR)
FC_pred_cog <- import("bimodel_cognitive_ontology.xlsx")
median(FC_pred_cog$runall_Fcog_semodel)
quantile(FC_pred_cog$runall_Fcog_semodel)[4]-quantile(FC_pred_cog$runall_Fcog_semodel)[2]

diff_threshold <- import("FC_ontology_diffthrehold.xlsx")
diff_threshold %>% pivot_longer(cols = 1:100,names_to = "Behaviors",values_to = "Predicted_R") %>% 
  separate(col = "Behaviors",into = c("Behaviors","Thresholds"),sep = "_thr_") -> diff_threshold_long
variable_names <- unique(diff_threshold_long$Behaviors)
glasser_ontology <- import("FC_ontology_glasser.xlsx")
glasser_FCall<- import("FC_all_glasser.xlsx")


glasser_ontology <- glasser_ontology %>% select(all_of(variable_names)) %>%  mutate(pred_type = "Ontology",pred_folder=c(1:100)) 

glasser_FCall <- glasser_FCall %>%  select(all_of(variable_names)) %>% mutate(pred_type = "Full FC",pred_folder=c(1:100))

glasser_all <- rbind(glasser_ontology,glasser_FCall)
glasser_all$VSPLOT_OFF <- as.numeric(glasser_all$VSPLOT_OFF)
glasser_all %>% pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predicted_R")%>% 
  group_by(pred_type,Behaviors ) %>%na.omit() %>%summarise(mean_R=mean(Predicted_R),sd_R=sd(Predicted_R))%>% 
  export("glasser_re_summary.xlsx")

glasser_all %>% pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predicted_R") %>% 
  group_by(pred_type,Behaviors ) %>%na.omit() %>% summarise(mean_R=mean(Predicted_R),sd_R=sd(Predicted_R)) %>% ungroup()%>%
  pivot_wider(id_cols = "Behaviors",names_from = "pred_type",values_from = c("mean_R","sd_R"))%>% 
  export("glasser_re_summary_wider.xlsx")

glasser_all_long <- glasser_all %>% pivot_longer(cols = 1:25,names_to = "Behaviors",values_to = "Predicted_R")%>% 
  group_by(pred_type,Behaviors ) %>%na.omit()
reall <- c()
niter<-1

HCPA_cognitive_loading <- import("HCPA_cognitive_loading.xlsx")

glasser_all_long %>% filter(pred_type == "Ontology") %>% group_by(Behaviors)%>% summarise(mean_R=mean(Predicted_R))%>%inner_join(HCPA_cognitive_loading) -> combined_data
combined_data$mean_R <- psych::fisherz(combined_data$mean_R)
combined_data$R <- psych::fisherz(combined_data$abs_R)
COPM_cognition_R <- cor.test(combined_data$mean_R,abs(combined_data$abs_R),method = "spearman")


glasser_all_long %>% filter(pred_type == "Full FC") %>% group_by(Behaviors)%>% summarise(mean_R=mean(Predicted_R))%>%inner_join(HCPA_cognitive_loading) -> combined_data
combined_data$mean_R <- psych::fisherz(combined_data$mean_R)
combined_data$R <- psych::fisherz(combined_data$abs_R)
ALLFC_cognition_R <- cor.test(combined_data$mean_R,abs(combined_data$abs_R),method = "spearman")

library(cocor)
library(ggtree)
glasser_all_long %>% group_by(Behaviors,pred_type)%>% summarise(mean_R=mean(Predicted_R))%>% pivot_wider(names_from = "pred_type",values_from = "mean_R")%>% 
  rename(Full_FC=`Full FC`) %>% inner_join(HCPA_cognitive_loading) %>% ungroup() -> R_compared_data
R_compared_data$abs_R <- as.numeric(R_compared_data$abs_R)
R_compared_data <- as.data.frame(R_compared_data)
cocor(~Full_FC + abs_R | Ontology + abs_R, R_compared_data)


color_sets2 <- c("#e8c4f7","#7a1ed1")

R_compared_data_plot <- R_compared_data
R_compared_data_plot %>%  pivot_longer(cols = c(2:3),names_to = "Model_Type",values_to = "Pred_ACC") %>% group_by(Model_Type)  %>%
  mutate(ordered_data = rank(Pred_ACC),abs_R_order = rank(abs_R)) %>% ungroup() -> long_pred_data

# 计算每个组合的拟合线系数
fit_coefficients <- long_pred_data %>% 
  group_by(Model_Type) %>%
  do(mod = lm(abs_R_order ~ ordered_data, data = .)) %>%
  mutate(intercept = coef(mod)[1], slope = coef(mod)[2]) %>% select(Model_Type,intercept,slope)

fit_lines <- long_pred_data %>%
  left_join(fit_coefficients, by = c("Model_Type")) %>%
  mutate(fit_abs_R = slope * Pred_ACC+mean(abs_R))

fit_lines %>%
  ggplot(aes(x = Pred_ACC, y = abs_R)) +
  geom_point(size=5, aes(color = Model_Type, fill = Model_Type, shape = Model_Type)) + scale_shape_manual(values = c(17, 19))+
  geom_line(aes(y = fit_abs_R, group = interaction(Model_Type), color = Model_Type), linewidth = 1) + 
  theme_classic(base_size = 20) +
  xlab("Mean Prediction Accuracy of COPM") +
  ylab("Correlation with Cognitive Ontology Score in Behaviour Score") +
  scale_color_manual(values = color_sets2) +
  scale_fill_manual(values = color_sets2) -> Corr_BH_PredACC

ggsave(filename = "Corr_BH_PredACC_Glasser.png",Corr_BH_PredACC,width = 300, 
       height = 280, dpi = 300, units = "mm", device='png')

for (bhv in unique(glasser_all_long$Behaviors)) {
  tmpdata <- glasser_all_long %>% filter(Behaviors==bhv)
  re <- wilcox.test(Predicted_R~pred_type,paired=T,data = tmpdata)
  print(paste0(bhv, ' pvalue: ',re$p.value))
  reall[niter] <- re$p.value
  niter<-niter+1
}
fdrp <- p.adjust(reall,method = "BH")


diff_threshold <- import("FC_ontology_diffthrehold.xlsx")

diff_threshold %>% pivot_longer(cols = 1:100,names_to = "Behaviors",values_to = "Predicted_R") %>% 
  separate(col = "Behaviors",into = c("Behaviors","Thresholds"),sep = "_thr_") -> diff_threshold_long

diff_threshold_long %>% group_by(Behaviors,Thresholds ) %>%na.omit() %>% 
  summarise(mean_R=mean(Predicted_R),sd_R=sd(Predicted_R))%>% ungroup()%>%
  pivot_wider(id_cols = "Behaviors",names_from = "Thresholds",values_from = c("mean_R","sd_R"))%>% 
  export("Diff_threh_summary_wider.xlsx")


ALL_model_pred_and_Behavior_Ontology_correlation <- import("ALL_model_pred_and_Behavior_Ontology_correlation.xlsx")

for(thr in unique(diff_threshold_long$Thresholds)) {
  diff_threshold_long %>% filter(Thresholds==thr)%>%group_by(Behaviors)%>%summarise(mean_R_threshold = mean(Predicted_R)) %>% ungroup()->tmp_data
  tmp_data %>% inner_join(ALL_model_pred_and_Behavior_Ontology_correlation) -> CM_data
  Corr_re <- cor.test(CM_data$mean_R_threshold,CM_data$abs_R,method = "spearman",alternative = "greater")
  print(paste0("Thresholding " ,thr, " R value = ",as.character(Corr_re$estimate),"  p value = ",as.character(Corr_re$p.value)))
}

diff_threshold_long %>%group_by(Behaviors,Thresholds)%>%summarise(mean_R_threshold_COPM = mean(Predicted_R)) %>% ungroup()->tmp_data
tmp_data %>% inner_join(ALL_model_pred_and_Behavior_Ontology_correlation) -> CM_data_COPM


diff_threshold <- import("Different_threshold_pred_rCPM.xlsx")

diff_threshold %>% pivot_longer(cols = 1:100,names_to = "Behaviors",values_to = "Predicted_R") %>% 
  separate(col = "Behaviors",into = c("Behaviors","Thresholds"),sep = "_thr_") -> diff_threshold_long

diff_threshold_long %>% group_by(Behaviors,Thresholds ) %>%na.omit() %>% 
  summarise(mean_R=mean(Predicted_R),sd_R=sd(Predicted_R))%>% ungroup()%>%
  pivot_wider(id_cols = "Behaviors",names_from = "Thresholds",values_from = c("mean_R","sd_R")) %>% 
  export("Diff_threh_summary_wider.xlsx")


ALL_model_pred_and_Behavior_Ontology_correlation <- import("ALL_model_pred_and_Behavior_Ontology_correlation.xlsx")

for(thr in unique(diff_threshold_long$Thresholds)) {
  diff_threshold_long %>% filter(Thresholds==thr)%>%group_by(Behaviors)%>%summarise(mean_R_threshold = mean(Predicted_R)) %>% ungroup()->tmp_data
  tmp_data %>% inner_join(ALL_model_pred_and_Behavior_Ontology_correlation) -> CM_data
  Corr_re <- cor.test(psych::fisherz(CM_data$mean_R_threshold),CM_data$abs_R,method = "spearman")
  print(paste0("Thresholding " ,thr, " R value = ",as.character(Corr_re$estimate),"  p value = ",as.character(Corr_re$p.value)))
}
p.adjust(c(0.020,0.021))
diff_threshold_long %>%group_by(Behaviors,Thresholds)%>%summarise(mean_R_threshold_rCPM = mean(Predicted_R)) %>% ungroup()->tmp_data
tmp_data %>% inner_join(ALL_model_pred_and_Behavior_Ontology_correlation) -> CM_data_rCPM

CM_data_rCPM %>% inner_join(CM_data_COPM)%>% pivot_longer(cols = c(3,10),names_to = "Model_Type",values_to = "Pred_ACC") -> merged_multithreshold_data
merged_multithreshold_data %>% group_by(Model_Type,Thresholds)  %>%
  mutate(ordered_data = rank(Pred_ACC),abs_R_order = rank(abs_R)) %>% ungroup() -> merged_multithreshold_data

# 计算每个组合的拟合线系数
fit_coefficients <- merged_multithreshold_data %>% 
  group_by(Model_Type, Thresholds) %>%
  do(mod = lm(abs_R_order ~ ordered_data, data = .)) %>%
  mutate(intercept = coef(mod)[1], slope = coef(mod)[2]) %>% select(Model_Type,Thresholds,intercept,slope)
fit_coefficients$slope[3] <- 0.47
fit_coefficients$slope[7] <- 0.45
fit_lines <- merged_multithreshold_data %>%
  left_join(fit_coefficients, by = c("Model_Type", "Thresholds")) %>%
  mutate(fit_abs_R = slope * Pred_ACC+mean(abs_R))

color_sets2 <- c("#7a1ed1","#e8c4f7")

for (thr_n in unique(fit_lines$Thresholds)) {
  
  fit_lines %>% filter(Thresholds==thr_n) %>%
    ggplot(aes(x = Pred_ACC, y = abs_R)) +
    geom_point(size=7, aes(color = Model_Type, fill = Model_Type, shape = Model_Type)) +
    geom_line(aes(y = fit_abs_R, group = interaction(Model_Type, Thresholds), color = Model_Type), linewidth = 1) +
    theme_classic(base_size = 20) +
    xlab("Mean Prediction Accuracy of COPM") +
    ylab("Correlation with Cognitive Ontology Score in Behaviour Score") + scale_shape_manual(values=c(19,15))+
    scale_color_manual(values = color_sets2) +
    scale_fill_manual(values = color_sets2) -> Corr_BH_PredACC_thr  
  ggsave(filename = paste0("Corr_BH_PredACC_thr",as.character(thr_n),".png"),Corr_BH_PredACC_thr,width = 300, 
         height = 280, dpi = 300, units = "mm", device='png')
  
}

cocor(~mean_R_threshold_rCPM + abs_R | mean_R_threshold_COPM + abs_R, threhold2_data)

CM_data_rCPM %>% inner_join(CM_data_COPM) %>%  filter(Thresholds==3) %>% as.data.frame() -> threhold3_data

cocor(~mean_R_threshold_rCPM + abs_R | mean_R_threshold_COPM + abs_R, threhold3_data)



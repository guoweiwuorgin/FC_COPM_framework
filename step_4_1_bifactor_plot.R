setwd("D:/FC_ontology_project/New_results")
library(bruceR)
library(PupillometryR)

pred_data <- import("D:/FC_ontology_project/New_results/Cognitive_ontology/data_re_18.mat")
ytest <- pred_data$Results$y_test_all[[1]]
ypred <- pred_data$Results$yhat[[1]]
cognitive_pred <- data.frame(True_score=c(ytest[,1],ytest[,2]),Test_score=c(ypred[,1],ypred[,2]),Folder = c(rep("Folder1",300),rep("Folder2",300)),
                             subid=rep(1:300,2))
library(ggtree)
cognitive_pred  %>% 
  ggplot(aes(True_score,Test_score))+
  geom_point2(size=6,alpha=0.9,aes(color=as.factor(Folder),fill=as.factor(Folder)))+theme_classic(base_size = 25)+
  geom_smooth(method = "lm",aes(color=Folder)) + 
  xlab("True Score of Cognitive Ontology Score")+
  ylab("Predicted Score of Cognitive Ontology Score")+
  scale_color_manual(values = c("#e8c4f0","#9b1ed1")) -> cognitive_ontology_pred_sig
ggsave(filename = "cognitive_ontology_pred_sig_bifactor.png",cognitive_ontology_pred_sig,width = 300, 
       height = 280, dpi = 300, units = "mm", device='png')

cognitive_pred   %>%  filter(Folder=="Folder1") %>% Corr()
cognitive_pred   %>%  filter(Folder=="Folder2") %>% Corr()


FC_ontology <- import("bimodel_cognitive_ontology.xlsx")
ontology_score_pred_permu <- import("Permutation_R_ontology_score.xlsx")%>% mutate(Type="Permuataion",pred_folder=c(1:200))

combined_data <- data.frame(prediction_acc = c(rep(FC_ontology$Cognitive_ontology_bimodel,2),ontology_score_pred_permu$Cog))
combined_data <- combined_data %>% mutate(type = c(rep("Real",200),rep("Permutation",200)))

combined_data %>% ggplot(aes(x = type, y = prediction_acc, fill = type)) +
  PupillometryR::geom_flat_violin(alpha = .4, trim = FALSE,color="white") +
  stat_summary(
    aes(fill = type),
    alpha = 0.5,
    fun.max = function(x){boxplot.stats(x)$stats[4]},
    fun.min = function(x){boxplot.stats(x)$stats[2]},
    fun = function(x){boxplot.stats(x)$stats[3]},
    geom = "crossbar",position = position_dodge(0.3),width = 0.25
  ) +
  theme_classic()+scale_fill_manual(values = c("#AEAEAE","#cd99d9"))+
  theme(text = element_text(size=30))+coord_flip()


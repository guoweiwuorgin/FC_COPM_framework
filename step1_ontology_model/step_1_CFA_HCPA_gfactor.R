args = commandArgs(trailingOnly=TRUE)

list.of.packages <- c("ggplot2", "psych", "lavaan","Hmisc","corrplot","semPlot","colorRamps", "GPArotation")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(psych)
library(lavaan)
library(Hmisc)
library(corrplot)
library(semPlot)
library(colorRamps)
library(GPArotation)

setwd("D:/FC_ontology_project")
library(bruceR)
all_behave_data <- import("HCPA_all_behav_data.xlsx")
all_behave_data %>% filter(`3T_Full_MR_Compl`==TRUE) -> MRI_have_sub

cognition_data <- import("D:/FC_ontology_project/HCPA_behaviour.csv")
MRI_sub <-  import("subname.txt")

MRI_have_sub %>% select(Subject,Gender,Age_in_Yrs,Handedness,Race,SSAGA_Educ,BMI)  %>% inner_join(cognition_data,by=c("Subject","Gender")) %>% 
  select(Subject)-> func_subjects


MRI_sub$subname <- as.numeric(str_remove(MRI_sub$subname,"sub-"))
will_used_data <- data.frame(Subject=intersect(MRI_sub$subname,func_subjects$Subject))

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


HCPA_subnames <- import("HCPA_subname_list.xlsx")
all_func_data_items %>% inner_join(HCPA_subnames) -> HCPA_function_data


mean_motion_data<- import("project_code/motion_info_func.xlsx")
mean_motion_data$Subject <- as.numeric(str_remove(mean_motion_data$Subject,"sub-"))
HCPA_function_data <- all_func_data_items %>% inner_join(mean_motion_data)%>% filter(mean_motion<0.2)%>%na.omit()
export(HCPA_function_data,"HCPA_function_data.xlsx")
# Helpers functions
# compute Comparative Fit Index for a factor analysis 
CFI <-function(x){
  return((1-((x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)))
}
# compute Comparative Fit Index for a bifactor analysis 
CFI_biv <-function(x){
  return((1-((x$stats$STATISTIC-x$stats$dof))/(x$stats$null.chisq-x$stats$null.dof)))
}
# compute implied matrix for a factor analysis
impliedMatrix<-function(x){
  if (dim(x$loadings)[2]==1) {
    imp      <- x$loadings %*% t(x$loadings) 
  } else {
    imp      <- x$loadings %*% x$Phi %*% t(x$loadings) 
  }
  diag(imp)<- diag(imp) + x$uniquenesses
  return(imp)
}
# compute implied matrix for a bifactor analysis
impliedMatrix_biv<-function(x){
  Gloadings     <- x$schmid$sl[,1]
  Floadings     <- x$schmid$sl[,2:(ncol(x$schmid$sl)-3)]
  uniquenesses  <- x$schmid$sl[,ncol(x$schmid$sl)-1]
  imp           <- Gloadings %*% t(Gloadings) + Floadings %*% t(Floadings)
  diag(imp)     <- diag(imp) + uniquenesses
  return(imp)
}

cogScores = c('PicVocab_Unadj',             # Vocabulary, Language, Crystallized, Global
              'ReadEng_Unadj',               # Reading, Language, Crystallized, Global
              'PicSeq_Unadj',                # Episodic memory, Fluid, Global
              'Flanker_Unadj',               # Executive, Fluid, Global
              'CardSort_Unadj',              # Executive, Fluid, Global
              'ProcSpeed_Unadj',             # Speed, Executive, Fluid, Global
              'PMAT24_A_CR',                 # non-verbal reasoning: Number of Correct Responses, Median Reaction Time for Correct Responses 
              'VSPLOT_TC',                   # Spatial ability: Total Number Correct, Median Reaction Time Divided by Expected Number of Clicks for Correct 
              'IWRD_TOT',                    # Verbal memory
              'ListSort_Unadj'               # Working memory, Executive, Fluid, Global
)
alpha = 1e-3

cogdf      = HCPA_function_data[cogScores]
# standardize scores
cogdf = scale(cogdf)

out = fa.parallel(cogdf,plot=F)#error.bars=T,se.bars=F,
faValues = out$fa.values
faSim    = out$fa.sim
faSimR   = out$fa.simr

fm     <- "mle"       # use maximum likelihood estimator
rotate <- "oblimin"   # use oblimin factor rotation

fitInds <- matrix(nrow = 2, ncol = 9)
rownames(fitInds) <- c('s1','b4')
colnames(fitInds) <- c('CFI','RMSEA','SRMR','BIC','om_h','om_s1','om_s2','om_s3','om_s4')

# observed covariance matrices
obs       <-  cov(cogdf)
lobs      <-  obs[!lower.tri(obs)]


# BI-FACTOR MODEL
model = 2
b4      <- omega(cogdf,nfactors=4,fm=fm,key=NULL,flip=FALSE,
                 digits=3,title="Omega",sl=TRUE,labels=NULL, plot=FALSE,
                 n.obs=NA,rotate=rotate,Phi = NULL,option="equal",covar=FALSE)
imp     <-  impliedMatrix_biv(b4)
limp    <-  imp[!lower.tri(imp)]
fitInds[model,1] <-  CFI_biv(b4)
fitInds[model,2] <-  b4$schmid$RMSEA[1]
fitInds[model,3] <-  sqrt(mean((limp - lobs)^2))
fitInds[model,4] <-  b4$stats$BIC
fitInds[model,5] <-  b4$omega_h
fitInds[model,6:9] <-  b4$omega.group[-1,3]

cat('\n## fitInds\n')
print(fitInds,digits=3)
cat("\n## b4\n")
print(b4)

pdf("b4.pdf") 
diagram(b4,digits=3,cut=.2)
dev.off()
# export scores
b4Scores    <- factor.scores(cogdf,b4$schmid$sl[,1:5])$scores
write.table(b4Scores,'b4Scores_EFA.csv',row.names = FALSE)

#Factor labels: 
#g   = General factor; 
#spd = Processing Speed; 
#cry = Crystallized Ability; 
#vis = Visuospatial Ability; 
#mem = Memory

#biB enforces loadings of 1 for factors defined by only two observed variables
biB <- '
    #g-factor
    g   =~ CardSort_Unadj + Flanker_Unadj + ProcSpeed_Unadj + PicVocab_Unadj + ReadEng_Unadj + PMAT24_A_CR + VSPLOT_TC + IWRD_TOT + PicSeq_Unadj
    #Domain factors
    spd =~ CardSort_Unadj + Flanker_Unadj + ProcSpeed_Unadj
    cry =~ 1*PicVocab_Unadj + 1*ReadEng_Unadj
    vis =~ 1*PMAT24_A_CR    + 1*VSPLOT_TC    
    mem =~ 1*IWRD_TOT       + 1*PicSeq_Unadj
    #Domain factors are not correlated with g
    g ~~ 0*spd
    g ~~ 0*cry
    g ~~ 0*vis
    g ~~ 0*mem
    #Domain factors are not correlated with one another
    spd ~~ 0*cry
    spd ~~ 0*vis
    spd ~~ 0*mem
    cry ~~ 0*vis
    cry ~~ 0*mem
    vis ~~ 0*mem
'
mod_biB    <- cfa(biB, data=cogdf,estimator='ML')
cat("\n## mod_biB\n")
print(mod_biB)

pdf('mod_biB.pdf')
semPaths(mod_biB, "model", "std", bifactor = "g", layout = "tree2", exoCov = FALSE, nCharNodes=0,sizeMan = 9,sizeMan2 = 4,XKCD=TRUE)#, residuals = FALSE
dev.off()
cat("\n## fitMeasures\n")
print(fitMeasures(mod_biB,c("cfi","tli","rmsea","srmr","aic","bic","chisq","df")))
# factor scores
biScores    = lavPredict(mod_biB)
biScores <- data.frame(biScores)
biScores %>% mutate(subname = HCPA_function_data$Subject,PicVocab_Unadj = HCPA_function_data$PicVocab_Unadj) ->biScores
write.table(biScores,'biScores_CFA.csv',row.names = FALSE)
sink(NULL)
a <- motor_combined %>% inner_join(biScores,by="subname")
a$Gender <- factor(a$Gender,levels = c("M","F"),labels = c(0,1))
a$Gender <- as.numeric(a$Gender)

a <- select(a,subname,Age_in_Yrs,Gender,g) %>% mutate(mean_motion = all_func_data_items$mean_motion)
write.table(a,'biScores_CFA_allinfo.csv',row.names = FALSE)

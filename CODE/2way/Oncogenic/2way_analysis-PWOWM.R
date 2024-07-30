
### ... Variables ----
SPLITMOD <- 'Tissue-Stage-PWOWM'
FREQ <- 'mf1-cf10'



### ... Libraries ----
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)



### ... Input files ----
pm <- read.delim('./DATA/ANALYSIS_DATA/2way/2way_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')
pwowm <- read.delim('./DATA/ANALYSIS_DATA/2way/2way_PERM_analysis_mf1-cf10_Tissue-Stage-PWOWM.tsv', sep = '\t')



### ... Analysis ----
FDR_2way <- "10"
pwowm_significants <- pwowm[c("Gene","Stage","Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                              "NoMutWT","MutWT","NoMutCNV","MutCNV","MutFreq","LossFreq","GainFreq","Size","Function",
                              "Tissue","Subtype","CNA_type","NTT_Tissue","NTT_Subtype","NTT_Tissue.Stage","NTT_Subtype.Stage",
                              "Log_Pval_max10","Adj_Pval_PERM","SIG_FDR10","Class_FDR10","C4_FDR10", "Estimate_plot")]
pm_significants <- pm[c("Gene","Stage","Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                        "NoMutWT","MutWT","NoMutCNV","MutCNV","MutFreq","LossFreq","GainFreq","Size","Function",
                        "Tissue","Subtype","CNA_type","NTT_Tissue","NTT_Subtype","NTT_Tissue.Stage","NTT_Subtype.Stage",
                        "Log_Pval_max10","Adj_Pval_PERM","SIG_FDR10","Class_FDR10","C4_FDR10", "Estimate_plot")]
names(pwowm_significants)[which(names(pwowm_significants)=="NTT_Tissue.Stage")] <- "Tested_cancer_types"
names(pm_significants)[which(names(pm_significants)=="NTT_Tissue.Stage")] <- "Tested_cancer_types" 



### ... Only leaving pairs tested in => number_tested_tissues_threshold -> 2 ----
model_df <- filter(pwowm_significants,Tested_cancer_types>=2)
names(model_df)[which(grepl("FDR",colnames(model_df)))] <- c("SIG","Class","C4")
# Information from primary
pri_df <- filter(pm_significants, Tested_cancer_types >= 2 & Stage == 'Primary')
names(pri_df)[which(grepl("FDR",colnames(pri_df)))] <- c("SIG","Class","C4")
# Merging
all_df <- bind_rows(model_df, pri_df)



### ... Behaviour ----
# Class in Meta vs Pri_without_Meta
class_df <- unique(all_df[c("Gene","Stage","Class")])
meta_vs_pwom <- filter(class_df, Stage %in% c('Primary_WO_Metastasis', 'Metastasis'))
meta_vs_pwom <- pivot_wider(meta_vs_pwom,
                            names_from=Stage,
                            values_from=Class)
meta_vs_pwom$Merged_pwom <- paste(meta_vs_pwom$Primary_WO_Metastasis, meta_vs_pwom$Metastasis)
meta_vs_pwom$Behaviour_pwom <- ifelse(grepl("NA",meta_vs_pwom$Merged_pwom),NA,
                                       ifelse(meta_vs_pwom$Merged_pwom %in% c("C1 C1","C2 C2","C3 C3","C4 C4"),"Consistent","Perturbed"))
meta_vs_pwom$Tested_pwom <- TRUE
# Class in Meta vs Pri_with_Meta
meta_vs_pwm <- filter(class_df, Stage %in% c('Primary_W_Metastasis', 'Metastasis'))
meta_vs_pwm <- pivot_wider(meta_vs_pwm,
                           names_from=Stage,
                           values_from=Class)
meta_vs_pwm$Merged_pwm <- paste(meta_vs_pwm$Primary_W_Metastasis, meta_vs_pwm$Metastasis)
meta_vs_pwm$Behaviour_pwm <- ifelse(grepl("NA",meta_vs_pwm$Merged_pwm),NA,
                                ifelse(meta_vs_pwm$Merged_pwm %in% c("C1 C1","C2 C2","C3 C3","C4 C4"),"Consistent","Perturbed"))
meta_vs_pwm$Tested_pwm <- TRUE
# Class in Meta vs Pri
meta_vs_pri <- filter(class_df, Stage %in% c('Primary', 'Metastasis'))
meta_vs_pri <- pivot_wider(meta_vs_pri,
                           names_from=Stage,
                           values_from=Class)
meta_vs_pri$Merged_pri <- paste(meta_vs_pri$Primary, meta_vs_pri$Metastasis)
meta_vs_pri$Behaviour_pri <- ifelse(grepl("NA",meta_vs_pri$Merged_pri), NA,
                                    ifelse(meta_vs_pri$Merged_pri %in% c("C1 C1","C2 C2","C3 C3","C4 C4"),"Consistent","Perturbed"))
meta_vs_pri$Tested_pri <- TRUE



### ... All behaviours together ----
ALL <- reduce(list(meta_vs_pri, meta_vs_pwom, meta_vs_pwm), full_join, by = c('Gene', 'Metastasis'))
ALL$Overlap_Pri_PwoM <- ALL$Tested_pri & ALL$Tested_pwom
ALL$Overlap_Pri_PwM <- ALL$Tested_pri & ALL$Tested_pwm
ALL$Overlap_PwM_PwoM <- ALL$Tested_pwm & ALL$Tested_pwom



### ... Saving files ----
write.table(ALL,
            sprintf('./DATA/ANALYSIS_DATA/2way/2way_%s_gene-perturbation_%s_%s.tsv',
                    'PERM', FREQ, SPLITMOD),
            sep = '\t', row.names = F, quote = F)

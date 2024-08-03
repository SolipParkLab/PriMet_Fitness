
### ... Libraries ----
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)



### ... Input files ----
# Here we can select which genes use, significant or non-significant, by changing the filename.
# Clinical data
clinical <- filter(read.delim('./DATA/PROCESSED_DATA/p_clinical_data_2genes-classification.tsv', check.names = F), N_ONCOGENIC_ALTERATIONS > 0)
# Twogenes results
twogenes <- filter(readRDS('./DATA/ANALYSIS_DATA/3way/3way_PERM_analysis_mf1-cf10_Tissue-Stage-PM_FDR10.RDS'),
                   SIG_FDR10_3way == F)



### ... Processing clinical data ----
clinical_red <- filter(clinical, os_days < 365.25 * 5)
clinical_long <- pivot_longer(clinical_red,
                              cols = colnames(clinical_red)[!colnames(clinical_red) %in% c('SAMPLE_ID', 'CANC_TYPE', 'CANC_SUBTYPE', 'STAGE_PM',
                                                                                           'sex', 'AGE_AT_SEQUENCING', 'os_days', 'os_status')],
                              names_to = 'Pair',
                              values_to = 'Class') %>% 
  na.omit() %>% 
  arrange(SAMPLE_ID)
clinical_long[,c('GeneA', 'GeneB', 'CNA_type')] <- str_split_fixed(clinical_long$Pair, '_', 4)[,c(1, 2, 4)]
clinical_long$CNA_type <- str_sub(clinical_long$CNA_type, 2, 2)
clinical_long$STAGE_PM <- factor(clinical_long$STAGE_PM, levels = c('Primary', 'Metastasis'))
clinical_long$os_status <- ifelse(clinical_long$os_status == 'dead', 1, 0)
clinical_long$sex <- factor(clinical_long$sex, levels = c('Male', 'Female'))
clinical_long <- clinical_long %>% rename('AGE' = 'AGE_AT_SEQUENCING')
clinical_long$os_days <- clinical_long$os_days / (365.25/12)



### ... Preparing sign information ----
twogenes$Estimate <- ifelse(twogenes$CNA_type == 'Loss', twogenes$Estimate * -1, twogenes$Estimate)
twogenes$Estimate_sign <- ifelse(twogenes$Estimate > 0, 'Pos', 'Neg')
twogenes$CNA_type <- ifelse(twogenes$CNA_type == 'Loss', 'L', 'G')



### ... Adding sign information to clinical data´
clinical_long <- merge(clinical_long,
                       twogenes[c('Tissue', 'Stage', 'CNA_type', 'GeneA', 'GeneB', 'Estimate_sign')],
                       by.x = c('CANC_TYPE', 'STAGE_PM', 'CNA_type', 'GeneA', 'GeneB'),
                       by.y = c('Tissue', 'Stage', 'CNA_type', 'GeneA', 'GeneB'),
                       all.x = T)
clinical_long <- clinical_long %>% 
  mutate('Class' = ifelse(Class == 'NoMutA_NoCNA_NoMutB', 'WT',
                          ifelse(Class == 'MutA_NoCNA_MutB', 'MutA_and_MutB',
                                 ifelse(Class == 'MutA_CNA_NoMutB', '2wayA_NoMutB',
                                        ifelse(Class == 'MutA_NoCNA_NoMutB', 'MutA_only',
                                               ifelse(Class == 'NoMutA_NoCNA_MutB', 'MutB_only',
                                                      ifelse(Class == 'NoMutA_CNA_NoMutB', 'CNA_only',
                                                             ifelse(Class == 'NoMutA_CNA_MutB', 'CNA_MutB', '2wayA_MutB'))))))))
clinical_long$Class <- factor(clinical_long$Class, levels = c('WT', 'MutA_only', 'CNA_only', 'MutB_only',
                                                              'CNA_MutB', 'MutA_and_MutB', '2wayA_NoMutB', '2wayA_MutB'))



### ... Running cox regression ----
# Getting multivariate model results
walds_pvals_3w <- lapply(unique(clinical_long$Pair), function(pair){
  df <- filter(clinical_long, Pair == pair)
  df$CANC_SUBTYPE <- as.factor(df$CANC_SUBTYPE)
  subt_levels <- levels(df$CANC_SUBTYPE)
  coxfit <- coxph(Surv(os_days, os_status) ~ Class + AGE + sex + CANC_SUBTYPE,
                  data = df)
  summary_cox <- summary(coxfit)
  return_df <- data.frame('GeneA' = unique(df$GeneA),
                          'GeneB' = unique(df$GeneB),
                          'Tissue' = unique(df$CANC_TYPE),
                          'Stage' = unique(df$STAGE_PM),
                          'CNA_type' = unique(df$CNA_type),
                          'Sign' = unique(df$Estimate_sign),
                          'CONTROL_SUBT' = subt_levels[1],
                          'Parameter' = names(summary_cox$coefficients[,5]),
                          'Coefficient' = as.numeric(summary_cox$coefficients[,1]),
                          'Hazard.Ratio' = as.numeric(summary_cox$coefficients[,2]),
                          'Z_value' = as.numeric(summary_cox$coefficients[,4]),
                          'Wald.Pval' = as.numeric(summary_cox$coefficients[,5]),
                          'LogRank.Pval' = as.numeric(summary_cox$sctest['pvalue']),
                          'HR.Lower.95' = as.numeric(summary_cox$conf.int[,3]),
                          'HR.Upper.95' = as.numeric(summary_cox$conf.int[,4]))
  return(return_df)
})
walds_pvals_3w <- bind_rows(walds_pvals_3w)
walds_pvals_3w$Parameter <- str_replace_all(walds_pvals_3w$Parameter, c('sex' = '', 'Class' = '', 'CANC_SUBTYPE' = ''))

### ... Saving file ----
if (!file.exists('./DATA/ANALYSIS_DATA/Survival')){
  dir.create('./DATA/ANALYSIS_DATA/Survival')
}
if (!file.exists('./DATA/ANALYSIS_DATA/Survival/3way')){
  dir.create('./DATA/ANALYSIS_DATA/Survival/3way')
}
setwd('./DATA/ANALYSIS_DATA/Survival/3way')
write.table(walds_pvals_3w, './OS_3way.tsv',
            sep = '\t', row.names = F, quote = F)



### ... If we want to compare 2way vs 2wayMutB and 2wayNoMutB ----
# FIRST, we have to create input_2w in the OS_2way.R script previously
input_2way <- readRDS('../2way/INPUT_2way.RDS')
input_2way <- filter(input_2way, Significant == 'sig')
input_combined <- merge(clinical_long,
                        input_2way[c('SAMPLE_ID', 'CANC_TYPE', 'STAGE_PM', 'CNA_type', 'Gene', 'Class')],
                        by.x = c('SAMPLE_ID', 'CANC_TYPE', 'STAGE_PM', 'CNA_type', 'GeneA'),
                        by.y = c('SAMPLE_ID', 'CANC_TYPE', 'STAGE_PM', 'CNA_type', 'Gene'),
                        suffixes = c('_3w', '_2w')) %>% 
  pivot_longer(cols = c('Class_3w', 'Class_2w'),
               names_to = 'Classification',
               values_to = 'Class')
input_combined <- filter(input_combined, !Class %in% c('WT_Both', 'MutA_Only', 'CNA_Only'))
input_combined$Class <- factor(input_combined$Class, levels = c('2way', 'WT', 'MutA_only', 'CNA_only', 'MutB_only',
                                                                'CNA_MutB', 'MutA_and_MutB', '2wayA_NoMutB', '2wayA_MutB'))
# Running cox regression
walds_pvals_BOTH <- lapply(unique(input_combined$Pair), function(pair){
  df <- filter(input_combined, Pair == pair)
  df$CANC_SUBTYPE <- as.factor(df$CANC_SUBTYPE)
  subt_levels <- levels(df$CANC_SUBTYPE)
  coxfit <- coxph(Surv(os_days, os_status) ~ Class + AGE + sex + CANC_SUBTYPE,
                  data = df)
  summary_cox <- summary(coxfit)
  return_df <- data.frame('GeneA' = unique(df$GeneA),
                          'GeneB' = unique(df$GeneB),
                          'Tissue' = unique(df$CANC_TYPE),
                          'Stage' = unique(df$STAGE_PM),
                          'CNA_type' = unique(df$CNA_type),
                          'Sign' = unique(df$Estimate_sign),
                          'CONTROL_SUBT' = subt_levels[1],
                          'Parameter' = names(summary_cox$coefficients[,5]),
                          'Coefficient' = as.numeric(summary_cox$coefficients[,1]),
                          'Hazard.Ratio' = as.numeric(summary_cox$coefficients[,2]),
                          'Z_value' = as.numeric(summary_cox$coefficients[,4]),
                          'Wald.Pval' = as.numeric(summary_cox$coefficients[,5]),
                          'LogRank.Pval' = as.numeric(summary_cox$sctest['pvalue']),
                          'HR.Lower.95' = as.numeric(summary_cox$conf.int[,3]),
                          'HR.Upper.95' = as.numeric(summary_cox$conf.int[,4]))
  return(return_df)
})
walds_pvals_BOTH <- bind_rows(walds_pvals_BOTH)
walds_pvals_BOTH$Parameter <- str_replace_all(walds_pvals_BOTH$Parameter, c('sex' = '', 'Class' = '', 'CANC_SUBTYPE' = ''))
### ... Saving file ----
write.table(walds_pvals_BOTH, './OS_3way_control-2way.tsv',
            sep = '\t', row.names = F, quote = F)

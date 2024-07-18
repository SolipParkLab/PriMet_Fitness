
### ... Libraries ----
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)



### ... Input files ----
clinical <- filter(read.delim('./DATA/PROCESSED_DATA/p_clinical_data_2way-classification.tsv'), N_ONCOGENIC_ALTERATIONS > 0)
# Removing columns we don't need
clinical_red <- filter(clinical, os_days < 365.25 * 5)
clinical_red[c('PATIENT_ID', 'MET_COUNT', 'PRIM_SITE', 'MET_SITE', 'IS_DIST_MET_MAPPED', 'TUMOR_PURITY', 'SAMPLE_COVERAGE',
               'YEARS_BTW_METAS', 'DIAGNOSES', 'N_ONCOGENIC_ALTERATIONS', 'N_clonal.Variants', 'N_subclonal.Variants', 'BIOPSY_LOCATION',
               'TREATMENT', 'met_count', 'met_site_count', 'surgical_procedure_age', 'evidence_of_mets_age', 'seq_report_age',
               'death_age', 'last_contact_age', 'death_or_last_contact_age', 'AGE_AT_EVIDENCE_OF_METS',
               'AGE_AT_DEATH', 'AGE_AT_SURGERY', 'AGE_AT_LAST_CONTACT', 'STAGE_PWOWM')] <- NULL
clinical_long <- pivot_longer(clinical_red,
                              cols = colnames(clinical_red)[!colnames(clinical_red) %in% c('SAMPLE_ID', 'CANC_TYPE', 'CANC_SUBTYPE', 'STAGE_PM',
                                                                                           'sex', 'AGE_AT_SEQUENCING', 'os_days', 'os_status')],
                              names_to = 'Pair',
                              values_to = 'Class') %>% 
  na.omit() %>% 
  arrange(SAMPLE_ID)
clinical_long[,c('Gene', 'CNA_type', 'Significant')] <- str_split_fixed(clinical_long$Pair, '_', 4)[,c(2, 3, 4)]
clinical_long$Pair <- NULL
clinical_long$CNA_type <- str_sub(clinical_long$CNA_type, 2, 2)
clinical_long$STAGE_PM <- factor(clinical_long$STAGE_PM, levels = c('Primary', 'Metastasis'))
clinical_long$os_status <- ifelse(clinical_long$os_status == 'dead', 1, 0)
clinical_long$sex <- factor(clinical_long$sex, levels = c('Male', 'Female'))
clinical_long$Class <- ifelse(clinical_long$Class == 'Mut_Only', 'MutA_Only', clinical_long$Class)



### ... Single genes analysis ----
input_2way <- clinical_long
input_2way$Class <- factor(input_2way$Class, levels = c('WT_Both', 'MutA_Only', 'CNA_Only', '2way'))
input_2way$Label <- paste(input_2way$Gene, input_2way$CANC_TYPE, input_2way$STAGE_PM, input_2way$CNA_type, sep = '_')
input_2way <- input_2way %>% rename('AGE' = 'AGE_AT_SEQUENCING')
input_2way$os_days <- input_2way$os_days / (365.25/12)
##  Wald P-values
pvals_2v2 <- lapply(unique(input_2way$Label), function(label){
  df <- filter(input_2way, Label == label)
  df$CANC_SUBTYPE <- as.factor(df$CANC_SUBTYPE)
  subt_levels <- levels(df$CANC_SUBTYPE)
  coxfit <- coxph(Surv(os_days, os_status) ~ Class + AGE + sex + CANC_SUBTYPE,
                  data = df)
  summary_2w <- summary(coxfit)
  return_df <- data.frame('Tissue' = unique(df$CANC_TYPE),
                          'Stage' = unique(df$STAGE_PM),
                          'Gene' = unique(df$Gene),
                          'CNA_type' = unique(df$CNA_type),
                          'Significant' = unique(df$Significant),
                          'CONTROL_SUBT' = subt_levels[1],
                          'Parameter' = rownames(summary_2w$coefficients),
                          'Coefficient' = as.numeric(summary_2w$coefficients[,1]),
                          'Hazard.Ratio' = as.numeric(summary_2w$coefficients[,2]),
                          'Z_value' = as.numeric(summary_2w$coefficients[,4]),
                          'Wald.Pval' = as.numeric(summary_2w$coefficients[,5]),
                          'LogRank.Pval' = as.numeric(summary_2w$sctest['pvalue']),
                          'HR.Lower.95' = as.numeric(summary_2w$conf.int[,3]),
                          'HR.Upper.95' = as.numeric(summary_2w$conf.int[,4]))
  return(return_df)
})
pvals_2v2_together <- bind_rows(pvals_2v2)
pvals_2v2_together$Parameter <- str_replace_all(pvals_2v2_together$Parameter, c('sex' = '', 'Class' = '', 'CANC_SUBTYPE' = ''))



### ... Saving files ----
if (!file.exists('./DATA/ANALYSIS_DATA/Survival')){
  dir.create('./DATA/ANALYSIS_DATA/Survival')
}
if (!file.exists('./DATA/ANALYSIS_DATA/Survival/2w')){
  dir.create('./DATA/ANALYSIS_DATA/Survival/2w')
}
setwd('./DATA/ANALYSIS_DATA/Survival/2w')
saveRDS(input_2way, './INPUT_2WAY.RDS')
write.table(pvals_2v2_together, './OS_2W.tsv',
            sep = '\t', row.names = F, quote = F)

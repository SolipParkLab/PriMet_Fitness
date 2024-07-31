
### ... Libraries
library(dplyr)
library(purrr)



### ... Input files
# All pairs in 2-way model
two_way <- readRDS('./DATA/ANALYSIS_DATA/2way/2way_PERM_analysis_mf1-cf10_Tissue-Stage-PM.RDS')
two_way <- two_way[c('Tissue', 'Stage', 'CNA_type', 'Gene', 'SIG_FDR10')] %>% 
  arrange(Tissue, Stage, CNA_type)
# Loading binary matrix
binary_matrix <- read.delim('./DATA/Binary_Matrix.tsv', check.names = F)
rownames(binary_matrix) <- binary_matrix$SAMPLE_ID



### ... Subsetting matrix according to model ----
binary_matrix <- binary_matrix %>% 
	mutate(splitmod = paste(CANC_TYPE, STAGE_PM, sep = '.'))
binary_matrix <- binary_matrix[,c("splitmod", colnames(binary_matrix)[grepl('_mutation|_Gain|_Loss', colnames(binary_matrix))])]
ready_bm <- split(binary_matrix[,c(2:ncol(binary_matrix))], binary_matrix$splitmod)



df_classification <- data.frame('Combination' = c('00','01','10','11'),
                                'Classification' = c('WT_Both', 'CNA_Only', 'Mut_Only', '2way'))
### ... Making classification ----
classification_2way <- apply(two_way, 1, function(row){
  df <- ready_bm[[paste0(row['Tissue'],'.',row['Stage'])]][c(paste0(row['Gene'],'_mutation'),
                                                             paste0(row['Gene'],'_',row['CNA_type']))]
  df$SAMPLE_ID <- row.names(df)
  row.names(df) <- NULL
  df$Combination <- paste0(df[,1], df[,2])
  df <- merge(df,
              df_classification,
              by = 'Combination')
  df[paste0(row['Tissue'], '_', row['Gene'] ,'_',ifelse(row['Stage'] == 'Primary', 'P', 'M'), ifelse(row['CNA_type'] == 'Loss', 'L', 'G'), '_',
            ifelse(row['SIG_FDR10'] == 'TRUE', 'sig', 'nosig'))] <- df$Classification
  return(df[,c(4, 6)])
})
classification_2way <- classification_2way %>% reduce(full_join, by = 'SAMPLE_ID')



### ... Adding this classification to the clinical data
full_clinical <- read.delim('./DATA/surv-clinical-data.tsv', sep = '\t')
full_clinical <- merge(full_clinical,
                       classification_2way,
                       by = 'SAMPLE_ID',
                       all = T)
write.table(full_clinical, './DATA/PROCESSED_DATA/p_clinical_data_2way-classification.tsv',
            sep = '\t', row.names = F, quote = F)

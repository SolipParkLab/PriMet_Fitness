
### ... Libraries
library(dplyr)
library(purrr)



### ... Input files
# Binary matrices
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
# Clinical data
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
clinical_data <- filter(clinical_data,N_ONCOGENIC_ALTERATIONS>0)
# Sig pairs in 3way model
two_way <- readRDS('./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_analysis_mf1-cf10_Tissue-Stage-PM.RDS')
two_way <- two_way[c('Tissue', 'Stage', 'CNA_type', 'Gene', 'SIG_FDR10')] %>% 
  arrange(Tissue, Stage, CNA_type)



### ... Subsetting clinical data according to model ----
model_split_data <- c("CANC_TYPE","STAGE_PM")
clinical_data <- clinical_data %>%
  group_by_at(model_split_data) %>%
  mutate(G = cur_group_id(),
         name = paste(CANC_TYPE,STAGE_PM,sep="."))



### ... Subsetting matrixes according to model ----
ready_bm <- unlist(mapply(function(binary_matrix,clinical_table){
  bm <- subset(binary_matrix,rownames(binary_matrix)%in%clinical_table$SAMPLE_ID)
  bm <- merge(bm,clinical_table[c("SAMPLE_ID","name")],by.x="row.names",by.y="SAMPLE_ID")
  rownames(bm) <- bm$Row.names
  bm$Row.names <- NULL
  split_bm <- split(bm[1:(ncol(bm)-1)],bm$name)
  return(split_bm)},
  unname(binary_mats),
  split(clinical_data,clinical_data$CANC_TYPE),
  SIMPLIFY=F),recursive=F)



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
full_clinical <- read.delim('./DATA/PROCESSED_DATA/p_clinical-data.tsv', sep = '\t')
full_clinical <- merge(full_clinical,
                       classification_2way,
                       by = 'SAMPLE_ID',
                       all = T)
write.table(full_clinical, './DATA/PROCESSED_DATA/p_clinical_data_2way-classification.tsv',
            sep = '\t', row.names = F, quote = F)

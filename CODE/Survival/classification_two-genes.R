
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
# Can change SIG_FDR10_3way == T TO SIG_FDR10_3way == F for non significant pairs
two_genes <- filter(readRDS('./DATA/ANALYSIS_DATA/3way/3wTG-PERM_analysis_mf1-cf10_Tissue-Stage-PM_FDR10.RDS'), SIG_FDR10_3way == T)
two_genes <- two_genes[c('Tissue', 'Stage', 'CNA_type', 'Gene_Pair', 'GeneA', 'GeneB')] %>% 
  arrange(Tissue, Stage, CNA_type)



## Subsetting clinical data according to model ----
model_split_data <- c("CANC_TYPE","STAGE_PM")
clinical_data <- clinical_data %>%
  group_by_at(model_split_data) %>%
  mutate(G = cur_group_id(),
         name = paste(CANC_TYPE,STAGE_PM,sep="."))



## Subsetting matrixes according to model ----
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

ready_bm[[paste0('Breast','.','Metastasis')]]

classification_2genes <- apply(two_genes, 1, function(row){
  df_classification <- data.frame('Combination' = c('000','001','010','011','100','101','110','111'),
                                  'Classification' = c('NoMutA_NoCNA_NoMutB', 'NoMutA_NoCNA_MutB', 'NoMutA_CNA_NoMutB', 'NoMutA_CNA_MutB',
                                                       'MutA_NoCNA_NoMutB', 'MutA_NoCNA_MutB', 'MutA_CNA_NoMutB', 'MutA_CNA_MutB'))
  df <- ready_bm[[paste0(row['Tissue'],'.',row['Stage'])]][c(paste0(row['GeneA'],'_mutation'),
                                                             paste0(row['GeneA'],'_',row['CNA_type']),
                                                             paste0(row['GeneB'],'_mutation'))]
  df$SAMPLE_ID <- row.names(df)
  row.names(df) <- NULL
  df$Combination <- paste0(df[,1], df[,2], df[,3])
  df <- merge(df,
              df_classification,
              by = 'Combination')
  df[paste0(row['Gene_Pair'],'_',row['Tissue'],'_',ifelse(row['Stage'] == 'Primary', 'P', 'M'), ifelse(row['CNA_type'] == 'Loss', 'L', 'G'))] <- df$Classification
  return(df[,c(5, 7)])
})

classification_2genes <- classification_2genes %>% reduce(full_join, by = 'SAMPLE_ID')



### ... Adding this classification to the clinical data
full_clinical <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
full_clinical <- merge(full_clinical,
                       classification_2genes,
                       by = 'SAMPLE_ID',
                       all = T)
write.table(full_clinical, './DATA/PROCESSED_DATA/p_clinical_data_2genes-classification_v2.tsv',
            sep = '\t', row.names = F, quote = F)

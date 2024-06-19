
### ... Libraries ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)



### ... Functions/Variables ----
# Changeable
SPLITMOD <- "Subtype-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
FDR_2way <- "10"



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Uploading input files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- filter(read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv", sep = "\t", header=T), N_ONCOGENIC_ALTERATIONS > 0)
clinic <- unique(clinical_data[c('CANC_TYPE','CANC_SUBTYPE')]) %>% 
  rename('Tissue' = 'CANC_TYPE',
         'Subtype' = 'CANC_SUBTYPE')
mod_res <- filter(read.delim('./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t'),
                  SIG_FDR10 == T)
mod_res$Subtype <- NULL
mod_res_subtypes <- merge(mod_res,
                          clinic,
                          by = 'Tissue')



### ... Subsetting clinical data according to model ----
model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PM")

clinical_data <-
  clinical_data %>%
  group_by_at(model_split_data) %>%
  mutate(G = cur_group_id(),
         name = ifelse(SPLITMOD == "Subtype-Stage-PWOWM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PWOWM,sep="."),
                       ifelse(SPLITMOD == "Subtype-Stage-PM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PM,sep="."),
                              ifelse(SPLITMOD == "Tissue-Stage-PWOWM",paste(CANC_TYPE,NA,STAGE_PWOWM,sep="."),
                                     ifelse(SPLITMOD == "Tissue-Stage-PM",paste(CANC_TYPE,NA,STAGE_PM,sep="."),
                                            ifelse(SPLITMOD == "Subtype",paste(CANC_TYPE,CANC_SUBTYPE,NA,sep="."),
                                                   paste(CANC_TYPE,NA,NA,sep=".")))))))



### ... Subsetting matrixes according to model
ready_bm <- unlist(mapply(function(binary_matrix,clinical_table){
  bm <- subset(binary_matrix,rownames(binary_matrix)%in%clinical_table$SAMPLE_ID)
  bm <- merge(bm,clinical_table[c("SAMPLE_ID","name")],by.x="row.names",by.y="SAMPLE_ID")
  rownames(bm) <- bm$Row.names
  bm$Row.names <- NULL
  split_bm <- split(bm[1:(ncol(bm)-1)],bm$name)
  return(split_bm)},
  unname(binary_mats),
  split(clinical_data, clinical_data$CANC_TYPE),
  SIMPLIFY=F),recursive=F)



### ... Compute mut freq and cnv freq for each pair in the 2way model
subtype_freqs <- bind_rows(apply(mod_res_subtypes, 1, function(row){
  mat <- ready_bm[[paste0(row['Tissue'],'.',row['Subtype'],'.',row['Stage'])]]
  Mut_freq <- sum(mat[paste0(row['Gene'],'_mutation')]) / nrow(mat)
  if (row['CNA_type'] == 'Loss') {
    CNA_freq <- sum(mat[paste0(row['Gene'],'_Loss')]) / nrow(mat)
  }
  if (row['CNA_type'] == 'Gain') {
    CNA_freq <- sum(mat[paste0(row['Gene'],'_Gain')]) / nrow(mat)
  }
  return_row <- data.frame('Gene' = row['Gene'],
                           'Tissue' = row['Tissue'],
                           'Subtype' = row['Subtype'],
                           'Stage' = row['Stage'],
                           'CNA_type' = row['CNA_type'],
                           'Mut_freq' = Mut_freq,
                           'CNA_freq' = CNA_freq)
  return(return_row)
}))
# Selecting genes according to our frequency filters
subtype_freqs_totest <- filter(subtype_freqs, (Mut_freq > as.numeric(freqmut_threshold)) & (CNA_freq > as.numeric(freqcnv_threshold)))



### ... Computing model inputs
gene_list <- unique(subtype_freqs_totest[c('Gene', 'Tissue', 'Subtype', 'Stage')])
twogenes_subtype_input <- lapply(names(ready_bm), function(name){
  bm_matrix <- ready_bm[[name]]
  tiss <- str_split_fixed(name, '\\.', 3)[,1] # Tissue
  subtype <- str_split_fixed(name, '\\.', 3)[,2] # Subtype
  stage <- str_split_fixed(name, '\\.', 3)[,3] # Stage
  genes_A <- filter(gene_list, Tissue == tiss & Subtype == subtype & Stage == stage)$Gene # List of genes A
  # Loop for each pemutation matrix
  to_return <- bind_rows(lapply(genes_A, function(gene_A){
    # Creating list of genes B
    genes_B <- str_replace_all(names(bm_matrix)[grep("_mutation",names(bm_matrix))],"_mutation","")
    genes_B <- genes_B[genes_B != gene_A]
    twogenes_list <- bind_rows(lapply(genes_B,function(gene_B){
      ### ... Matrix with only Gene A and Gene B information
      two_genes_bm <- bm_matrix[c(paste0(gene_A,"_mutation"),paste0(gene_A,"_Gain"),paste0(gene_A,"_Loss"),paste0(gene_B,"_mutation"))]
      names(two_genes_bm) <- c("GA_mut","GA_Gain","GA_Loss","GB_mut")
            ### ... Calculating alteration frequencies 
      merged_gene_mat <- two_genes_bm %>% group_by(GB_mut) %>% summarise(Size = n(),
                                                                         FreqMutA = mean(GA_mut),
                                                                         FreqGainA = mean(GA_Gain),
                                                                         FreqLossA = mean(GA_Loss))
      merged_gene_mat <- merge(data.frame(GB_mut=c(0,1)),merged_gene_mat,all.x=T)
      merged_gene_mat[is.na(merged_gene_mat)] <- 0
      muta <- mean(two_genes_bm$GA_mut)
      mutb <- mean(two_genes_bm$GB_mut)
      mutgain <- mean(two_genes_bm$GA_Gain)
      mutloss <- mean(two_genes_bm$GA_Loss)
      ### ... Merging CNA information and turning loss into -1, wt into 0, gain into 1
      two_genes_bm$cnvs <- paste0(two_genes_bm$GA_Loss,two_genes_bm$GA_Gain)
      two_genes_bm$cnvs <- ifelse(two_genes_bm$cnvs=="10",-1,
                                  ifelse(two_genes_bm$cnvs=="01",1,0))
      ### ... Calculating occurrences (mutA*cnv*mutB)
      possib_occur <- c("1-11","0-11","111","011","101","001",
                        "1-10","0-10","110","010","100","000")
      ### ... Count of each occurrence
      N <- unlist(lapply(seq_along(possib_occur),function(x) sum(paste0(two_genes_bm$GA_mut,two_genes_bm$cnvs,two_genes_bm$GB_mut) == possib_occur[x])))
      ### ... Adding pseudocounts
      N <- N + 1
      N_df <- data.frame(matrix(N,
                                ncol=6,
                                nrow=2,
                                byrow=T))
      names(N_df) <-c("MutALoss","NoMutALoss","MutAGain","NoMutAGain","MutAWT","NoMutAWT")
      N_df$GB_mut <- c(1,0)
      input_row <- merge(merged_gene_mat,N_df,all=T)
      input_row <- bind_cols(input_row,
                             data.frame(GeneA = gene_A,
                                        GeneB = gene_B,
                                        Gene_Pair = paste0(gene_A,"_",gene_B),
                                        FreqMutA_wholeGene = muta,
                                        FreqMutB = mutb))
      return(input_row)}))
    return(twogenes_list)
  }))
  if(nrow(to_return) == 0) { return(NULL) } 
  else { return(filter(to_return, FreqMutB >= as.numeric(freqmut_threshold))) }
})
names(twogenes_subtype_input) <- names(ready_bm)
twogenes_subtype_input <- twogenes_subtype_input %>% 
  discard(is.null)



### ... Saving files ----
if (!file.exists('./DATA/MODEL_INPUTS/')){
  dir.create('./DATA/MODEL_INPUTS/')
}
if (!file.exists('./DATA/MODEL_INPUTS/3way__TG')){
  dir.create('./DATA/MODEL_INPUTS/3way__TG')
}
setwd('./DATA/MODEL_INPUTS/3way__TG')
saveRDS(subtype_freqs_totest, sprintf('./3w__TG__sig-2w-pairs-to-test_%s_%s.RDS', FREQ, SPLITMOD))
saveRDS(twogenes_subtype_input, sprintf('./3w__TG__inputs_%s_%s.RDS', FREQ, SPLITMOD))

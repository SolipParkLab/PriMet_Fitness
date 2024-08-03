
### ... Load libraries and other stuff ----
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
source('./CODE/common_reg-model-functions.R')



### ... Parameter definition ----
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"



### ... Uploading input files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/Pos-Specific_oncokb_binary-matrices.RDS")
# Clinical data
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
# Treatment results
treat_res_sig <- filter(read.delim('./DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t'), SIG_FDR10 == T)
treat_res_sig <- unique(treat_res_sig[c('Gene', 'Tissue', 'Stage', 'Treatment')])



### ... Subsetting clinical data according to model ----
clinical_data <- filter(clinical_data,N_ONCOGENIC_ALTERATIONS>0)
model_split_data <- c("CANC_TYPE","STAGE_PM","TREATMENT")
clinical_data <-
  clinical_data %>% 
  group_by_at(model_split_data) %>% 
  mutate(G = cur_group_id(),
         name = paste(CANC_TYPE,NA,STAGE_PM,TREATMENT,sep="."))



### ... Subsetting matrices according to model ----
ready_bm <- unlist(mapply(function(binary_matrix,clinical_table){
  bm <- subset(binary_matrix,rownames(binary_matrix) %in% clinical_table$SAMPLE_ID)
  bm <- merge(bm,
              clinical_table[c("SAMPLE_ID","name")],
              by.x="row.names",
              by.y="SAMPLE_ID")
  rownames(bm) <- bm$Row.names
  bm$Row.names <- NULL
  split_bm <- split(bm[1:(ncol(bm)-1)],bm$name)
  return(split_bm)},
  unname(binary_mats),
  split(clinical_data, clinical_data$CANC_TYPE),
  SIMPLIFY=F),recursive=F)


### ... Function to create position specific inputs (base function needs to be changed) ----
pos_specific_inputs <- function(df, freqmut_threshold, freqcnv_threshold){
  # Changing tibble to df, if needed
  if(is_tibble(df)){df <- as.data.frame(df)}
  # Row number indicates the size (number of pairs)
  size <- nrow(df)
  # Retrieving all gene names (from column names) by eliminating the "_mutation" part
  mut_pos <- str_replace_all(colnames(df)[grep("_mutation",colnames(df))],"_mutation","")
  # Looping over mut_pos vector
  inputs_table <- lapply(mut_pos,function(mut_pos){
    # Defining gene and position
    gene <- str_split_fixed(mut_pos, '_', 3)[,1]
    pos <- str_split_fixed(mut_pos, '_', 3)[,2]
    # Retrieving correspondent columns
    gene_bm <- df[,colnames(df)[grepl(paste0(gene,'_'), colnames(df))]]
    # Merging Loss and Gain counts into "cnvs"
    gene_bm$cnvs <- paste0(gene_bm[,colnames(gene_bm)[grepl('Loss',colnames(gene_bm))]],
                           gene_bm[,colnames(gene_bm)[grepl('Gain',colnames(gene_bm))]])
    # Transforming loss and gain counts to -1, 0 or 1
    gene_bm$cnvs <- ifelse(gene_bm$cnvs=="10",-1,
                           ifelse(gene_bm$cnvs=="01",1,0))
    ### ... Vector with possible occurrences
    possib_occur <- c("0-1","1-1","00","10","01","11")
    ### ... Count of each occurrence
    other_positions <- colnames(gene_bm)[grepl('_mutation', colnames(gene_bm))]
    other_positions <- other_positions[other_positions != paste0(gene,'_',pos,'_mutation')]
    other_positions_df <- gene_bm %>% select(all_of(other_positions))
    gene_bm[,other_positions] <- NULL
    gene_bm$muts <- ifelse(gene_bm[,paste0(gene,'_',pos,'_mutation')] == 1, 1,
                           ifelse(rowSums(other_positions_df) == 0, 0, NA))
    N <- unlist(lapply(1:6,function(x) sum(paste0(gene_bm$muts,gene_bm$cnvs)==possib_occur[x]))) ## THIS STEP IS CRUCIAL
    ### ... Adding pseudocounts
    N <- N + 1
    N_df <- data.frame(t(N))
    names(N_df) <- c("NoMutLoss","MutLoss","NoMutWT","MutWT","NoMutGain","MutGain")
    ### ... Calculating mutation frequency from binary vector
    mutfreq <- sum(gene_bm[,colnames(gene_bm)[grepl('mutation', colnames(gene_bm))]])/nrow(gene_bm)
    ### ... Calculating CNA loss frequency from binary vector
    lossfreq <- sum(gene_bm[,colnames(gene_bm)[grepl('Loss', colnames(gene_bm))]])/nrow(gene_bm)
    ### ... Calculating CNA gain frequency from binary vector
    gainfreq <- sum(gene_bm[,colnames(gene_bm)[grepl('Gain', colnames(gene_bm))]])/nrow(gene_bm)
    input_row <- bind_cols(data.frame(Gene = gene,
                                      Position = pos,
                                      Alteration = paste0(gene,'_',pos),
                                      Size = size,
                                      MutFreq = mutfreq,
                                      LossFreq = lossfreq,
                                      GainFreq = gainfreq,
                                      Mutation_filter = mutfreq >= freqmut_threshold,
                                      Loss_filter = lossfreq >= freqcnv_threshold,
                                      Gain_filter = gainfreq >= freqcnv_threshold),
                           N_df)
    return(input_row)})
  inputs_table <- bind_rows(inputs_table)
  return(inputs_table)
}



### ... Running input with significant genes in 2w__TREATMENT ----
# Subsetting matrixes to reduce time
ready_bm_subset <- lapply(names(ready_bm), function(matname){
  mat <- ready_bm[[matname]]
  tiss <- str_split_fixed(matname, '\\.', 4)[,1]
  stage <- str_split_fixed(matname, '\\.', 4)[,3]
  treatment <- str_split_fixed(matname, '\\.', 4)[,4]

  gene_list <- filter(treat_res_sig, Tissue == tiss & Stage == stage & Treatment == treatment)$Gene
  if(length(gene_list) > 0){
    pattern <- paste0(paste(gene_list, collapse = '_|'),'_')
    return(mat[,colnames(mat)[grepl(pattern, colnames(mat))]])
  } else {return()}
})
names(ready_bm_subset) <- names(ready_bm)
ready_bm_subset <- ready_bm_subset %>%
  discard(is.null)
# Computing inputs
model_inputs_subset <- lapply(ready_bm_subset, pos_specific_inputs, freqmut_threshold, freqcnv_threshold)
model_inputs_df_subset <- map_df(model_inputs_subset, ~as.data.frame(.x), .id="id")
model_inputs_df_subset[c("Tissue","Subtype","Stage","Treatment")] <- str_split_fixed(model_inputs_df_subset$id,"\\.",4)
model_inputs_df_subset$id <- NULL
saveRDS(model_inputs_subset, sprintf('./DATA/GLM_INPUTS/2way_Treatment/2way_Treatment__glm-inputs_Pos-Specific_%s_%s.RDS', 'mf05-cf10', 'Tissue-Stage-PM'))



### ... Computing outputs ----
model_inputs <- readRDS(sprintf('./DATA/GLM_INPUTS/2way_Treatment/2way_Treatment__glm-inputs_Pos-Specific_%s_%s.RDS', 'mf05-cf10', 'Tissue-Stage-PM'))
# Function to create result table for the 2way model 
pos_specific_results <- function(df) {
  if(is.null(df)){return(NULL)}
  ### ... Filtering by mutation frequency threshold
  mut_df <- filter(df,Mutation_filter==TRUE)
  to_return <- list()
  for (cna in c("Loss","Gain")){
    ### ... Filtering by cnv frequency threshold
    cnv_df <- filter(mut_df,get(paste0(cna,"_filter"))==T)
    if (nrow(cnv_df)==0){
      next
    }else{
      ### ... Running regression model over every gene in cnv_df
      regmod_df <- as.data.frame(t(sapply(cnv_df$Alteration, function(pos){
        ### ... Retrieving the correspondent row
        datarow <- cnv_df[cnv_df$Alteration==pos,]
        ### ... Creating input dataframe for the regression model
        if (cna=="Loss"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutLoss,datarow$MutLoss)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=N)
        }else if (cna=="Gain"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutGain,datarow$MutGain)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=N)
        }
        ### ... Running the regression model
        result <- reg_model_2w(dataf)
        return(result)
      })))
      colnames(regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                               "NoMutWT","MutWT","NoMutCNV","MutCNV")
      regmod_df$Alteration <- rownames(regmod_df)
      to_return[[cna]] <- regmod_df
    }
  }
  return(to_return)
}
# Running the output
model_results <- lapply(model_inputs, pos_specific_results)
model_results <- lapply(model_results,function(x){return(map_df(x,~as.data.frame(.x),.id="CNA_type"))})
model_results_df <- map_df(model_results, ~as.data.frame(.x), .id="id")
model_results_df[c("Tissue","Subtype","Stage","Treatment")] <- str_split_fixed(model_results_df$id,"\\.",4)
model_results_df$id <- NULL

model_results_df <- merge(model_results_df,
                          model_inputs_df_subset[c("Gene", "Alteration","Tissue","Subtype","Stage","Treatment","MutFreq","LossFreq","GainFreq","Size")],
                          by = intersect(names(model_results_df),
                                         c("Gene", "Alteration", "Tissue","Subtype","Stage","Treatment","MutFreq","LossFreq","GainFreq")))
saveRDS(model_results_df, sprintf('./DATA/GLM_OUTPUTS/2way_Treatment/2way_Treatment__glm-outputs_Pos-Specific_%s_%s.RDS', 'mf05-cf10', 'Tissue-Stage-PM'))
write.table(model_results_df, sprintf('./DATA/GLM_OUTPUTS/2way_Treatment/2way_Treatment__glm-outputs_Pos-Specific_%s_%s.tsv', 'mf05-cf10', 'Tissue-Stage-PM'),
            sep = '\t', row.names = F, quote = F)


### ... Analysing results
total <- model_results_df
total$Position <- str_split_fixed(total$Alteration,'_',2)[,2]
total$Estimate_plot <- ifelse(total$CNA_type == 'Gain', total$Estimate, total$Estimate * -1) # Correcting models' estimates
# Removing some columns we don't need
total[c('Estimate', 'Subtype', 'Size')] <- NULL
# Saving
write.table(total, sprintf('./DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_analysis_Pos-Specific_%s_%s.tsv', 'mf05-cf10', 'Tissue-Stage-PM'),
            sep = '\t', row.names = F, quote = F)

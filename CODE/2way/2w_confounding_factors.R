
### ... Libraries ----
library(tidyr)
library(reshape2)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
options(dplyr.summarise.inform = FALSE)



### ... Functions/Variables ----
# Not changeable 
SPLITMOD <- "Tissue-Stage-PM"
### ... Changeable
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Uploading files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
clonality <- read.table('./DATA/ANALYSIS_DATA/Sample_Clonality/clonality_results.tsv',
                        sep = '\t', header = T)
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)



### ... Subsetting clinical data according to model ----
clinical_data <- filter(clinical_data,N_ONCOGENIC_ALTERATIONS>0)
if (SPLITMOD == "Tissue"){model_split_data <- "CANC_TYPE"
}else if (SPLITMOD == "Subtype"){model_split_data <- c("CANC_TYPE","CANC_SUBTYPE")
}else if (SPLITMOD == "Tissue-Stage-PM"){model_split_data <- c("CANC_TYPE","STAGE_PM")
}else if (SPLITMOD == "Subtype-Stage-PM"){model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PM")
}else if (SPLITMOD == "Tissue-Stage-PWOWM"){model_split_data <- c("CANC_TYPE","STAGE_PWOWM")
}else if (SPLITMOD == "Subtype-Stage-PWOWM"){model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PWOWM")}

clinical_data <-
  clinical_data %>%
  group_by_at(model_split_data) %>%
  mutate(G = cur_group_id(),
         name = ifelse(SPLITMOD == "Subtype-Stage-PWOWM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PWOWM,sep="."),
                       ifelse(SPLITMOD == "Subtype-Stage-PM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PM,sep="."),
                              ifelse(SPLITMOD == "Tissue-Stage-PWOWM",paste(CANC_TYPE,NA,STAGE_PWOWM,sep="."),
                                     ifelse(SPLITMOD == "Tissue-Stage-PM",paste(CANC_TYPE,NA,STAGE_PM,sep="."),
                                            ifelse(SPLITMOD == "Subtype",paste(CANC_TYPE,CANC_SUBTYPE,NA,sep="."),paste(CANC_TYPE,NA,NA,sep="."))))))
  )


### ... Subsetting matrices according to model ----
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



### ... Calculating Mut/CNVfreq per sample ----
freqs_list <- lapply(ready_bm,function(mat){
  freqs <- apply(mat,1,function(row){
    MutFreq <- sum(as.integer(row[grep("mutation",names(row))]))/length(row[grep("mutation",names(row))])
    GainFreq <- sum(as.integer(row[grep("Gain",names(row))]))/length(row[grep("Gain",names(row))])
    LossFreq <- sum(as.integer(row[grep("Loss",names(row))]))/length(row[grep("Loss",names(row))])
    return(c(MutFreq,GainFreq,LossFreq))
  })
  freqs <- setNames(data.frame(t(freqs)),
                    c("MutFreq","GainFreq","LossFreq"))
  return(freqs)
})
freqs_df <- map_df(freqs_list,~as.data.frame(.x),.id="Tissue")
freqs_df[c("Tissue","Subtype","Stage")] <- str_split_fixed(freqs_df$Tissue,"\\.",3)
confounding_factors_df <- merge(clinical_data[c("SAMPLE_ID", 'CANC_TYPE', "STAGE_PM", 'CANC_SUBTYPE', 'STAGE_PWOWM',
                                                "TUMOR_PURITY","SAMPLE_COVERAGE","AGE_AT_SEQUENCING")],
                                freqs_df[c("Tissue","MutFreq","LossFreq","GainFreq")],
                                by.y="row.names", by.x="SAMPLE_ID", all.y=T)
confounding_factors_df <- merge(confounding_factors_df,
                                clonality[c('Sample', 'N_CLONAL', 'N_SUBCLONAL', 'Clonal_Fraction')],
                                by.x = 'SAMPLE_ID', by.y = 'Sample', all.x = T)



### ... Assigning CF levels ----
a <- lapply(c('TUMOR_PURITY', 'SAMPLE_COVERAGE', 'AGE_AT_SEQUENCING', 'Clonal_Fraction'), function(cf){
  column <- ifelse(cf == 'TUMOR_PURITY', 'TP_Level',
                   ifelse( cf == 'SAMPLE_COVERAGE', 'SC_Level',
                           ifelse(cf == 'AGE_AT_SEQUENCING', 'AGE_Level','CL_Level')))
  sorted <- confounding_factors_df[c('SAMPLE_ID','Tissue', cf)] %>%
    na.omit() %>% 
    group_by(Tissue) %>% 
    arrange(desc(.data[[cf]]),
            bygroup = T)
  # Top 50%
  top <- sorted %>%
    na.omit() %>%
    group_by(Tissue) %>% 
    slice_head(prop = 0.5) %>% 
    mutate(TEMP = 'HIGH')
  # Bottom 50%
  bot <- sorted %>%
    na.omit() %>%
    group_by(Tissue) %>% 
    slice_tail(prop = 0.5) %>% 
    mutate(TEMP = 'LOW')
  # Merging
  top_n_bot <- rbind(top, bot)
  names(top_n_bot)[names(top_n_bot) == 'TEMP'] <- column
  return(top_n_bot)
})



### ... Merging CF levels information ----
confounding_factors_df <- merge(confounding_factors_df,
                                a[[1]][c('SAMPLE_ID','Tissue','TP_Level')],
                                by = c('SAMPLE_ID','Tissue'),
                                all.x = T)
confounding_factors_df <- merge(confounding_factors_df,
                                a[[2]][c('SAMPLE_ID','Tissue','SC_Level')],
                                by = c('SAMPLE_ID','Tissue'),
                                all.x = T)
confounding_factors_df <- merge(confounding_factors_df,
                                a[[3]][c('SAMPLE_ID','Tissue','AGE_Level')],
                                by = c('SAMPLE_ID','Tissue'),
                                all.x = T)
confounding_factors_df <- merge(confounding_factors_df,
                                a[[4]][c('SAMPLE_ID','Tissue','CL_Level')],
                                by = c('SAMPLE_ID','Tissue'),
                                all.x = T)
# Saving table for sample classification according to each confounding factor
write.table(confounding_factors_df[colnames(confounding_factors_df)[colnames(confounding_factors_df) %!in% c('MutFreq', 'GainFreq', 'LossFreq',
                                                                                                             'N_CLONAL', 'N_SUBCLONAL', 'Tissue_color')]],
            './DATA/CONF_FACTORS/Confounding-Factors_Samples-Levels.tsv',
            sep = '\t', row.names = F, quote = F)



### ... Generating inputs for each CF ----
readying_bm <- function(cf){
  ready_datasets <- lapply(ready_bm,function(tiss_matrix){
    # Subsetting matrixes
    bm <- merge(tiss_matrix,na.omit(confounding_factors_df[c("SAMPLE_ID",cf)]),by.x="row.names",by.y="SAMPLE_ID")
    rownames(bm) <- bm$Row.names
    bm$Row.names <- NULL
    split_bm <- split(bm[1:(ncol(bm)-1)],bm[cf])
    return(split_bm)
  })
  ready_datasets <- unlist(ready_datasets,recursive=F)
  return(ready_datasets)
}

generate_inputs <- function(cf){
  bm <- readying_bm(cf)
  model_inputs <- lapply(bm,oncogenic_model_input,freqmut_threshold,freqcnv_threshold)
  model_inputs_df <- map_df(model_inputs, ~as.data.frame(.x), .id="id")
  model_inputs_df[c("Tissue","Subtype","Stage",cf)] <- str_split_fixed(model_inputs_df$id,"\\.",4)
  model_inputs_df$id <- NULL
  return(list(List = model_inputs,
              DF = model_inputs_df))
}
TPinputs <- generate_inputs('TP_Level')
SC_inputs <- generate_inputs('SC_Level')
AGE_inputs <- generate_inputs('AGE_Level')
CL_inputs <- generate_inputs('CL_Level')
# Saving inputs
if (!file.exists('./DATA/GLM_INPUTS/')){
  dir.create('./DATA/GLM_INPUTS/')
}
if (!file.exists('./DATA/GLM_INPUTS/2way__CONF_FACTORS/')){
  dir.create('./DATA/GLM_INPUTS/2way__CONF_FACTORS/')
}
saveRDS(TP_inputs$List, './DATA/GLM_INPUTS/2way__CONF_FACTORS/TP_glm-inputs.RDS')
saveRDS(SC_inputs$List, './DATA/GLM_INPUTS/2way__CONF_FACTORS/SC_glm-inputs.RDS')
saveRDS(AGE_inputs$List, './DATA/GLM_INPUTS/2way__CONF_FACTORS/AGE_glm-inputs.RDS')
saveRDS(CL_inputs$List, './DATA/GLM_INPUTS/2way__CONF_FACTORS/CL_glm-inputs.RDS')



### ... Generating outputs for the models ----
# List of dfs with regression results
generate_outputs <- function(cf_df,cf){
  ### ... Retrieving results
  model_res_list <- lapply(cf_df,oncogenic_retrieve_results)
  ### ... Merging results with input information 
  model_res <- mapply(function(results,inputs){
    cna_list <- lapply(results,function(cna_results){
      cna_results <- merge(cna_results,inputs[c("Gene","MutFreq","LossFreq","GainFreq","Size")],by="Gene")
      return(cna_results)
    })
    cna_df <- map_df(cna_list, ~as.data.frame(.x), .id="CNA_type")
    return(cna_df)
  },model_res_list,cf_df,SIMPLIFY = F)
  model_res <- map_df(model_res, ~as.data.frame(.x), .id="id")
  model_res[c("Tissue","Subtype","Stage",cf)] <- str_split_fixed(model_res$id,"\\.",4)
  model_res$id <- NULL
  return(list(List = model_res_list,
              DF = model_res))
}
TP_outputs <- generate_outputs(TP_inputs$List,'TP_Level')
SC_outputs <- generate_outputs(SC_inputs$List,'SC_Level')
AGE_outputs <- generate_outputs(AGE_inputs$List,'AGE_Level')
CL_outputs <- generate_outputs(CL_inputs$List,'CL_Level')
# Saving outputs
if (!file.exists('./DATA/GLM_OUTPUTS/2way__CONF_FACTORS/')){
  dir.create('./DATA/GLM_OUTPUTS/2way__CONF_FACTORS/')
}
saveRDS(TP_outputs$DF, './DATA/GLM_OUTPUTS/2way__CONF_FACTORS/TP_glm-outputs.RDS')
saveRDS(SC_outputs$DF, './DATA/GLM_OUTPUTS/2way__CONF_FACTORS/SC_glm-outputs.RDS')
saveRDS(AGE_outputs$DF, './DATA/GLM_OUTPUTS/2way__CONF_FACTORS/AGE_glm-outputs.RDS')
saveRDS(CL_outputs$DF, './DATA/GLM_OUTPUTS/2way__CONF_FACTORS/CL_glm-outputs.RDS')

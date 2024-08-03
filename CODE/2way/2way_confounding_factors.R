
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


# INPUT FILES
binary_matrix <- read.delim('./PriMet_Fitness/DATA/Binary_Matrix.tsv', check.names = F)
rownames(binary_matrix) <- binary_matrix$SAMPLE_ID

if (SPLITMOD == "Tissue-Stage-PM"){binary_matrix <- binary_matrix %>% 
  mutate(splitmod = paste(CANC_TYPE, NA, STAGE_PM, sep = '.'))
}else if (SPLITMOD == "Subtype-Stage-PM"){binary_matrix <- binary_matrix %>% 
  mutate(splitmod = paste(CANC_TYPE, CANC_SUBTYPE, STAGE_PM, sep = '.'))
}else if (SPLITMOD == "Tissue-Stage-PWOWM"){binary_matrix <- binary_matrix %>% 
  mutate(splitmod = paste(CANC_TYPE, NA, STAGE_PWOWM, sep = '.'))
}else if (SPLITMOD == "Subtype-Stage-PWOWM"){binary_matrix <- binary_matrix %>% 
  mutate(splitmod = paste(CANC_TYPE, CANC_SUBTYPE, STAGE_PWOWM, sep = '.'))}
binary_matrix <- binary_matrix[,c("splitmod", "AGE_Level", "TP_Level", "SC_Level", "CL_Level",
                                  colnames(binary_matrix)[grepl('_mutation|_Gain|_Loss', colnames(binary_matrix))])]
ready_bm <- split(binary_matrix[,c(2:ncol(binary_matrix))], binary_matrix$splitmod)



### ... Generating inputs for each CF ----
readying_bm <- function(cf){
  ready_datasets <- lapply(ready_bm,function(tiss_matrix){
    # Subsetting matrixes
    split_bm <- split(tiss_matrix[5:(ncol(tiss_matrix))],tiss_matrix[cf])
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
TP_inputs <- generate_inputs('TP_Level')
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

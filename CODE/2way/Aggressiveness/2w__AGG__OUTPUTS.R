
### ... Libraries ----
library(dplyr)
library(purrr)
library(reshape2)
library(stringr)
options(dplyr.summarise.inform = FALSE)



### ... Variables ----
# Not changeable
SPLITMOD <- "Tissue-Stage-PWOWM"
# Changeable
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



## Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



## Uploading files ----
model_inputs <- readRDS(sprintf('./DATA/MODEL_INPUTS/2way__AGG/2w__AGG__inputs_%s_%s.RDS', FREQ, SPLITMOD))



## Generating model outputs ----
model_results <- lapply(model_inputs,oncogenic_retrieve_results)
model_results <- mapply(function(results,inputs){
  cna_list <- lapply(results,function(cna_results){
    cna_results <- merge(cna_results,inputs[c("Gene","MutFreq","LossFreq","GainFreq","Size")],by="Gene")
    return(cna_results)
    return(cna_list)
  })
},model_results,model_inputs,SIMPLIFY = F)
mod_res <- unlist(model_results,recursive=FALSE)
mod_res <- map_df(mod_res, ~as.data.frame(.x), .id="Tissue")
mod_res[c("Tissue","Stage","AGGGROUP","CNA_type")] <- str_split_fixed(mod_res$Tissue,"\\.",4)



### ... Saving files
if (!file.exists('./DATA/MODEL_OUTPUTS/')){
  dir.create('./DATA/MODEL_OUTPUTS/')
}
if (!file.exists('./DATA/MODEL_OUTPUTS/2way__AGG')){
  dir.create('./DATA/MODEL_OUTPUTS/2way__AGG')
}
setwd('./DATA/MODEL_OUTPUTS/2way__AGG')
saveRDS(mod_res, sprintf('./2w__AGG__outputs_%s_%s.RDS', FREQ, SPLITMOD))
write.table(mod_res, sprintf('./2w__AGG__outputs_%s_%s.tsv', FREQ, SPLITMOD),
            sep="\t", quote=FALSE, row.names = FALSE)

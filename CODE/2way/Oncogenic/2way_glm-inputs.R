
### ... Loading libraries ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)



### ... Variables ----
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions
source("./CODE/common_reg-model-functions.R",local=T)



### ... Loading binary matrix ----
binary_matrix <- read.delim('./DATA/Binary_Matrix.tsv', check.names = F)
rownames(binary_matrix) <- binary_matrix$SAMPLE_ID



### ... Subsetting matrixes according to model ----
if (SPLITMOD == "Tissue-Stage-PM"){binary_matrix <- binary_matrix %>% 
	mutate(splitmod = paste(CANC_TYPE, NA, STAGE_PM, sep = '.'))
}else if (SPLITMOD == "Subtype-Stage-PM"){binary_matrix <- binary_matrix %>% 
	mutate(splitmod = paste(CANC_TYPE, CANC_SUBTYPE, STAGE_PM, sep = '.'))
}else if (SPLITMOD == "Tissue-Stage-PWOWM"){binary_matrix <- binary_matrix %>% 
	mutate(splitmod = paste(CANC_TYPE, NA, STAGE_PWOWM, sep = '.'))
}else if (SPLITMOD == "Subtype-Stage-PWOWM"){binary_matrix <- binary_matrix %>% 
	mutate(splitmod = paste(CANC_TYPE, CANC_SUBTYPE, STAGE_PWOWM, sep = '.'))}
binary_matrix <- binary_matrix[,c("splitmod", colnames(binary_matrix)[grepl('_mutation|_Gain|_Loss', colnames(binary_matrix))])]
ready_bm <- split(binary_matrix[,c(2:ncol(binary_matrix))], binary_matrix$splitmod)



### ... Generating model inputs ----
model_inputs <- lapply(ready_bm,oncogenic_model_input,freqmut_threshold,freqcnv_threshold)
model_inputs_df <- map_df(model_inputs, ~as.data.frame(.x), .id="id")
model_inputs_df <- model_inputs_df %>% relocate("id", .after = last_col())
model_inputs_df[c("Tissue","Subtype","Stage")] <- str_split_fixed(model_inputs_df$id,"\\.",3)
if (grepl("Tissue",SPLITMOD)){model_inputs_df$Subtype <- NA}
if (!grepl("Stage",SPLITMOD)){model_inputs_df$Stage <- NA}
model_inputs_df[c("id")] <- NULL



### ... Storing files ----
if (!file.exists('./DATA/GLM_INPUTS/')){
  dir.create('./DATA/GLM_INPUTS/')
}
if (!file.exists('./DATA/GLM_INPUTS/2way')){
  dir.create('./DATA/GLM_INPUTS/2way')
}
setwd('./DATA/GLM_INPUTS/2way')
# Saving
saveRDS(model_inputs,sprintf("2way_glm-inputs_%s_%s.RDS",FREQ,SPLITMOD))
write.table(model_inputs_df,
            sprintf("2way_glm-inputs_%s_%s.tsv",FREQ,SPLITMOD),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

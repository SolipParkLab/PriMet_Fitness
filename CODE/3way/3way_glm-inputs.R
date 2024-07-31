
### ... Loading libraries ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)



### ... Functions/Variables ----
# Changeable
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
FDR_2way <- "10"



### ... Loading functions ----
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



### ... To save files ----
if (!file.exists('./DATA/GLM_INPUTS/')){
  dir.create('./DATA/GLM_INPUTS/')
}
if (!file.exists('./DATA/GLM_INPUTS/3way')){
  dir.create('./DATA/GLM_INPUTS/3way')
}
setwd('./DATA/GLM_INPUTS/3way')



### ... Generating twogenes binary matrixes ----
# If we want to do all 20 matrices in a row, uncomment next two rows. I recommend going 1 by 1.
# saveRDS(model_inputs, sprintf('./3w__TG__inputs_%s_%s.RDS', FREQ, SPLITMOD))
# model_inputs <- lapply(ready_bm,generating_twogenes_binary_mat)

# To go 1 by 1. Save inputs first, then merge into a single file.
args <- "1"
model_inputs <- generating_twogenes_binary_mat(ready_bm[[as.numeric(args[1])]])
if (!file.exists(sprintf('./inputs_by_%s/', SPLITMOD))){
  dir.create(sprintf('./inputs_by_%s/', SPLITMOD))
}
saveRDS(model_inputs, sprintf('./inputs_by_%s/3way_bm_%s.RDS', SPLITMOD, names(ready_bm[as.numeric(args[1])])))
# The loop/function would end here.
files <- list.files(path = sprintf('./inputs_by_%s/', SPLITMOD), pattern="3way_bm_*")
input_files <- lapply(files,function(name){
  return(readRDS(sprintf("./inputs_by_%s/%s",SPLITMOD, name)))
})
inputnames <- str_remove_all(str_remove_all(files,"3way_bm_"),".RDS")
names(input_files) <- inputnames



### ... Saving files ----
saveRDS(input_files, sprintf('./3way__glm-inputs_%s_%s.RDS', FREQ, SPLITMOD))

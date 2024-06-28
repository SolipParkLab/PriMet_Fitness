
### ... Loading libraries ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)
library(eulerr)
library(gridExtra)
library(ggplot2)
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



### ... Uploading files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)



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



### ... To save files ----
if (!file.exists('./DATA/GLM_INPUTS/')){
  dir.create('./DATA/GLM_INPUTS/')
}
if (!file.exists('./DATA/GLM_INPUTS/3way__TG')){
  dir.create('./DATA/GLM_INPUTS/3way__TG')
}
setwd('./DATA/GLM_INPUTS/3way__TG')



### ... Generating twogenes binary matrixes ----
# If we want to do all 20 matrices in a row, uncomment next two rows. I recommend going 1 by 1.
# saveRDS(model_inputs, sprintf('./3w__TG__inputs_%s_%s.RDS', FREQ, SPLITMOD))
# model_inputs <- lapply(ready_bm,generating_twogenes_binary_mat)

# To go 1 by 1. Matrices have to be saved first, then they are merged into a single file.
args <- "1"
model_inputs <- generating_twogenes_binary_mat(ready_bm[[as.numeric(args[1])]])
if (!file.exists(sprintf('./inputs_by_%s/', SPLITMOD))){
  dir.create(sprintf('./inputs_by_%s/', SPLITMOD))
}
saveRDS(model_inputs, sprintf('./inputs_by_%s/3wTG_bm_%s.RDS', SPLITMOD, names(ready_bm[as.numeric(args[1])])))
# The loop/function would end here.
files <- list.files(path = sprintf('./inputs_by_%s/', SPLITMOD), pattern="3wTG_bm_*")
input_files <- lapply(files,function(name){
  return(readRDS(sprintf("./inputs_by_%s/%s",SPLITMOD, name)))
})
inputnames <- str_remove_all(str_remove_all(files,"3wTG_bm_"),".RDS")
names(input_files) <- inputnames



### ... Saving files ----
saveRDS(input_files, sprintf('./3w__TG__glm-inputs_%s_%s.RDS', FREQ, SPLITMOD))

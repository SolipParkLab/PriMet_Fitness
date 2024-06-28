
### ... Libraries ----
library(cluster)
library(plotrix)
library(ggsignif)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(stringr)



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
model_results <- readRDS(sprintf('./DATA/GLM_OUTPUTS/2way__AGG/2w__AGG__outputs_%s_%s.RDS', FREQ, SPLITMOD))



### ... Correcting the estimate ----
total <- model_results
total$Estimate_plot <- ifelse(total$CNA_type=="Gain",total$Estimate,total$Estimate*(-1))
total[c('Estimate', 'Size')] <- NULL



### ... Saving
if (!file.exists('./DATA/ANALYSIS_DATA/')){
  dir.create('./DATA/ANALYSIS_DATA/')
}
if (!file.exists('./DATA/ANALYSIS_DATA/2way__AGG')){
  dir.create('./DATA/ANALYSIS_DATA/2way__AGG')
}
setwd('./DATA/ANALYSIS_DATA/2way__AGG')
saveRDS(total, sprintf('./2wAGG_analysis_%s_%s.RDS', FREQ, SPLITMOD))
write.table(total, sprintf('./2wAGG_analysis_%s_%s.tsv', FREQ, SPLITMOD),
            sep="\t", quote=FALSE, row.names = FALSE)


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
FDR_2way <- 10
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
###### WHICH FDR CORRECTION METHOD ARE WE USING?? ######
FDR_method <- "PERM"



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Loading files ----
# Glm outputs
results_twogenes <- readRDS(sprintf('./DATA/MODEL_OUTPUTS/3way/3way_glm-outputs_%s_%s.RDS', FREQ, SPLITMOD))
# FDR table obtained from permutations
FDR_conversion_table <- read.delim(sprintf("./DATA/ANALYSIS_DATA/3way/3way_FDR-conversion-table_%s_%s.tsv",FREQ,SPLITMOD),
                                   sep="\t", header=T)



### ... Analysing results ----
total <- results_twogenes
total$P_value_temp <- ifelse(total$P_value<.2,round(total$P_value,6),round(total$P_value,3))
total$P_value_temp[total$P_value_temp==0] <- 0.000001
total <- merge(total,FDR_conversion_table,by.x=c("CNA_type","Stage","P_value_temp"),by.y=c("CNA_type","Stage","P_cut"),all.x=T)
names(total)[which(names(total)=="Mean_FDR")] <- "Adj_Pval_3way_PERM"
total$SIG_FDR10_3way <- total$Adj_Pval_3way_PERM <= 0.1
total$Log_Pval <- -log10(total$P_value + 0.0000001)
total$Log_Pval_3way_PERM <- -log10(total$Adj_Pval_3way_PERM + 0.0000001)
total <- bind_cols(total[1:29], as.data.frame(apply(total[30:ncol(total)],2,as.numeric)))
# Correction of 2-way estimates
total$Estim_2w_Plot <- ifelse(total$CNA_type=="Loss",total$Estimate_2w*(-1),total$Estimate_2w)
total$Estim_2w_NoMutB_Plot <- ifelse(total$CNA_type=="Loss",total$Estimate_2w_NoMutB*(-1),total$Estimate_2w_NoMutB)
total$Estim_2w_MutB_Plot <- ifelse(total$CNA_type=="Loss",total$Estimate_2w_MutB*(-1),total$Estimate_2w_MutB)
# Correction of 3-way estimate
total$Estimate <- ifelse(total$CNA_type == 'Loss', total$Estimate * -1, total$Estimate)



### ... Saving files ----
if (!file.exists('./DATA/ANALYSIS_DATA/')){
  dir.create('./DATA/ANALYSIS_DATA/')
}
if (!file.exists('./DATA/ANALYSIS_DATA/3way')){
  dir.create('./DATA/ANALYSIS_DATA/3way')
}
setwd('./DATA/ANALYSIS_DATA/3way')
saveRDS(total, sprintf("3way_%s_analysis_%s_%s_FDR%s.RDS",FDR_method,FREQ,SPLITMOD,FDR_2way))
write.table(total, sprintf("3way_%s_analysis_%s_%s_FDR%s.tsv",FDR_method,FREQ,SPLITMOD,FDR_2way),
            sep="\t", quote=FALSE, row.names = FALSE)

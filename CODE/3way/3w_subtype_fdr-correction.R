
### ... Loading libraries ----
library(stringr)
library(tidyr)
library(tibble)
library(plyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
options(dplyr.summarise.inform = FALSE)



### ... Loading functions and parameters
source("./CODE/common_reg-model-functions.R",local=T)
# Parameters: CHANGEABLE
SPLITMOD <- 'Subtype-Stage-PM'
FREQ <- 'mf1-cf10'



### ... Loading input files ----
twogenes_results <- readRDS(sprintf('./DATA/GLM_OUTPUTS/3way__TG/3w__TG__glm-outputs_%s_%s.RDS', FREQ, SPLITMOD))
FDR_bytiss <- readRDS(sprintf('./DATA/ANALYSIS_DATA/3way__TG/3wTG_FDR-conversion-table_%s_%s.RDS', FREQ, SPLITMOD))



### .... Adjusting FDR by tissue
twogenes_results$P_cut <- ifelse(twogenes_results$P_value < .2,
                                 round_any(twogenes_results$P_value, accuracy = 0.000001, f = ceiling),
                                 round_any(twogenes_results$P_value, accuracy = 0.001, f = ceiling))
twogenes_results <- merge(twogenes_results,
                          FDR_bytiss,
                          by = c('Tissue', 'Stage', 'CNA_type', 'P_cut'))
names(twogenes_results)[which(names(twogenes_results)=="Mean_FDR")] <- "Adj_Pval_3w_PERM_bytiss"
twogenes_results$Adj_Pval_3w_PERM_bytiss <- ifelse(twogenes_results$Adj_Pval_3w_PERM_bytiss < 0, 0, twogenes_results$Adj_Pval_3w_PERM_bytiss)



### ... Adding SIG columns ----
twogenes_results$SIG_FDR10_3way <- twogenes_results$Adj_Pval_3w_PERM_bytiss <= 0.1



### ... Saving files ----
setwd('./DATA/ANALYSIS_DATA/3way__TG/')
write.table(twogenes_results, sprintf('./3wTG_PERM_analysis_%s_%s_FDR10.tsv', FREQ, SPLITMOD),
            sep = '\t', row.names = F, quote = F)
saveRDS(twogenes_results, sprintf('./3wTG_PERM_analysis_%s_%s_FDR10.RDS', FREQ, SPLITMOD))

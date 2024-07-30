
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



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Loading files ----
sig_genes_2way <- read.csv(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_%s_analysis_%s_%s.tsv", 'PERM', FREQ, SPLITMOD),
                           sep="\t", header=T)
sig_genes_2way <- filter(sig_genes_2way,get(paste0("SIG_FDR", FDR_2way))==T)
inputs_twogenes <- readRDS(sprintf('./DATA/GLM_INPUTS/3way/3w__glm-inputs_%s_%s.RDS', FREQ, SPLITMOD))



### ... Running model ----
# Running two-genes model
sig_genes_2way$Model <- paste0(sig_genes_2way$Tissue,".",NA,".",sig_genes_2way$Stage)
sig_genes_2way_list <- split(sig_genes_2way,sig_genes_2way$Model)
### ... Must be true
all(names(sig_genes_2way_list)==names(inputs_twogenes))
twogenes_outputs <- mapply(function(sig_genes_2way,twog_bm){
  twogenes_results <- apply(sig_genes_2way,1,function(row){
    gene <- row["Gene"]
    cna <- row["CNA_type"]
    twogenes_counts <- filter(twog_bm,GeneA==gene)
    twogenes_counts <- filter(twogenes_counts,FreqMutB>=freqmut_threshold)
    results_2way <- twogenes_retrieve_2way_results(twogenes_counts,cna)
    wide_twogenes_inputs <- pivot_wider(twogenes_counts,names_from=GB_mut,values_from=2:11)
    names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),"_1","MutB")
    names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),"_0","NoMutB")
    if (cna=="Loss"){
      notcna <- "Gain"
      dataf <-  data.frame(mutA=c(1,0,1,0,1,0,1,0),
                           cnvA=c(-1,-1,0,0,-1,-1,0,0),
                           mutB=c(1,1,1,1,0,0,0,0))
    }else if (cna=="Gain"){
      notcna <- "Loss"
      dataf <-  data.frame(mutA=c(1,0,1,0,1,0,1,0),
                           cnvA=c(1,1,0,0,1,1,0,0),
                           mutB=c(1,1,1,1,0,0,0,0))
    }
    wide_twogenes_inputs <- wide_twogenes_inputs[-grep(paste0("A",notcna),names(wide_twogenes_inputs))]
    names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),paste0("A",cna),"ACNV")
    wide_twogenes_inputs$CNA_type <- cna
    reg_model <- as.data.frame(t(apply(wide_twogenes_inputs,1,function(datarow){
      N <- c(datarow["MutACNVMutB"],datarow["NoMutACNVMutB"],datarow["MutAWTMutB"],datarow["NoMutAWTMutB"],
             datarow["MutACNVNoMutB"],datarow["NoMutACNVNoMutB"],datarow["MutAWTNoMutB"],datarow["NoMutAWTNoMutB"])
      dataf$N <- as.numeric(N)
      result <- reg_model_3w_twog(dataf)
      return(result)
    },simplify=T)))
    colnames(reg_model) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance")
    wide_twogenes_inputs <- bind_cols(wide_twogenes_inputs,reg_model)
    wide_twogenes_inputs <- bind_cols(wide_twogenes_inputs,results_2way)
    return(wide_twogenes_inputs)
  })
  twogenes_results <- bind_rows(twogenes_results)
  return(twogenes_results)
},
sig_genes_2way_list,
inputs_twogenes,
SIMPLIFY = F)

twogenes_outputs <- map_df(twogenes_outputs,~as.data.frame(.x),.id="Tissue")
twogenes_outputs[c("Tissue","Subtype","Stage")] <- str_split_fixed(twogenes_outputs$Tissue,"\\.",3)

results_2way_all <- setNames(sig_genes_2way[c("Tissue","Gene","Stage", "CNA_type","Estimate","P_value","Adj_Pval_PERM")],
                             c(c("Tissue","GeneA","Stage", "CNA_type","Estimate_2w","P_value_2w","Adj_Pval_2w")))
twogenes_outputs <- merge(twogenes_outputs,results_2way_all)



### ... Saving files ----
if (!file.exists('./DATA/GLM_OUTPUTS/')){
  dir.create('./DATA/GLM_OUTPUTS/')
}
if (!file.exists('./DATA/GLM_OUTPUTS/3way')){
  dir.create('./DATA/GLM_OUTPUTS/3way')
}
setwd('./DATA/GLM_OUTPUTS/3way')
saveRDS(twogenes_outputs, sprintf('./3w__glm-outputs_%s_%s.RDS', FREQ, SPLITMOD))

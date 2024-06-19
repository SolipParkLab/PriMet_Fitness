
### ... Loading libraries ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)




### ... Functions/Variables ----
# Changeable
SPLITMOD <- "Subtype-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
FDR_2way <- "10"



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Input files
model_inputs <- readRDS(sprintf('./DATA/MODEL_INPUTS/3way__TG/3w__TG__inputs_%s_%s.RDS', FREQ, SPLITMOD))
gene_list <- readRDS(sprintf('./DATA/MODEL_INPUTS/3way__TG/3w__TG__sig-2w-pairs-to-test_%s_%s.RDS', FREQ, SPLITMOD)) %>% 
  mutate(model = paste0(Tissue,'.',Subtype,'.',Stage))
gene_list_split <- gene_list[c('Gene', 'Tissue', 'Subtype', 'Stage', 'CNA_type', 'model')]
gene_list_split <- split(gene_list_split, gene_list$model)


### ... Running output
# Function
twogenes_2way_results <- function(df,cna) {
  if(is.null(df)){return(NULL)}
  to_return <- list()
  regmod_df <- sapply(unique(df$Gene_Pair), function(gene_pair){
    datarows <- df[df$Gene_Pair==gene_pair,]
    resultrows <- apply(datarows,1,function(datarow){
      if (cna=="Loss"){
        N <- c(datarow[["NoMutAWT"]],datarow[["MutAWT"]],datarow[["NoMutALoss"]],datarow[["MutALoss"]])
        dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=as.numeric(N))
      }else if (cna=="Gain"){
        N <- c(datarow[["NoMutAWT"]],datarow[["MutAWT"]],datarow[["NoMutAGain"]],datarow[["MutAGain"]])
        dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=as.numeric(N))
      }
      result <- reg_model_2w(dataf)
      return(result)
    })
    resultrow <- as.data.frame(t(c(resultrows[1:4],resultrows[1:4,2])))
    colnames(resultrow) <- c("Estimate_2w_NoMutB","Std_Error_2w_NoMutB","z_value_2w_NoMutB","P_value_2w_NoMutB",
                             "Estimate_2w_MutB","Std_Error_2w_MutB","z_value_2w_MutB","P_value_2w_MutB")
    return(resultrow)})
  regmod_df <- as.data.frame(t(regmod_df))
  return(regmod_df)
}
# Computing the output
# all(names(gene_list_split)==names(model_inputs)) # Must be true
twogenes_outputs <- mapply(function(mod_res_subtype,twog_bm){
  twogenes_results <- apply(mod_res_subtype,1,function(row){
    gene_A <- row["Gene"]
    cna <- row["CNA_type"]
    twogenes_counts <- filter(twog_bm, GeneA == gene_A)
    twogenes_counts <- filter(twogenes_counts,FreqMutB >= freqmut_threshold)
    results_2way <- twogenes_2way_results(twogenes_counts,cna)
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
    if (nrow(wide_twogenes_inputs) == 0) {
      return(NULL)
    } else {
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
      }
    })
  twogenes_results <- bind_rows(twogenes_results)
  return(twogenes_results)
  },
  gene_list_split,
  model_inputs,
  SIMPLIFY = F)

model_outputs <- map_df(twogenes_outputs,~as.data.frame(.x),.id="Tissue")
model_outputs[c("Tissue","Subtype","Stage")] <- str_split_fixed(model_outputs$Tissue,"\\.",3)
model_outputs$Estimate_plot <- ifelse(model_outputs$CNA_type == 'Loss', model_outputs$Estimate * -1, model_outputs$Estimate)

a <- model_outputs
for (name in colnames(a)) {
  if (class(a[,name]) == 'list') {
    a[, name] <- unlist(a[, name])
  } else {
    a[,name] <- a[,name]}}



### ... Saving files ----
if (!file.exists('./DATA/MODEL_OUTPUTS/')){
  dir.create('./DATA/MODEL_OUTPUTS//')
}
if (!file.exists('./DATA/MODEL_OUTPUTS/3way__TG')){
  dir.create('./DATA/MODEL_OUTPUTS/3way__TG')
}
setwd('./DATA/MODEL_OUTPUTS/3way__TG')
saveRDS(a, sprintf('./3w__TG__outputs_%s_%s.RDS', FREQ, SPLITMOD))
write.table(a, sprintf('./3w__TG__outputs_%s_%s.tsv', FREQ, SPLITMOD),
            sep = '\t', row.names = F, quote = F)

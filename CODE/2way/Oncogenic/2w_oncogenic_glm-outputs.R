
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



### ... Variables ----
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Uploading files ----
model_inputs <- readRDS(sprintf("./DATA/MODEL_INPUTS/2way__OG/2w__OG__glm-inputs_%s_%s.RDS",FREQ,SPLITMOD))



### ... Generating model outputs ----
model_results <- lapply(model_inputs,oncogenic_retrieve_results)
model_results <- mapply(function(results,inputs){
  cna_list <- lapply(results,function(cna_results){
    cna_results <- merge(cna_results,inputs[c("Gene","MutFreq","LossFreq","GainFreq","Size")],by="Gene")
    return(cna_results)
  return(cna_list)
    })
},model_results,model_inputs,SIMPLIFY = F)
mod_res <- unlist(model_results,recursive=FALSE)
if (grepl("Stage",SPLITMOD)){mod_res <- unlist(model_results,recursive=F)}
mod_res <- map_df(mod_res, ~as.data.frame(.x), .id="id")
mod_res <- mod_res %>% relocate("id", .after = last_col())
mod_res[c("Tissue","Subtype","Stage","CNA_type")] <- str_split_fixed(mod_res$id,"\\.",4)
if (grepl("Tissue",SPLITMOD)){mod_res$Subtype <- NA}
if (!grepl("Stage",SPLITMOD)){mod_res$Stage <- NA}
mod_res[c("id")] <- NULL



### ... Counting number of tested tissues under each condition (by tissue, by subtype, by tissue and stage, by subtype and stage) ----
ntt_tiss <- melt(table(unique(as.data.frame(cbind(mod_res$Gene,mod_res$Tissue)))$V1))
mod_res <- merge(mod_res,ntt_tiss,by.x="Gene",by.y="Var1",all.x=T)
colnames(mod_res)[which(colnames(mod_res)=="value")] <- "NTT_Tissue"
ntt_subtype <- melt(table(unique(as.data.frame(cbind(mod_res$Gene,paste0(mod_res$Tissue,mod_res$Subtype))))$V1))
mod_res <- merge(mod_res,ntt_subtype,by.x="Gene",by.y="Var1",all.x=T)
colnames(mod_res)[which(colnames(mod_res)=="value")] <- "NTT_Subtype"
ntt_tisstage <- melt(table(unique(as.data.frame(cbind(paste0(mod_res$Gene,"-",mod_res$Stage),mod_res$Tissue)))$V1))
ntt_tisstage[c("Gene","Stage")] <- str_split_fixed(ntt_tisstage$Var1,"-",2)
ntt_tisstage[c("Var1")] <- NULL
mod_res <- merge(mod_res,ntt_tisstage[],by=c("Gene","Stage"),all.x=T)
colnames(mod_res)[which(colnames(mod_res)=="value")] <- "NTT_Tissue-Stage"
ntt_substage <- melt(table(unique(as.data.frame(cbind(paste0(mod_res$Gene,"-",mod_res$Stage),paste0(mod_res$Tissue,mod_res$Subtype))))$V1))
ntt_substage[c("Gene","Stage")] <- str_split_fixed(ntt_substage$Var1,"-",2)
ntt_substage[c("Var1")] <- NULL
mod_res <- merge(mod_res,ntt_substage[],by=c("Gene","Stage"),all.x=T)
colnames(mod_res)[which(colnames(mod_res)=="value")] <- "NTT_Subtype-Stage"



### ... Storing files ----
if (!file.exists('./DATA/MODEL_OUTPUTS/')){
  dir.create('./DATA/MODEL_OUTPUTS/')
}
if (!file.exists('./DATA/MODEL_OUTPUTS/2way__OG')){
  dir.create('./DATA/MODEL_OUTPUTS/2way__OG')
}
setwd('./DATA/MODEL_OUTPUTS/2way__OG')
# Saving
saveRDS(mod_res,sprintf("2w__OG__glm-outputs_%s_%s.RDS",FREQ,SPLITMOD))
write.table(mod_res, sprintf("2w__OG__glm-outputs_%s_%s.tsv",FREQ,SPLITMOD),
            sep="\t", quote=FALSE, row.names = FALSE)

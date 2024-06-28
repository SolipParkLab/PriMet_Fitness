
### ... TREATMENT 2WAY MODEL ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(tibble)
options(dplyr.summarise.inform = FALSE)


### ... functions/Variables ----
### ... Changeable
SPLITMOD <- "Subtype-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)


## Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



## Uploading files ----
model_inputs <- readRDS(sprintf("./DATA/MODEL_INPUTS/2way__T/2w__T__glm-inputs_%s_%s.RDS",FREQ,SPLITMOD))



## Model implementation ----
model_results <- lapply(model_inputs,oncogenic_retrieve_results)
model_results <- lapply(model_results,function(x){return(map_df(x,~as.data.frame(.x),.id="CNA_type"))})
model_results_df <- map_df(model_results, ~as.data.frame(.x), .id="id")
model_results_df[c("Tissue","Subtype","Stage","Treatment")] <- str_split_fixed(model_results_df$id,"\\.",4)
model_results_df$id <- NULL
model_results_df <- merge(model_results_df,
                          model_inputs_df[c("Gene",c("Tissue","Subtype","Stage","Treatment"),"MutFreq","LossFreq","GainFreq","Size")],
                          by=intersect(names(model_results_df),c("Gene",c("Tissue","Subtype","Stage","Treatment"),"MutFreq","LossFreq","GainFreq")))
if(grepl("Tissue",SPLITMOD)){
  ntt_tisstagetreat <- melt(table(unique(as.data.frame(cbind(paste0(model_results_df$Gene,
                                                                    "-",
                                                                    model_results_df$Stage,
                                                                    "-",
                                                                    model_results_df$Treatment),
                                                             model_results_df$Tissue)))$V1))
}else if(grepl("Subtype",SPLITMOD)){
  ntt_tisstagetreat <- melt(table(unique(as.data.frame(cbind(paste0(model_results_df$Gene,
                                                                    "-",
                                                                    model_results_df$Stage,
                                                                    "-",
                                                                    model_results_df$Treatment),
                                                             model_results_df$Tissue,
                                                             model_results_df$Subtype)))$V1))
}
ntt_tisstagetreat[c("Gene","Stage","Treatment")] <- str_split_fixed(ntt_tisstagetreat$Var1,"-",3)
ntt_tisstagetreat$Var1 <- NULL
model_results_df <- merge(model_results_df,ntt_tisstagetreat,by=c("Gene","Stage","Treatment"),all.x=T)
colnames(model_results_df)[which(colnames(model_results_df)=="value")] <- "NTT_Tissue-Stage-Treatment"
model_results_df$Size <- NULL



### ... Storing files ----
if (!file.exists('./DATA/GLM_OUTPUTS/')){
  dir.create('./DATA/GLM_OUTPUTS/')
}
if (!file.exists('./DATA/GLM_OUTPUTS/2way__T')){
  dir.create('./DATA/GLM_OUTPUTS/2way__T')
}
setwd('./DATA/GLM_OUTPUTS/2way__T')
saveRDS(model_results_df,sprintf("2w__T__glm-outputs_%s_%s.RDS",FREQ,SPLITMOD))
write.table(model_results_df,
            sprintf("2w__T__glm-outputs_%s_%s.tsv",FREQ,SPLITMOD),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

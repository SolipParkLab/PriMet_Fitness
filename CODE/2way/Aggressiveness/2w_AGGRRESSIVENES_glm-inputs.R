
### ... Libraries ----
library(dplyr)
library(reshape2) 
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tibble)
library(purrr)
options(dplyr.summarise.inform = FALSE)



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
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep="\t", header=T)
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv", 
                       sep = "\t",header = TRUE, stringsAsFactors = FALSE)



## Preparing Clinical Table ----
clinical_data_filtered <- filter(clinical_data, N_ONCOGENIC_ALTERATIONS != 0)
clinical_data_filtered$TEMP <- ifelse((clinical_data_filtered$MET_COUNT) == 0 & (clinical_data_filtered$STAGE_PWOWM == 'Metastasis'),
                                      'REMOVE', 'KEEP')
clinical_data_filtered <- filter(clinical_data_filtered, TEMP == 'KEEP')
clinical_data_filtered$TEMP <- NULL



### ... Stablishing aggressiveness levels distribution with summary tables and plots ----
agg_distribution <- lapply(binary_mats,function(mat){
  ### ... Removing patients without met_count information
  agg_mat <- merge(mat,clinical_data_filtered[c("SAMPLE_ID")],by.x="row.names",by.y="SAMPLE_ID")
  ### ... Retrieving patient list and sorting by number of metastases
  redmat <- filter(clinical_data_filtered[c("SAMPLE_ID","MET_COUNT","STAGE_PWOWM","CANC_TYPE")],SAMPLE_ID%in%agg_mat$Row.names)
  redmat <- redmat[order(redmat$MET_COUNT),]
  ### ... Levels of aggressiveness
  levels <- c("LOW","HIGH")
  med_pwm <- median(filter(redmat,STAGE_PWOWM=="Primary_W_Metastasis")$MET_COUNT)
  med_m <- median(filter(redmat,STAGE_PWOWM=="Metastasis")$MET_COUNT)
  # Primary with metastasis
  patients_pwom <- filter(redmat,STAGE_PWOWM=="Primary_W_Metastasis") %>%
    mutate(AGGGROUP = ifelse(filter(redmat,STAGE_PWOWM=="Primary_W_Metastasis")$MET_COUNT>0 & filter(redmat,STAGE_PWOWM=="Primary_W_Metastasis")$MET_COUNT<=med_pwm,
                             "LOW",
                             "HIGH"))
  agg_mat_pwom_high <- filter(agg_mat,Row.names%in%filter(patients_pwom,AGGGROUP=="HIGH")$SAMPLE_ID)
  rownames(agg_mat_pwom_high) <- agg_mat_pwom_high$Row.names
  agg_mat_pwom_high$Row.names <- NULL
  agg_mat_pwom_low <- filter(agg_mat,Row.names%in%filter(patients_pwom,AGGGROUP=="LOW")$SAMPLE_ID)
  rownames(agg_mat_pwom_low) <- agg_mat_pwom_low$Row.names
  agg_mat_pwom_low$Row.names <- NULL
  # Metastasis
  patients_m <- filter(redmat,STAGE_PWOWM=="Metastasis") %>%
    mutate(AGGGROUP = ifelse(filter(redmat,STAGE_PWOWM=="Metastasis")$MET_COUNT>0 & filter(redmat,STAGE_PWOWM=="Metastasis")$MET_COUNT<=med_m,
                             "LOW",
                             "HIGH"))
  agg_mat_m_high <- filter(agg_mat,Row.names%in%filter(patients_m,AGGGROUP=="HIGH")$SAMPLE_ID)
  rownames(agg_mat_m_high) <- agg_mat_m_high$Row.names
  agg_mat_m_high$Row.names <- NULL
  agg_mat_m_low <- filter(agg_mat,Row.names%in%filter(patients_m,AGGGROUP=="LOW")$SAMPLE_ID)
  rownames(agg_mat_m_low) <- agg_mat_m_low$Row.names
  agg_mat_m_low$Row.names <- NULL
  #Primary without metastasis
  patients_pwOUTm <- filter(redmat,STAGE_PWOWM=="Primary_WO_Metastasis") %>%
    mutate(AGGGROUP = 'LOW')
  agg_mat_pwOUTm <- filter(agg_mat, Row.names %in% patients_pwOUTm$SAMPLE_ID)
  rownames(agg_mat_pwOUTm) <- agg_mat_pwOUTm$Row.names
  agg_mat_pwOUTm$Row.names <- NULL
  summ <-
    data.frame(Tissue=rep(unique(redmat$CANC_TYPE), 5),
               Stage=c(rep("Primary_W_Metastasis",2), rep("Metastasis",2), 'Primary_WO_Metastasis'),
               GROUP=c(rep(c("LOW","HIGH"),2), 'LOW'),
               From=c(0,med_pwm+0.0001,0,med_m+0.0001,0),
               To=c(med_pwm, max(filter(redmat,STAGE_PWOWM=="Primary_W_Metastasis")$MET_COUNT),
                    med_m, max(filter(redmat,STAGE_PWOWM=="Metastasis")$MET_COUNT), 0.0001))
  # return(list(melt(table(agg_mat$AGGGROUP,agg_mat$STAGE2)),agg_mat,summ,levels))
  return(list(summ,
              list(Primary_W_Metastasis.LOW = agg_mat_pwom_low,
                   Primary_W_Metastasis.HIGH = agg_mat_pwom_high,
                   Metastasis.LOW = agg_mat_m_low,
                   Metastasis.HIGH = agg_mat_m_high,
                   Primary_WO_Metastasis.LOW = agg_mat_pwOUTm),
              bind_rows(patients_pwom, patients_m, patients_pwOUTm)))
})



### ... Updating clinical data table with aggressiveness level information ----
clinical_data_filtered <- merge(clinical_data_filtered,bind_rows(lapply(agg_distribution,function(x) return(x[[3]]))))



## Subsetting matrixes according to model ----
agg_mutcnv_mats <- unlist(lapply(agg_distribution,function(x) return(x[[2]])),recursive=F)



## Generating model inputs ----
model_inputs <- lapply(agg_mutcnv_mats,function(tissue_dataset) return(oncogenic_model_input(tissue_dataset,freqmut_threshold,freqcnv_threshold)))
model_inputs_df <- map_df(model_inputs, ~as.data.frame(.x),.id="Tissue")
model_inputs_df[c("Tissue","Stage","AGGGROUP")] <- str_split_fixed(model_inputs_df$Tissue,"\\.",3)



### ... Saving files ----
if (!file.exists('./DATA/MODEL_INPUTS/')){
  dir.create('./DATA/MODEL_INPUTS/')
}
if (!file.exists('./DATA/MODEL_INPUTS/2way__AGG')){
  dir.create('./DATA/MODEL_INPUTS/2way__AGG')
}
setwd('./DATA/MODEL_INPUTS/2way__AGG')
saveRDS(model_inputs, sprintf('./2w__AGG__inputs_%s_%s.RDS', FREQ, SPLITMOD))
write.table(model_inputs_df, sprintf('./2w__AGG__inputs_%s_%s.tsv', FREQ, SPLITMOD),
            sep="\t", quote=FALSE, row.names = FALSE)

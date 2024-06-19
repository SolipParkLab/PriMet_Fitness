
### ... TREATMENT 2WAY MODEL ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(tibble)
library(gridExtra)
options(dplyr.summarise.inform = FALSE)



### ... functions/Variables ----
### ... Changeable
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
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t",
                            header=T)
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv", 
                       sep = "\t", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
color_codes <- read.csv("./DATA/PROCESSED_DATA/p_color-codes_cancer-types.tsv",
                        sep="\t",
                        header=T)



### ... Subsetting clinical data according to model ----
clinical_data <- filter(clinical_data, N_ONCOGENIC_ALTERATIONS > 0)
if (SPLITMOD == "Tissue-Stage-PM"){model_split_data <- c("CANC_TYPE","STAGE_PM","TREATMENT")
}else if (SPLITMOD == "Subtype-Stage-PM"){model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PM","TREATMENT")
}else if (SPLITMOD == "Tissue-Stage-PWOWM"){model_split_data <- c("CANC_TYPE","STAGE_PWOWM","TREATMENT")
}else if (SPLITMOD == "Subtype-Stage-PWOWM"){model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PWOWM","TREATMENT")}

clinical_data <-
  clinical_data %>% 
  group_by_at(model_split_data) %>% 
  mutate(G = cur_group_id(),
         name = ifelse(SPLITMOD == "Subtype-Stage-PWOWM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PWOWM,TREATMENT,sep="."),
                       ifelse(SPLITMOD == "Subtype-Stage-PM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PM,TREATMENT,sep="."),
                              ifelse(SPLITMOD == "Tissue-Stage-PWOWM",paste(CANC_TYPE,NA,STAGE_PWOWM,TREATMENT,sep="."),paste(CANC_TYPE,NA,STAGE_PM,TREATMENT,sep=".")))))



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



### ...Creating model inputs list and df ----
model_inputs <- lapply(ready_bm,oncogenic_model_input,freqmut_threshold,freqcnv_threshold)
model_inputs_df <- map_df(model_inputs, ~as.data.frame(.x), .id="id")
model_inputs_df[c("Tissue","Subtype","Stage","Treatment")] <- str_split_fixed(model_inputs_df$id,"\\.",4)
if (!grepl("Stage",SPLITMOD)){model_inputs_df$Stage <- NA}
model_inputs_df$id <- NULL



### ... Storing files ----
if (!file.exists('./DATA/MODEL_INPUTS/')){
  dir.create('./DATA/MODEL_INPUTS/')
}
if (!file.exists('./DATA/MODEL_INPUTS/2way__T')){
  dir.create('./DATA/MODEL_INPUTS/2way__T')
}
setwd('./DATA/MODEL_INPUTS/2way__T')
saveRDS(model_inputs,sprintf("2w__T__inputs_%s_%s.RDS",FREQ,SPLITMOD))
write.table(model_inputs_df,
            sprintf("2w__T__inputs_%s_%s.tsv",FREQ,SPLITMOD),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

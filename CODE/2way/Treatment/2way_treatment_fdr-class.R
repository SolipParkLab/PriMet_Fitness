
### ... TREATMENT 2WAY MODEL ----
library(stringr)
library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)
library(tibble)
options(dplyr.summarise.inform = FALSE)



### ... functions/Variables ----
### ... Changeable
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
number_tested_tissues_threshold <- 2
FDR_2way <- "10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Loading files ----
# Glm outputs
model_results_df <- readRDS(sprintf("./DATA/GLM_OUTPUTS/2way_Treatment/2way_Treatment_glm-outputs_%s_%s.RDS", FREQ, SPLITMOD))
# Cancer gene list with gene function in tumors
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# FDR conversion table obtained from permutations
FDR_conversion_table <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_FDR-conversion-table_%s_%s.tsv",FREQ,SPLITMOD),
                                   sep="\t", header=T)



### ... Setting repository to store files ----
if (!file.exists('./DATA/ANALYSIS_DATA/')){
  dir.create('./DATA/ANALYSIS_DATA/')
}
if (!file.exists('./DATA/ANALYSIS_DATA/2way_Treatment')){
  dir.create('./DATA/ANALYSIS_DATA/2way_Treatment')
}
setwd('./DATA/ANALYSIS_DATA/2way_Treatment')



## Analysing results ----
### ... FDR correction to assess significance
total <- model_results_df
if(grepl("Stage",SPLITMOD)){
  totalcna <- split(total,list(total$CNA_type,total$Stage,total$Treatment))
}else{
  totalcna <- split(total,list(total$CNA_type,total$Treatment))
}
totalcna <- lapply(totalcna, function(total){
  total$P_value_temp <- ifelse(total$P_value<.2,round(total$P_value,6),round(total$P_value,3))
  total$P_value_temp[total$P_value_temp==0] <- 0.000001
  total <- merge(total,FDR_conversion_table,by.x=c("CNA_type","Stage","Treatment","P_value_temp"),by.y=c("CNA_type","Stage","Treatment","P_cut"),all.x=T)
  names(total)[which(names(total)=="Mean_FDR")] <- "Adj_Pval_PERM"
  total$SIG_FDR10 <- total$Adj_Pval_PERM <= 0.1 & ifelse(total$CNA_type=="Gain",total$Estimate>0,total$Estimate<0)
  return(total)
})
total <- bind_rows(totalcna)
# NTT filtering to assess class
elim <- filter(total,`NTT_Tissue-Stage-Treatment` < number_tested_tissues_threshold)
total <- filter(total,`NTT_Tissue-Stage-Treatment` >= number_tested_tissues_threshold)
fdrs <- c("10")
for (fdr_index in match(c("SIG_FDR10"),names(total))){
    total[,ncol(total)+1] <- ifelse(total[,fdr_index]==T & total$CNA_type=="Gain","C3",
                                  ifelse(total[,fdr_index]==T & total$CNA_type=="Loss","C2","C1"))
}
colnames(total)[(ncol(total))] <- c("PC_FDR10")
total <- total %>% group_by(Gene,Stage,Treatment) %>% 
  mutate(!!str_glue("Class_FDR{fdrs[1]}") := ifelse("C3" %in% PC_FDR10 & "C2" %in% PC_FDR10, "C4",
                                                    ifelse("C3" %in% PC_FDR10 & "C2" %!in% PC_FDR10,"C3",
                                                           ifelse("C3" %!in% PC_FDR10 & "C2" %in% PC_FDR10,"C2","C1"))))
total <- as.data.frame(total)
total[,match(c("PC_FDR10"),names(total))] <- NULL
total <- bind_rows(total,elim)
total <- merge(total,cancgenedf)
total$Estimate_plot <- ifelse(total$CNA_type=="Gain",total$Estimate,total$Estimate*(-1))
total[c('Estimate', 'P_value_temp', 'Size')] <- NULL
# Saving table
write.table(total,
            sprintf("2way_Treatment_%s_analysis_%s_%s.tsv",'PERM',FREQ,SPLITMOD),
            sep="\t",
            row.names = FALSE)

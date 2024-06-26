
### ... Loading libraries ----
library(dplyr)
library(stringr)
library(tibble)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(ppcor)
library(corrr)
library(ggcorrplot)
library(factoextra)
library(gridExtra)
library(rafalib)
library(plotrix)
library(ggalluvial)
library(purrr)




### ... Variables ----
number_tested_tissues_threshold <- 2
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FDR_2way <- "10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Uploading files ----
model_results <- readRDS(sprintf("./DATA/MODEL_OUTPUTS/2way__OG/2w__OG__glm-outputs_%s_%s.RDS",FREQ,SPLITMOD))
model_inputs <- read.delim(sprintf("./DATA/MODEL_INPUTS/2way__OG/2w__OG__glm-inputs_%s_%s.tsv",FREQ,SPLITMOD),
                           sep="\t",
                           header=T)
color_codes <- read.csv("./DATA/PROCESSED_DATA/p_color-codes_cancer-types.tsv",
                        sep="\t",
                        header=T)
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
                       sep = "\t",
                       header = TRUE, 
                       stringsAsFactors = FALSE)
clonal_functional_clinical_maf <- read.delim("./DATA/PROCESSED_DATA/p_clonal-MAF.tsv",
                                             sep="\t",
                                             header=T)
FDR_conversion_table <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_FDR-conversion-table_%s_%s.tsv",FREQ,SPLITMOD),
                                   sep="\t",
                                   header=T)



### ... Setting repository to store files ----
if (!file.exists('./DATA/ANALYSIS_DATA/')){
  dir.create('./DATA/ANALYSIS_DATA/')
}
if (!file.exists('./DATA/ANALYSIS_DATA/2way__OG')){
  dir.create('./DATA/ANALYSIS_DATA/2way__OG')
}
setwd('./DATA/ANALYSIS_DATA/2way__OG')



### ... Analysing results ----
# FDR correction to assess significance
if(grepl("Stage",SPLITMOD)){
  totalcna <- split(model_results,list(model_results$CNA_type,model_results$Stage))
}else{
  totalcna <- split(model_results,model_results$CNA_type)
}
if(grepl("Tissue",SPLITMOD)){
  totalcna <- lapply(totalcna, function(total){
    total$P_value_temp <- ifelse(total$P_value<.2, round(total$P_value,6), round(total$P_value,3))
    total$P_value_temp[total$P_value_temp==0] <- 0.000001
    total <- merge(total,FDR_conversion_table,by.x=c("CNA_type","Stage","P_value_temp"),by.y=c("CNA_type","Stage","P_cut"),all.x=T)
    names(total)[which(names(total)=="Mean_FDR")] <- "Adj_Pval_PERM"
    ### ... Defining significants (pval + estimate (positive/negative))
    ### ... Assessing significance according to PERMUTATION-Adjusted FDR
    total$SIG_FDR10 <- total$Adj_Pval_PERM <= 0.1 & ifelse(total$CNA_type=="Gain",total$Estimate>0,total$Estimate<0)
    return(total)
    })
}
total <- bind_rows(totalcna)
if(grepl("Subtype",SPLITMOD)){
  ### ... FDR correction to assess significance BY TISSUE
  totalcna <- split(total,list(total$CNA_type,total$Stage,total$Tissue))
  totalcna <- lapply(totalcna, function(total){
    total$P_value_temp <- ifelse(total$P_value<.2,round(total$P_value,6),round(total$P_value,3))
    total$P_value_temp[total$P_value_temp==0] <- 0.000001
    total <- merge(total,FDR_conversion_table,by.x=c("CNA_type","Stage","Tissue","P_value_temp"),by.y=c("CNA_type","Stage","Tissue","P_cut"),all.x=T)
    names(total)[which(names(total)=="Mean_FDR")] <- "Adj_Pval_PERM"
  total$SIG_FDR10 <- total$Adj_Pval_PERM <= 0.1 & ifelse(total$CNA_type=="Gain",total$Estimate>0,total$Estimate<0)
  return(total)
  })
  total <- bind_rows(totalcna)
}
### ... NTT (number of tested tissues) filtering to assess class
if(grepl("Stage",SPLITMOD)){
  elim <- filter(total,get(paste0("NTT_",str_split(SPLITMOD,"-")[[1]][1],"-Stage"))<as.numeric(number_tested_tissues_threshold))
  total <- filter(total,get(paste0("NTT_",str_split(SPLITMOD,"-")[[1]][1],"-Stage"))>=as.numeric(number_tested_tissues_threshold))
} else{
  elim <- filter(total,get(paste0("NTT_",str_split(SPLITMOD,"-")[[1]][1]))<as.numeric(number_tested_tissues_threshold))
  total <- filter(total,get(paste0("NTT_",str_split(SPLITMOD,"-")[[1]][1]))>=as.numeric(number_tested_tissues_threshold))
}
fdrs <- c("10")
### ... For Tissue-Stage-PM model
### ... C1 (no interaction), C2 (Loss), C3 (Gain).
if(grepl("Tissue",SPLITMOD)){
  for (fdr_index in match(c("SIG_FDR10"),names(total))){
    total[,ncol(total)+1] <- ifelse(total[,fdr_index]==T & total$CNA_type=="Gain","C3",
                                    ifelse(total[,fdr_index]==T & total$CNA_type=="Loss","C2","C1"))
    }
colnames(total)[(ncol(total))] <- c("PC_FDR10")
### ... Changing to C4 if it were the case, revising that the other classes are correct
total <- total %>% group_by(Gene,Stage) %>% mutate(!!str_glue("Class_FDR{fdrs[1]}") := 
                                                     ifelse("C3" %in% PC_FDR10 & "C2" %in% PC_FDR10, "C4",
                                                            ifelse("C3" %in% PC_FDR10 & "C2" %!in% PC_FDR10,"C3",
                                                                   ifelse("C3" %!in% PC_FDR10 & "C2" %in% PC_FDR10,"C2","C1"))))
total <- as.data.frame(total)
total[,match(c("PC_FDR10"),names(total))] <- NULL
}
### ... Same for subtype
if(grepl("Subtype",SPLITMOD)){
  for (fdr_index in match(c("SIG_FDR10"),names(total))){
    total[,ncol(total)+1] <- ifelse(total[,fdr_index]==T & total$CNA_type=="Gain","C3",
                                    ifelse(total[,fdr_index]==T & total$CNA_type=="Loss","C2","C1"))
    }
  colnames(total)[(ncol(total))] <- c("PC_FDR10")
  total <- total %>% group_by(Gene,Stage) %>% mutate(!!str_glue("Class_FDR{fdrs[1]}") := 
                                                       ifelse("C3" %in% PC_FDR10 & "C2" %in% PC_FDR10, "C4",
                                                              ifelse("C3" %in% PC_FDR10 & "C2" %!in% PC_FDR10,"C3",
                                                                     ifelse("C3" %!in% PC_FDR10 & "C2" %in% PC_FDR10,"C2","C1"))))
  total <- as.data.frame(total)
  total[,match(c("PC_FDR10"),names(total))] <- NULL
  }
total <- bind_rows(total,elim)
total$Estimate_plot <- ifelse(total$CNA_type=="Gain",total$Estimate,total$Estimate*(-1))
total$Log_Pval_max10 <- ifelse(-log(total$P_value,10)>10,10,-log(total$P_value,10))
total <- merge(total,cancgenedf)
if(grepl("Tissue",SPLITMOD)){
  total <- merge(total,unique(color_codes[c("Tissue","Tissue_color")]))
  total <- total %>%
    group_by(Tissue,CNA_type,Stage) %>%
    mutate(N_pairs=n())
}else if(grepl("Subtype",SPLITMOD)){
  total <- merge(total,unique(color_codes[c("Subtype","Subtype_color")]))
  total <- total %>%
    group_by(Tissue,Subtype,CNA_type,Stage) %>%
    mutate(N_pairs=n())
}
# Saving file
total[c('P_value_temp', 'Estimate', 'Size', 'NTT_Tissue', 'NTT_Subtype', 'Tissue_color', 'N_pairs')] <- NULL
write.table(total, sprintf("./2wOG_%s_analysis_%s_%s.tsv",'PERM',FREQ,SPLITMOD),
            sep="\t", quote=FALSE, row.names = FALSE)

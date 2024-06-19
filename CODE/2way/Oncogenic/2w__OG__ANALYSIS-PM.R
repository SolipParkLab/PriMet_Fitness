
### ... Variables ----
SPLITMOD <- 'Tissue-Stage-PM'
FREQ <- 'mf1-cf10'



### ... Libraries ----
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)



### ... Input files ----
total <- read.delim('./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')



###############################################################################
# We're going to analyse class perturbation across stages
FDR_2way <- "10"
total_significants <- total[c("Gene","Stage","Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                              "NoMutWT","MutWT","NoMutCNV","MutCNV","MutFreq","LossFreq","GainFreq","Size","Function",
                              "Tissue","Subtype","CNA_type","NTT_Tissue","NTT_Subtype","NTT_Tissue.Stage","NTT_Subtype.Stage",
                              "Log_Pval_max10","Adj_Pval_PERM",paste0("SIG_FDR",FDR_2way),paste0("Class_FDR",FDR_2way),paste0("C4_FDR",FDR_2way),
                              "Estimate_plot","Pair")]

if(grepl("Tissue",SPLITMOD)){
  names(total_significants)[which(names(total_significants)=="NTT_Tissue.Stage")] <- "Tested_cancer_types" 
}else if(grepl("Subtype",SPLITMOD)){
  names(total_significants)[which(names(total_significants)=="NTT_Subtype.Stage")] <- "Tested_cancer_types" 
}
# Only leaving pairs tested in => number_tested_tissues_threshold -> 2
model_df <- filter(total_significants,Tested_cancer_types>=2)
names(model_df)[which(grepl("FDR",colnames(model_df)))] <- c("SIG","Class","C4")
# Creating class classification tables and plots
class_df <- pivot_wider(unique(model_df[c("Gene","Stage","Class","Function")]),names_from=Stage,values_from=Class)
class_df$Merged <- paste(class_df$Primary,class_df$Metastasis)
class_df$Behaviour <- ifelse(grepl("NA",class_df$Merged),NA,
                             ifelse(class_df$Merged%in%c("C1 C1","C2 C2","C3 C3","C4 C4"),"Consistent","Perturbed"))
write.table(class_df, sprintf('./2wOG_%s_gene-perturbation_%s_%s.tsv', 'PERM', FREQ, SPLITMOD),
            sep = '\t', quote = F, row.names = F)



## Odds ratio calculation ----
### ... All consistent / perturbed genes with all their interactions
behaviour_model_df <- merge(model_df,class_df[c("Gene","Merged","Behaviour")])
behaviour_model_df <- behaviour_model_df[!is.na(behaviour_model_df$Behaviour),]
filtered_clonal_oncogenic_maf <- merge(unique(behaviour_model_df[c("Gene","Tissue","Stage")]),
                                       clonal_functional_clinical_maf[c("Hugo_Symbol","CANC_TYPE","STAGE_PM","clonality")],
                                       by.x=c("Gene","Tissue","Stage"),
                                       by.y=c("Hugo_Symbol","CANC_TYPE","STAGE_PM"))
clonal_stage_counts <-
  filtered_clonal_oncogenic_maf %>% 
  group_by(Gene) %>% 
  summarise(CP=sum(Stage=="Primary"&clonality=="CLONAL")+.5,
            CM=sum(Stage=="Metastasis"&clonality=="CLONAL")+.5,
            SP=sum(Stage=="Primary"&clonality=="SUBCLONAL")+.5,
            SM=sum(Stage=="Metastasis"&clonality=="SUBCLONAL")+.5)

### ... Only C2 / C3 genes significant tissues, C1 should stay the same 
only_tissue_significants <- behaviour_model_df[c("Gene","Tissue","Stage","SIG","Class")]
only_tissue_significants <- only_tissue_significants %>% group_by(Gene) %>% mutate(Keep = ifelse(Class%in%c("C2","C3","C4")&SIG==F,F,T))
only_tissue_significants <- filter(only_tissue_significants,Keep==T)

more_filtered_clonal_oncogenic_maf <- merge(filtered_clonal_oncogenic_maf,
                                            unique(only_tissue_significants[c("Gene","Tissue","Stage")]),
                                            by=c("Gene","Tissue","Stage"))
clonal_stage_counts_sigtissue <-
  more_filtered_clonal_oncogenic_maf %>% 
  group_by(Gene) %>% 
  summarise(CP=sum(Stage=="Primary"&clonality=="CLONAL")+.5,
            CM=sum(Stage=="Metastasis"&clonality=="CLONAL")+.5,
            SP=sum(Stage=="Primary"&clonality=="SUBCLONAL")+.5,
            SM=sum(Stage=="Metastasis"&clonality=="SUBCLONAL")+.5)
clonal_stage_counts_sigtissue <- merge(clonal_stage_counts_sigtissue,unique(behaviour_model_df[c("Gene","Merged","Behaviour")]))

### ... If we want to only include as C1 the tissue where the c2 or c3 is significant, 
only_sig_tissues_significants <- only_tissue_significants %>% group_by(Gene,Tissue) %>%
  mutate(Lab = ifelse(length(unique(Stage))==2,"Y","N"))
only_sig_tissues_significants <- only_sig_tissues_significants %>% group_by(Gene) %>%
  mutate(Keep = ifelse(length(unique(Lab))==2,T,F))
only_sig_tissues_significants <- filter(only_sig_tissues_significants,Keep==T)
more_sig_filtered_clonal_oncogenic_maf <- merge(filtered_clonal_oncogenic_maf,
                                                unique(only_sig_tissues_significants[c("Gene","Tissue","Stage")]),
                                                by=c("Gene","Tissue","Stage"))
clonal_stage_counts_sigtissue_sig <-
  more_sig_filtered_clonal_oncogenic_maf %>% 
  group_by(Gene) %>% 
  summarise(CP=sum(Stage=="Primary"&clonality=="CLONAL")+.5,
            CM=sum(Stage=="Metastasis"&clonality=="CLONAL")+.5,
            SP=sum(Stage=="Primary"&clonality=="SUBCLONAL")+.5,
            SM=sum(Stage=="Metastasis"&clonality=="SUBCLONAL")+.5)
clonal_stage_counts_sigtissue_sig <- merge(clonal_stage_counts_sigtissue_sig,
                                           unique(behaviour_model_df[c("Gene","Merged","Behaviour")]))



### ... Taking into account the location of the metastatic biopsy ---- 
clonal_functional_clinical_maf$Location_plot <- ifelse(clonal_functional_clinical_maf$BIOPSY_LOCATION=="Unspecified"|clonal_functional_clinical_maf$BIOPSY_LOCATION=="Unknown primary"|clonal_functional_clinical_maf$BIOPSY_LOCATION=="Other",
                                                       "Unknown",
                                                       clonal_functional_clinical_maf$BIOPSY_LOCATION)
clonal_functional_clinical_maf$Location_plot <- ifelse(clonal_functional_clinical_maf$Location_plot=="Lypmh",
                                                       "Lymph",
                                                       clonal_functional_clinical_maf$Location_plot)
clonal_counts_location_persample <- filter(clonal_functional_clinical_maf,STAGE_PM!="Primary") %>% 
  group_by(Tumor_Sample_Barcode,CANC_TYPE,Location_plot) %>% 
  summarise(Clonal_Count=sum(clonality=="CLONAL"),
            Subclonal_Count=sum(clonality=="SUBCLONAL"),
            Clonal_Fraction = Clonal_Count/sum(Clonal_Count,Subclonal_Count))

write.table(clonal_counts_location_persample, 
            sprintf("2wOG_%s_Biopsy-Location-CF_%s_%s.tsv",'PERM',FREQ,SPLITMOD),
            sep="\t", row.names=F, quote=F)
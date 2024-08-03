
### ... Loading libraries ----
library(dplyr)
library(readxl)
library(stringr)
library(reshape2)
library(tidyr)
library(tibble)
library(purrr)



### ... Function to negate "%in%" ----
"%!in%" <- Negate("%in%")



### ... OncoKB cancer gene list ----
cancerGeneList <- read.csv("./DATA/RAW_DATA/cancerGeneList.tsv", 
                           sep = "\t",  header = TRUE, stringsAsFactors = FALSE)
cancgenedf <- filter(cancerGeneList, OncoKB.Annotated == "Yes")[c("Hugo.Symbol","Is.Oncogene","Is.Tumor.Suppressor.Gene")]
colnames(cancgenedf) <- c("Gene","Is.OG","Is.TSG")
cancgenedf$Function <- ifelse(cancgenedf$Is.OG=="Yes" & cancgenedf$Is.TSG=="No","OG",
                              ifelse(cancgenedf$Is.OG=="Yes" & cancgenedf$Is.TSG=="Yes","DFG",
                                     ifelse(cancgenedf$Is.OG=="No" & cancgenedf$Is.TSG=="Yes","TSG",NA)))
cancgenedf <- na.omit(cancgenedf)
write.table(cancgenedf[c(1,4)], "./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
            sep = "\t", row.names = F, quote = F)



### ... Gene Name Conversion Table Preparation ----
# cBioPortal data annotated with oncoKB
mut_annot_data <- read.delim("./DATA/PROCESSED_DATA/p_OKB-annotated_MAF_data_mutations.oncokb.txt",
                             sep="\t")
mut_annot_reduced <- unique(mut_annot_data[c("Hugo_Symbol","GENE_IN_ONCOKB")])
mut_annot_reduced$Is.in.cancerGeneList.AsAliases <- ifelse(mut_annot_reduced$Hugo_Symbol %in% unique(sort(unlist(str_split(cancerGeneList$Gene.Aliases,", ")))),
                                                           T,
                                                           F)
mut_annot_reduced$New.Hugo.Symbol <- unlist(apply(mut_annot_reduced,1,function(row){
  if (row["Is.in.cancerGeneList.AsAliases"]==T){
    newsymb <- cancerGeneList[grep(paste0("\\b",row["Hugo_Symbol"],"\\b"),cancerGeneList$Gene.Aliases),]$Hugo.Symbol
    return(newsymb)
  }else{
    return(NA)
  }
}))
conversion_df <- filter(mut_annot_reduced,Is.in.cancerGeneList.AsAliases==T)
write.table(conversion_df, "./DATA/PROCESSED_DATA/p_gene-names_conversion_table.tsv",
            sep="\t", quote=FALSE, row.names = FALSE)



### ... Binary matrix creation for oncogenic alterations ----
# Maf annotated with oncoKB.
mut_annot_data <- read.delim("./DATA/PROCESSED_DATA/p_OKB-annotated_MAF_data_mutations.oncokb.txt",
                             sep="\t")
# Conversion table created before in this script.
conversion_df <- read.delim("./DATA/PROCESSED_DATA/p_gene-names_conversion_table.tsv",
                            sep="\t")
# Raw clinical data table, downloaded from cBioPortal.
clinical_data <- read.table("./DATA/RAW_DATA/msk_met_2021/data_clinical_sample.txt",
                            sep="\t",
                            header=T)[c("SAMPLE_ID","ORGAN_SYSTEM", "SUBTYPE_ABBREVIATION", "SAMPLE_TYPE", "METASTATIC_SITE")]
clinical_data$ORGAN_SYSTEM <- str_replace_all(clinical_data$ORGAN_SYSTEM," ","-")
maf <- mut_annot_data
### ... Changing outdated names to new ones
maf$Hugo_Symbol <-  unlist(apply(maf,1,function(row){
  if(any(grepl(paste0("\\b",row["Hugo_Symbol"],"\\b"),conversion_df$Hugo_Symbol))){
    newname <- conversion_df[grep(paste0("\\b",row["Hugo_Symbol"],"\\b"),conversion_df$Hugo_Symbol),"New.Hugo.Symbol"]
    return(newname)
  }else{
    return(row["Hugo_Symbol"])
  }
}))
oncogenic_maf <- filter(maf,ONCOGENIC=="Likely Oncogenic" | ONCOGENIC=="Oncogenic")
write.table(oncogenic_maf, "./DATA/PROCESSED_DATA/p_oncogenic-MAF.tsv",
            sep="\t", quote=FALSE, row.names = FALSE)
# Raw Binary Matrix Generation ----
oncogenic_maf <- merge(oncogenic_maf,clinical_data,by.x="Tumor_Sample_Barcode","SAMPLE_ID")
matnames <- list.files(path="./DATA/PROCESSED_DATA/processed_cna_matrices/",
                       pattern="cnv_data*")
cancernames <- str_remove_all(str_remove_all(matnames,"cnv_data_cna_hg19_abs0.2_"),".txt")
cancernames <- str_replace_all(cancernames," ","-")
cnv_mats <- lapply(matnames,function(x) {
  cm_mat <- read.csv(sprintf("./DATA/PROCESSED_DATA/processed_cna_matrices/%s",x),
                     sep = "\t",
                     header = TRUE)
  rownames(cm_mat) <- cm_mat$Sample
  cm_mat$Sample <- NULL
  return(cm_mat)
})
names(cnv_mats) <- cancernames
# mutcnv_mats.genes <- lapply(mutcnv_mats,function(genenames){return(genenames[[2]])})
# mutcnv_mats <- lapply(mutcnv_mats,function(genenames){return(genenames[[1]])})

data_list <- split(oncogenic_maf,oncogenic_maf$ORGAN_SYSTEM)
names(data_list) <- str_replace_all(names(data_list)," ","-")
raw_binary_matrixes_list <- lapply(names(data_list),function(tiss){
  tissue_df <- data_list[[tiss]]
  clinical_samples <- filter(clinical_data,ORGAN_SYSTEM==tiss)$SAMPLE_ID
  mutated_genes_in_samples <- tissue_df[,c("Tumor_Sample_Barcode","Hugo_Symbol")]
  names(mutated_genes_in_samples) <- c("SAMPLE_ID","Gene")
  ### ... Storing MAF file gene names in vector (genes with oncogenic mutations)
  mutGenes <- unique(sort(mutated_genes_in_samples$Gene))
  mutated_genes_in_samples$Gene <- paste0(mutated_genes_in_samples$Gene,"_mutation")
  ### ... Storing SP mutcnv matrix for the correspondent tissue
  cnv_df <- cnv_mats[[tiss]]
  ### ... Elongating matrix (from wide to long)
  print(class(cnv_df))
  cnv_df_long <- cnv_df %>%
    tibble::rownames_to_column(var = "SAMPLE_ID") %>%
    pivot_longer(-SAMPLE_ID,
                 names_to="Gene",
                 values_to="count")
  cnv_df_long[c("Gene","Alter")] <- str_split_fixed(cnv_df_long$Gene,"_",2)
  ### ... Checking that we don't have both loss & gain cnas in the same patient and gene
  checking_cnas <- cnv_df_long %>%
    group_by(SAMPLE_ID,Gene) %>%
    summarise(L=sum(count))
  any(checking_cnas$L==2)
  ### ... Keeping only cna information (CNAS TAHT EXIST) to merge with mut information later
  cnv_df_long <- data.frame(filter(cnv_df_long,count>0))
  ### ... Checking if we need to update Hugo-Symbols (FALS)
  any(cnv_df_long$Gene %in% conversion_df$Hugo_Symbol)
  cnv_df_long$count <- NULL
  ### ... Storing genes that have at least 1 cna alteration in vector
  cnaGenes <- unique(sort(cnv_df_long$Gene))
  cnv_df_long$Gene <- paste0(cnv_df_long$Gene,"_",cnv_df_long$Alter)
  cnv_df_long$Alter <- NULL
  ### ... Merging oncogenic mutations with cna information (ALL ONCOGENIC FROM MAF)
  tissue_df_alterations <- bind_rows(mutated_genes_in_samples,cnv_df_long)
  ### ... List of all altered genes (at least 1 mut or 1 gain or 1 loss)
  altered_genes <- unique(sort(str_split_fixed(tissue_df_alterations$Gene,"_",2)[,1]))
  print(tiss)
  print("Total alterations")
  print(nrow(tissue_df_alterations))
  print("Total mutations")
  print(sum(grepl("_mutation",tissue_df_alterations$Gene)))
  print("Total Gain")
  print(sum(grepl("_Gain",tissue_df_alterations$Gene)))
  print("Total Loss")
  print(sum(grepl("_Loss",tissue_df_alterations$Gene)))
  ### ... Creating binary matrix (long to wide)
  binary_matrix_alterations <- pivot_wider(tissue_df_alterations,
                                           names_from = "Gene",
                                           values_from = "Gene",
                                           values_fill = 0,
                                           values_fn = function(x) 1)
  binary_matrix_alterations <- as.data.frame(binary_matrix_alterations)
  rownames(binary_matrix_alterations) <- binary_matrix_alterations$SAMPLE_ID
  binary_matrix_alterations$SAMPLE_ID <- NULL
  ### ... List of column names of all columns that should be in the final binary matrix (each gene x3 (mut, gain, loss))
  final_alterations <- c(paste0(altered_genes,"_mutation"),paste0(altered_genes,"_Gain"),paste0(altered_genes,"_Loss"))
  missing_alterations <- setdiff(final_alterations,colnames(binary_matrix_alterations))
  binary_matrix_alterations[missing_alterations] <- 0
  binary_matrix_alterations <- binary_matrix_alterations[,order(colnames(binary_matrix_alterations))]
  ### ... Adding missing samples
  missing_samples <- setdiff(clinical_samples,rownames(binary_matrix_alterations))
  missing_tissue_df <- setNames(data.frame(matrix(0,ncol=ncol(binary_matrix_alterations),nrow=length(missing_samples)),row.names=missing_samples),
                                names(binary_matrix_alterations))
  binary_matrix_alterations <- rbind(binary_matrix_alterations,missing_tissue_df)
  return(binary_matrix_alterations)
})
names(raw_binary_matrixes_list) <- names(data_list)


saveRDS(raw_binary_matrixes_list,
        "./DATA/PROCESSED_DATA/000_raw_binary_matrixes.RDS")



### ... Binary matrix filtering by genes in OncoKB ----
# List with genes in oncoKB, created before.
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
raw_binary_matrixes_list <- readRDS("./DATA/PROCESSED_DATA/000_raw_binary_matrixes.RDS")
oncokb_binary_matrixes_list  <- lapply(raw_binary_matrixes_list,function(mat){
  okb_mat <- mat[-which(str_split_fixed(names(mat),"_",2)[,1]%!in%cancgenedf$Gene)]
  return(okb_mat)})
saveRDS(oncokb_binary_matrixes_list,
        "./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")



### .... Clinical data table cleaning, modifying, ----
## Adding information of use for all the models 
## STAGE_PWOWM, TREATMENT, METASTATIC AGGRESSIVENESS (MET_COUNT), 
## METASTATIC AGGRESSIVENESS (TIMING), from cBioPortal.
clinical_data <- read.table("./DATA/RAW_DATA/msk_met_2021/data_clinical_sample.txt",
                            sep="\t",
                            header=T)
# For treatment information, Supplementary Table 1B of MSK-MET publication.
suppl_table_1_2 <- read.xlsx("./DATA/RAW_DATA/msk_met_2021/supplementary_data_paper/1-s2.0-S0092867422000034-mmc1.xlsx",
                             sheet = 2, startRow = 2, colNames = T)
# for timing of metastases information, from cBioPortal.
diagnosis_data <- read.table("./DATA/RAW_DATA/msk_met_2021/data_timeline_diagnosis.txt",
                             sep="\t", header=T)
# For age at diagnosis, from cBioPortal.
patient_data <- read.table("./DATA/RAW_DATA/msk_met_2021/data_clinical_patient.txt",
                           sep="\t", header=T)
# For metastatic biopsy location, Supplementary Table 1C of MSK-MET publication.
suppl_table_1_3 <- read.xlsx("./DATA/RAW_DATA/msk_met_2021/supplementary_data_paper/1-s2.0-S0092867422000034-mmc1.xlsx",
                             sheet = 3, startRow = 2, colNames = T)
## For number of alterations per sample
oncokb_binary_matrixes_list <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
#  Changing spaces to - 
clinical_data$SUBTYPE_ABBREVIATION <- str_replace_all(clinical_data$SUBTYPE_ABBREVIATION," ","-")
clinical_data$ORGAN_SYSTEM <- str_replace_all(clinical_data$ORGAN_SYSTEM," ","-")
# STAGE2 column with 3-level stage classification (Pri_W/WO_Meta)
stage2 <- clinical_data[c("SAMPLE_ID","SAMPLE_TYPE","MET_COUNT")]
stage2$STAGE2 <- ifelse(stage2$SAMPLE_TYPE=="Primary" & stage2$MET_COUNT==0,"Primary_WO_Metastasis",
                        ifelse(stage2$SAMPLE_TYPE=="Primary" & stage2$MET_COUNT>0,"Primary_W_Metastasis","Metastasis"))
stage2 <- setNames(stage2[c(1,4)],c("SAMPLE_ID","STAGE_PWOWM"))
# Retrieving binary TREATMENT information
treatment_table <- setNames(suppl_table_1_2[,c("sample_id","prior_treatement")],c("EVENT","TREATMENT"))
treatment_table$TREATMENT <- ifelse(treatment_table$TREATMENT==1,"Treatment","Non_Treatment")
names(treatment_table) <- c("SAMPLE_ID","TREATMENT")
# Age at sequencing
patient_data <- merge(patient_data[c("PATIENT_ID","AGE_AT_SEQUENCING")],clinical_data[c("PATIENT_ID","SAMPLE_ID")])
# Modifying diagnosis date table for TIMING OF METASTASES INFO
# As in MET_COUNT model,  Removing NOT MAPPED LYMPH patients
long_diagnosis_data <- merge(merge(diagnosis_data,clinical_data[c("SAMPLE_ID","PATIENT_ID")],by="PATIENT_ID"),stage2,by="SAMPLE_ID")
long_diagnosis_data <- filter(long_diagnosis_data,SUBTYPE!="Lymph"&SUBTYPE!="Regional Lymph")
# Eliminating same-date diagnosis (leaving only one per date)
clean_diagnosis_data <- unique(long_diagnosis_data[c("SAMPLE_ID","START_DATE")])
# Checking that there are no PWOWM patients in the diagnoses data table
any(filter(clean_diagnosis_data,SAMPLE_ID%in%filter(stage2,STAGE_PWOWM=="Primary_WO_Metastasis")[c("SAMPLE_ID")]))
# Couting diagnosis of metastasis in order (1 - first (earliest), 2 - second ...)
clean_diagnosis_data <- clean_diagnosis_data %>% group_by(SAMPLE_ID) %>% mutate(DIAGNOSIS_ORDER=rank(START_DATE,ties.method = "first"))
# Transforming start date FOR EACH PATIENT to 0, ant the rest, proportionally
clean_diagnosis_data <- clean_diagnosis_data %>% group_by(SAMPLE_ID) %>% mutate(Zero_Start = START_DATE-min(START_DATE))
### ... Transforming dates (days) into years
clean_diagnosis_data$YEARS_BTW_METAS <- round(clean_diagnosis_data$Zero_Start/365,4)
clean_diagnosis_data <- merge(clean_diagnosis_data,melt(table(clean_diagnosis_data$SAMPLE_ID),varnames="SAMPLE_ID",value.name="DIAGNOSES"))
# Keeping only patients with Order = 2, (2nd metastasis)
# Their Date_Years will indicate the time between 1st and 2nd metastases
timing_diagnosis_data <- filter(clean_diagnosis_data,DIAGNOSIS_ORDER==2)



## ... Retrieving metastatic biopsy location
metastatic_location_data <- clinical_data[c("SAMPLE_ID","SAMPLE_TYPE","ORGAN_SYSTEM","ONCOTREE_CODE","PRIMARY_SITE","METASTATIC_SITE")]
metastatic_location_list <- split(filter(metastatic_location_data,SAMPLE_TYPE=="Metastasis"),
                                  filter(metastatic_location_data,SAMPLE_TYPE=="Metastasis")$METASTATIC_SITE)
ehr_pathology_to_organs <- suppl_table_1_3[c(1,3)]
ehr_pathology_to_organs$`Free-text.Description` <- str_remove_all(ehr_pathology_to_organs$`Free-text.Description`,c("\\["))
ehr_pathology_to_organs$`Free-text.Description` <- str_remove_all(ehr_pathology_to_organs$`Free-text.Description`,c("\\]"))
ehr_pathology_to_organs$`Free-text.Description` <- str_remove_all(ehr_pathology_to_organs$`Free-text.Description`,c("'"))
pathology_to_organs_list <- split(ehr_pathology_to_organs$`Free-text.Description`,ehr_pathology_to_organs$`MSK-Met.Organ.Sites`)
pathology_to_organs_list <- lapply(pathology_to_organs_list,function(organ_list){return(data.frame(Description=str_split(organ_list,", ")[[1]]))})
pathology_to_organs_df <- map_df(pathology_to_organs_list,~as.data.frame(.x),.id="Organ_Site")
pathology_to_organs_df <- unique(pathology_to_organs_df)
metastatic_location_df <- filter(metastatic_location_data,SAMPLE_TYPE=="Metastasis")
metastatic_location_df$PRIM_CAP <- toupper(metastatic_location_df$PRIMARY_SITE)
metastatic_location_df <- merge(metastatic_location_df,
                                pathology_to_organs_df,
                                by.x = 'PRIM_CAP',
                                by.y = 'Description') %>% 
  rename('Prim.Organ' = 'Organ_Site')
metastatic_location_df$META_CAP <- toupper(metastatic_location_df$METASTATIC_SITE)
metastatic_location_df <- merge(metastatic_location_df,
                                pathology_to_organs_df,
                                by.x = 'META_CAP', by.y = 'Description', all.x = T) %>% 
  rename('Meta.Organ' = 'Organ_Site')
metastatic_location_df$Meta.Organ <- ifelse(
  !is.na(metastatic_location_df$Meta.Organ), metastatic_location_df$META_CAP,
  ifelse(metastatic_location_df$META_CAP == 'BLADDER/UT', 'BLADDER_OR_URINARY_TRACT',
         ifelse(metastatic_location_df$META_CAP == 'CNS/BRAIN', 'CNS_BRAIN',
                ifelse(metastatic_location_df$META_CAP == 'DISTANT LN', 'LYMPH',
                       ifelse(metastatic_location_df$META_CAP == 'FEMALE GENITAL', 'GENITAL_FEMALE',
                              ifelse(metastatic_location_df$META_CAP == 'INTRA-ABDOMINAL', 'OTHER',
                                     ifelse(metastatic_location_df$META_CAP == 'MALE GENITAL', 'GENITAL_MALE',
                                            ifelse(metastatic_location_df$META_CAP == 'UNSPECIFIED', 'UNSPECIFIED', 'MISTAKE_CHECK'))))))))
metastatic_location_df$BIOPSY_LOCATION <- ifelse(grepl('Unknown', metastatic_location_df$PRIMARY_SITE), 'Unknown primary',
                                                 ifelse(metastatic_location_df$Meta.Organ == 'LYMPH', 'Lymph',
                                                        ifelse(metastatic_location_df$Meta.Organ == 'OTHER', 'Other',
                                                               ifelse(metastatic_location_df$Meta.Organ == 'UNSPECIFIED', 'Unspecified',
                                                                      ifelse(metastatic_location_df$Meta.Organ == metastatic_location_df$Prim.Organ, 'Local',
                                                                             'Distal')))))
table(metastatic_location_df$BIOPSY_LOCATION)
samples_per_biopsy_location <- metastatic_location_df %>% 
  group_by(PRIMARY_SITE, METASTATIC_SITE, BIOPSY_LOCATION) %>%
  summarise(N_SAMPLES = n())
# Some samples have been duplicated
meta_location_unique <- unique(metastatic_location_df[,c("SAMPLE_ID", "PRIMARY_SITE", "METASTATIC_SITE", "BIOPSY_LOCATION")])
# Some samples are still duplicated
dupped_samples <- as.character(filter(as.data.frame(table(meta_location_unique$SAMPLE_ID)), Freq == 2)$Var1)
to_add <- filter(meta_location_unique, SAMPLE_ID %in% dupped_samples)
to_add <- filter(to_add, BIOPSY_LOCATION == 'Local') # Manually checked that duplicated samples are local biopsies
meta_location_unique <- meta_location_unique[-which(meta_location_unique$SAMPLE_ID %in% dupped_samples),] # Removing duplicated samples
# Adding the non-duplicated rows
meta_location_unique <- bind_rows(meta_location_unique,
                                  to_add)
samples_per_biopsy_location <- meta_location_unique %>% 
  group_by(PRIMARY_SITE, METASTATIC_SITE, BIOPSY_LOCATION) %>% 
  summarise('N_SAMPLES' = n())
write.table(samples_per_biopsy_location, "./DATA/PROCESSED_DATA/p_metastatic-biopsy-locations.tsv",
            sep = "\t", row.names = F, quote = F)



### ... Adding info about samples having oncogenic alterations
alterations_per_sample <- lapply(unname(oncokb_binary_matrixes_list),function(mat){return(rowSums(mat))})
alterations_per_sample <- data.frame(N_ONCOGENIC_ALTERATIONS=unlist(alterations_per_sample,recursive=F))


clean_clinical_data <- clinical_data[c("SAMPLE_ID","ORGAN_SYSTEM","SUBTYPE_ABBREVIATION","SAMPLE_TYPE","MET_COUNT",
                                       "PRIMARY_SITE","METASTATIC_SITE","IS_DIST_MET_MAPPED","TUMOR_PURITY","SAMPLE_COVERAGE")]
names(clean_clinical_data) <- c("SAMPLE_ID","CANC_TYPE","CANC_SUBTYPE","STAGE_PM","MET_COUNT","PRIM_SITE","MET_SITE",
                                "IS_DIST_MET_MAPPED","TUMOR_PURITY","SAMPLE_COVERAGE")
clean_clinical_data <- merge(clean_clinical_data,stage2,by="SAMPLE_ID")
clean_clinical_data <- merge(clean_clinical_data,treatment_table,by="SAMPLE_ID")
clean_clinical_data <- merge(clean_clinical_data,timing_diagnosis_data[c("SAMPLE_ID","YEARS_BTW_METAS","DIAGNOSES")],all.x=T)
clean_clinical_data <- merge(clean_clinical_data,patient_data[c("SAMPLE_ID","AGE_AT_SEQUENCING")])



### ... Creating maf with CLONALITY INFORMATION ----
maf <- read.delim("./DATA/RAW_DATA/MSK_MET_msk_impact_facets_annotated.ccf.maf",
                  sep="\t", header=T)
conversion_df <- read.delim("./DATA/PROCESSED_DATA/p_gene-names_conversion_table.tsv",
                            sep="\t")
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
### ... Changing outdated names to new ones
maf$Hugo_Symbol <-  unlist(apply(maf,1,function(row){
  if(any(grepl(paste0("\\b",row["Hugo_Symbol"],"\\b"),conversion_df$Hugo_Symbol))){
    newname <- conversion_df[grep(paste0("\\b",row["Hugo_Symbol"],"\\b"),conversion_df$Hugo_Symbol),"New.Hugo.Symbol"]
    return(newname)
  }else{
    return(row["Hugo_Symbol"])
  }
}))
clonal_proportion <- melt(table(maf$clonality),varnames="C",value.name="Variants")
clonal_proportion <- bind_rows(clonal_proportion,data.frame(C="NA",Variants=sum(is.na(maf$clonality))))
clonal_proportion$'%' <- round(prop.table(clonal_proportion$Variants)*100,2)
oncogenic_maf <- filter(maf,oncogenic=="Likely Oncogenic" | oncogenic=="Oncogenic")
clonal_maf <- filter(oncogenic_maf,clonality=="CLONAL"|clonality=="SUBCLONAL")
clonal_functional_maf <- filter(clonal_maf,Hugo_Symbol%in%cancgenedf$Gene)
clonal_functional_clinical_maf <- merge(clonal_functional_maf,clean_clinical_data,by.x="Tumor_Sample_Barcode",by.y="SAMPLE_ID",all.x=T)
clonal_functional_clinical_maf <- merge(clonal_functional_clinical_maf,cancgenedf,by.x="Hugo_Symbol",by.y="Gene",all.x=T)
write.table(clonal_functional_clinical_maf, './DATA/PROCESSED_DATA/p-clonal-functional-clinical-maf.tsv',
            sep = '\t', row.names = F, quote = F)


maf_for_models <- unique(clonal_functional_clinical_maf[c("Tumor_Sample_Barcode","STAGE_PM","Function",
                                                          "Hugo_Symbol","clonality","CANC_TYPE","TREATMENT")])
write.table(maf_for_models,
            "./DATA/PROCESSED_DATA/p_clonal-MAF_filt-for-model.tsv",
            sep="\t", quote=F, row.names=F)



### ... Table with proportion of clonal/subclonal variants across stages
clonal_fraction_bytissue_allvariants <- melt(table(clonal_functional_clinical_maf$Hugo_Symbol,clonal_functional_clinical_maf$CANC_TYPE,
                                                   clonal_functional_clinical_maf$clonality,clonal_functional_clinical_maf$STAGE_PM,
                                                   clonal_functional_clinical_maf$TREATMENT,clonal_functional_clinical_maf$Function),
                                             varnames=c("Gene","Tissue","Clonality","Stage","Treatment","Function"),
                                             value.name="All.Variants")
clonal_fraction_bytissue_modelvariants <- melt(table(maf_for_models$Hugo_Symbol,maf_for_models$CANC_TYPE,
                                                     maf_for_models$clonality,maf_for_models$STAGE_PM,maf_for_models$TREATMENT,
                                                     maf_for_models$Function),
                                               varnames=c("Gene","Tissue","Clonality","Stage","Treatment","Function"),
                                               value.name="Model.Variants")
clonal_fraction_bytissue <- merge(clonal_fraction_bytissue_allvariants,clonal_fraction_bytissue_modelvariants)
rm(clonal_fraction_bytissue_allvariants,clonal_fraction_bytissue_modelvariants)
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Gene,Tissue,Stage,Treatment,Function) %>% 
  mutate(Gene_CF_All.Variants_wStagewTreatwFunction = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         Gene_CF_Model.Variants_wStagewTreatwFunction = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Gene,Tissue,Stage,Function) %>% 
  mutate(Gene_CF_All.Variants_wStagewFunction = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         Gene_CF_Model.Variants_wStagewFunction = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Gene,Tissue,Stage,Treatment) %>% 
  mutate(Gene_CF_All.Variants_wStagewTreat = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         Gene_CF_Model.Variants_wStagewTreat = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Gene,Tissue,Stage) %>% 
  mutate(Gene_CF_All.Variants_wStage = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         Gene_CF_Model.Variants_wStage = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))

clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Tissue,Stage,Treatment,Function) %>% 
  mutate(CF_All.Variants_wStagewTreatwFunction = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         CF_Model.Variants_wStagewTreatwFunction = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Tissue,Stage,Function) %>% 
  mutate(CF_All.Variants_wStagewFunction = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         CF_Model.Variants_wStagewFunction = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Tissue,Stage,Treatment) %>% 
  mutate(CF_All.Variants_wStagewTreat = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         CF_Model.Variants_wStagewTreat = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
clonal_fraction_bytissue <- clonal_fraction_bytissue %>% 
  group_by(Tissue,Stage) %>% 
  mutate(CF_All.Variants_wStage = sum(All.Variants[Clonality=="CLONAL"])/sum(All.Variants),
         CF_Model.Variants_wStage = sum(Model.Variants[Clonality=="CLONAL"])/sum(Model.Variants))
write.table(clonal_fraction_bytissue, "./DATA/PROCESSED_DATA/p_clonality_by-tissue.tsv",
            sep="\t", quote=F, row.names=F)



################################################################################
############ ADDING MORE THINGS TO CLINICAL DATA TABLE #########################
clean_clinical_data <- merge(clean_clinical_data,
                             alterations_per_sample,
                             by.x="SAMPLE_ID",
                             by.y="row.names")
clean_clinical_data <- merge(clean_clinical_data,
                             no_altered_patients)
################################################################################


## Checking clonality & metastatic biopsy location things ----
metastatic_biopsy_locations <- read.delim("./DATA/PROCESSED_DATA/p_metastatic-biopsy-locations.tsv",
                                          sep="\t", header=T)
matching_primary_sites <- data.frame(Clinical_data_prims=sort(unique(filter(clean_clinical_data,STAGE_PM=="Metastasis")$PRIM_SITE)),
                                     Meta_sites_prims=sort(unique(metastatic_biopsy_locations$PRIMARY_SITE)))
matching_metastatic_sites <- data.frame(Clinical_data_metas=sort(unique(filter(clean_clinical_data,STAGE_PM=="Metastasis")$MET_SITE)),
                                        Meta_sites_metas=sort(unique(metastatic_biopsy_locations$METASTATIC_SITE)))
metastatic_biopsy_locations <- merge(metastatic_biopsy_locations,
                                     matching_primary_sites,
                                     by.x="PRIMARY_SITE",
                                     by.y="Meta_sites_prims")
metastatic_biopsy_locations <- merge(metastatic_biopsy_locations,
                                     matching_metastatic_sites,
                                     by.x="METASTATIC_SITE",
                                     by.y="Meta_sites_metas")

location_data <- merge(clean_clinical_data[c("SAMPLE_ID","PRIM_SITE","MET_SITE","STAGE_PM")],
                       metastatic_biopsy_locations[c("Clinical_data_prims","Clinical_data_metas","BIOPSY_LOCATION")],
                       by.x=c("PRIM_SITE","MET_SITE"),
                       by.y=c("Clinical_data_prims","Clinical_data_metas"),all.x=T)

clonal_functional_clinical_maf <- merge(clonal_functional_clinical_maf,
                                        location_data[c("SAMPLE_ID","BIOPSY_LOCATION")],
                                        by.x=c("Tumor_Sample_Barcode"),
                                        by.y=c("SAMPLE_ID"),
                                        all.x=T)
write.table(clonal_functional_clinical_maf, "./DATA/PROCESSED_DATA/p_clonal-MAF.tsv",
            sep="\t", quote=F, row.names=F)


################################################################################
############ ADDING MORE THINGS TO CLINICAL DATA TABLE #########################
clean_clinical_data <- merge(clean_clinical_data,location_data[c("SAMPLE_ID","BIOPSY_LOCATION")])
################################################################################


################################################################################
############ FINALLY SAVING CLINICAL DATA TABLE #########################
write.table(clean_clinical_data,"./DATA/PROCESSED_DATA/p_clinical-data.tsv",
            sep = "\t", row.names = F, quote = F)


### ... Loading libraries ----
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(tibble)
library(purrr)
"%!in%" <- Negate("%in%")



### ... Input files ----
clinical_data <- read.table('./DATA/PROCESSED_DATA/p_clinical-data.tsv',
                            sep = "\t", header=T)
cancgenedf <- read.csv('./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv', 
                       sep = "\t",header = TRUE, stringsAsFactors = FALSE)
conversion_df <- read.delim('./DATA/PROCESSED_DATA/p_gene-names_conversion_table.tsv',
                            sep = '\t')



### ... Processing MAF and CNV files ----
oncogenic_maf <- read.delim("./DATA/PROCESSED_DATA/p_oncogenic-MAF.tsv", sep="\t")
oncogenic_maf <- merge(oncogenic_maf,
                       clinical_data,
                       by.x="Tumor_Sample_Barcode",
                       by.y = "SAMPLE_ID")
matnames <- list.files(path="./DATA/PROCESSED_DATA/processed_cna_matrices/",
                       pattern="cnv_data*")
cancernames <- str_remove_all(str_remove_all(matnames,"cnv_data_cna_hg19_abs0.2_"),".txt")
cancernames <- str_replace_all(cancernames," ","-")
# CNV matrices
cnv_mats <- lapply(matnames,function(x) {
  cm_mat <- read.csv(sprintf("./DATA/PROCESSED_DATA/processed_cna_matrices/%s",x),
                     sep = "\t",
                     header = TRUE)
  rownames(cm_mat) <- cm_mat$Sample
  cm_mat$Sample <- NULL
  return(cm_mat)
})
names(cnv_mats) <- cancernames
# Mutations
oncogenic_maf$HGVSp_Short <- str_replace_all(oncogenic_maf$HGVSp_Short, "_", "-")
data_list <- split(oncogenic_maf,oncogenic_maf$CANC_TYPE)
names(data_list) <- str_replace_all(names(data_list)," ","-")



### ... Creating raw binary matrices ----
raw_binary_matrixes_list <- lapply(names(data_list),function(tiss){
  tissue_df <- data_list[[tiss]] # Select one cancer type
  clinical_samples <- filter(clinical_data,CANC_TYPE==tiss)$SAMPLE_ID
  mutated_genes_in_samples <- tissue_df[,c("Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short")]
  names(mutated_genes_in_samples) <- c("SAMPLE_ID","Gene","HGVSp_Short")
  ### ... Storing MAF file gene names in vector (genes with oncogenic mutations)
  mutGenes <- unique(sort(mutated_genes_in_samples$Gene))
  mutated_genes_in_samples$Gene <- paste0(mutated_genes_in_samples$Gene,"_",mutated_genes_in_samples$HGVSp_Short,"_mutation")
  mutated_genes_in_samples$HGVSp_Short <- NULL
  ### ... Storing mutcnv matrix for the correspondent tissue
  cnv_df <- cnv_mats[[tiss]]
  ### ... Elongating matrix (from wide to long)
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
  ### ... Keeping only cna information (CNAS THAT EXIST) to merge with mut information later
  cnv_df_long <- data.frame(filter(cnv_df_long,count>0))
  ### ... Checking if we need to update Hugo-Symbols (FALSE)
  any(cnv_df_long$Gene %in% conversion_df$Hugo_Symbol)
  cnv_df_long$count <- NULL
  
  ### ... Storing genes that have at least 1 cna alteration in vector
  cnaGenes <- unique(sort(cnv_df_long$Gene))
  print(paste0("For ",tiss," tissue, ",length(mutGenes)," genes have mutations, ",length(cnaGenes)," have CNAs. ",length(intersect(mutGenes,cnaGenes))," have both."))
  cnv_df_long$Gene <- paste0(cnv_df_long$Gene,"_",cnv_df_long$Alter)
  cnv_df_long$Alter <- NULL
  
  ### ... Merging oncogenic mutations with cna information (ALL ONCOGENIC FROM MAF)
  tissue_df_alterations <- bind_rows(mutated_genes_in_samples,cnv_df_long)
  ### ... List of all altered genes (at least 1 mut or 1 gain or 1 loss)
  altered_genes <- unique(sort(str_split_fixed(tissue_df_alterations$Gene,"_",2)[,1]))
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
  final_alterations <- c(paste0(altered_genes,"_Gain"),paste0(altered_genes,"_Loss"))
  missing_alterations <- setdiff(final_alterations,colnames(binary_matrix_alterations))
  binary_matrix_alterations[missing_alterations] <- 0
  binary_matrix_alterations <- binary_matrix_alterations[,order(colnames(binary_matrix_alterations))]
  ### ... Adding missing samples
  missing_samples <- setdiff(clinical_samples,rownames(binary_matrix_alterations))
  missing_tissue_df <- setNames(data.frame(matrix(0,ncol=ncol(binary_matrix_alterations),
                                                  nrow=length(missing_samples)),
                                           row.names=missing_samples),
                                names(binary_matrix_alterations))
  binary_matrix_alterations <- rbind(binary_matrix_alterations,missing_tissue_df)
  return(binary_matrix_alterations)
})
names(raw_binary_matrixes_list) <- names(data_list)



### ... Binary matrix filtering by genes in OncoKB ----
# raw_binary_matrixes_list <- readRDS("./inputs/Pos-Specific_RAW-binary-matrices.RDS")
oncokb_binary_matrixes_list  <- lapply(raw_binary_matrixes_list,function(mat){
  okb_mat <- mat[-which(str_split_fixed(names(mat),"_",2)[,1]%!in%cancgenedf$Gene)]
  return(okb_mat)})



### ... Removing in each matrix genes that have CNA (loss/gain) but no mutations ----
oncokb_bm_clean <- lapply(oncokb_binary_matrixes_list, function(mat){
  mut_genes <- colnames(mat)[grepl('_mutation', colnames(mat))]
  mut_genes <- unique(str_split_fixed(mut_genes, '_', 3)[,1])
  pattern <- paste(mut_genes, collapse = '|')
  return(mat[,colnames(mat)[grepl(pattern, colnames(mat))]])
})
saveRDS(oncokb_bm_clean, './DATA/PROCESSED_DATA/Pos-Specific_oncokb_binary-matrices.RDS')

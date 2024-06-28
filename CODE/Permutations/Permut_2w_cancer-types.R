
### ... PERMUTATION SCRIPT TO ASSESS SIGNIFICANCE OF INTERACTIONS ----
library(vegan)
library(stringr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
options(dplyr.summarise.inform = FALSE)



## Functions/Variables ----
### ... Changeable
n_permutations <- 100
# SPLITMOD <- "Tissue-Stage-PM"
SPLITMOD <- args[1]
# if(SPLITMOD==NA){
#   break
# }
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
## Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



##### 1. CREATING PERMUTATION MATRIXES #####
### ... Uploading files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)



### ... Subsetting clinical data according to model ----
clinical_data <- filter(clinical_data,N_ONCOGENIC_ALTERATIONS>0)
if (SPLITMOD == "Tissue"){
  model_split_data <- "CANC_TYPE"
}else if (SPLITMOD == "Subtype"){
  model_split_data <- c("CANC_TYPE","CANC_SUBTYPE")
}else if (SPLITMOD == "Tissue-Stage-PM"){
  model_split_data <- c("CANC_TYPE","STAGE_PM")
}else if (SPLITMOD == "Subtype-Stage-PM"){
  model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PM")
}else if (SPLITMOD == "Tissue-Stage-PWOWM"){
  model_split_data <- c("CANC_TYPE","STAGE_PWOWM")
}else if (SPLITMOD == "Subtype-Stage-PWOWM"){
  model_split_data <- c("CANC_TYPE","CANC_SUBTYPE","STAGE_PWOWM")
}
clinical_data <-
  clinical_data %>%
  group_by_at(model_split_data) %>%
  mutate(G = cur_group_id(),
         name = ifelse(SPLITMOD == "Subtype-Stage-PWOWM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PWOWM,sep="."),
                       ifelse(SPLITMOD == "Subtype-Stage-PM",paste(CANC_TYPE,CANC_SUBTYPE,STAGE_PM,sep="."),
                              ifelse(SPLITMOD == "Tissue-Stage-PWOWM",paste(CANC_TYPE,NA,STAGE_PWOWM,sep="."),
                                     ifelse(SPLITMOD == "Tissue-Stage-PM",paste(CANC_TYPE,NA,STAGE_PM,sep="."),
                                            ifelse(SPLITMOD == "Subtype",paste(CANC_TYPE,CANC_SUBTYPE,NA,sep="."),
                                                   paste(CANC_TYPE,NA,NA,sep=".")))))))



### ... Subsetting matrixes according to model ----
ready_bm <- unlist(
  mapply(function(binary_matrix,clinical_table){
    bm <- subset(binary_matrix,
                 rownames(binary_matrix)%in%clinical_table$SAMPLE_ID)
    bm <- merge(bm,
                clinical_table[c("SAMPLE_ID","name")],
                by.x="row.names",
                by.y="SAMPLE_ID")
    rownames(bm) <- bm$Row.names
    bm$Row.names <- NULL
    split_bm <- split(bm[1:(ncol(bm)-1)],bm$name)
    return(split_bm)
    },
    unname(binary_mats),
    split(clinical_data,clinical_data$CANC_TYPE),
    SIMPLIFY=F),
  recursive=F)



### ... Running permutation function ----
permuted_matrixes <- lapply(seq_along(ready_bm), function(binary_mat_index){
  binary_mat <- ready_bm[[binary_mat_index]]
  binary_mat_name <- names(ready_bm)[binary_mat_index]
  print(binary_mat_name)
  binary_muts <- binary_mat[grep("_mutation",names(binary_mat))]
  binary_gain <- binary_mat[grep("_Gain",names(binary_mat))]
  binary_loss <- binary_mat[grep("_Loss",names(binary_mat))]
  permutation_matrixes_muts <- permatswap(binary_muts,times=n_permutations,method = "quasiswap")
  print("Mutation matrixes done")
  permutation_matrixes_gain <- permatswap(binary_gain,times=n_permutations,method = "quasiswap")
  print("Gain matrixes done")
  permutation_matrixes_loss <- permatswap(binary_loss,times=n_permutations,method = "quasiswap")
  print("Loss matrixes done")
  permutation_matrixes <- lapply(1:n_permutations,function(times){
    return(bind_cols(permutation_matrixes_muts$perm[[times]],
                     permutation_matrixes_gain$perm[[times]],
                     permutation_matrixes_loss$perm[[times]]))})
  return(permutation_matrixes)
})
names(permuted_matrixes) <- names(ready_bm)
saveRDS(permuted_matrixes,
        sprintf("./DATA/PROCESSED_DATA/PERMUTED_matrixes_%s.RDS",SPLITMOD))



#### 2. RUNNING 2WAY REGRESSION MODEL OVER MATRICES ####
### ... Loading files ----
permuted_matrixes <- readRDS(sprintf("./DATA/PROCESSED_DATA/PERMUTED_matrixes_%s.RDS",SPLITMOD))
### ... Function to create the input table from PERMUTED MATRIX ----
pseudo_oncogenic_model_input <- function(df,freqmut_threshold,freqcnv_threshold,pairs){
  if(is_tibble(df)){df <- as.data.frame(df)}
  size <- nrow(df)
  inputs_table <- mapply(function(gene,cna){
    ### ... Retrieving gene columns
    gene_column_names <- paste0(gene,c("_mutation","_Loss","_Gain"))
    gene_bm <- df[gene_column_names]
    ### ... Preparing to count occurrences
    gene_bm$cnvs <- paste0(gene_bm[,2],gene_bm[,3])
    counts <- paste0(gene_bm[,1],gene_bm$cnvs)
    ### ... If gain and loss found together ("x11"), added both as gain AND loss
    counts <- c(counts,rep("110",length(counts[counts=="111"])))
    counts[counts=="111"] <- "101"
    counts <- c(counts,rep("010",length(counts[counts=="011"])))
    counts[counts=="011"] <- "001"
    possib_occur <- c("010","110","000","100","001","101")
    ### ... Count of each occurrence
    N <- unlist(lapply(1:6,function(x) sum(counts==possib_occur[x])))
    ### ... Adding pseudocounts
    N <- N + 1
    N_df <- data.frame(t(N))
    names(N_df) <- c("NoMutLoss","MutLoss","NoMutWT","MutWT","NoMutGain","MutGain")
    ### ... Calculating mutation frequency from binary vector
    mutfreq <- sum(gene_bm[,1])/nrow(gene_bm)
    ### ... Calculating CNA loss frequency from binary vector
    lossfreq <- sum(gene_bm[,2])/nrow(gene_bm)
    ### ... Calculating CNA gain frequency from binary vector
    gainfreq <- sum(gene_bm[,3])/nrow(gene_bm)
    input_row <- bind_cols(data.frame(Gene = gene,
                                      Size = size,
                                      MutFreq = mutfreq,
                                      LossFreq = lossfreq,
                                      GainFreq = gainfreq,
                                      Mutation_filter = mutfreq >= freqmut_threshold,
                                      Loss_filter = lossfreq >= freqcnv_threshold,
                                      Gain_filter = gainfreq >= freqcnv_threshold,
                                      CNA_type = cna),
                           N_df)
    return(input_row)},
    pairs$Gene,
    pairs$CNA_type,
    SIMPLIFY=F)
  inputs_table <- bind_rows(inputs_table)
  return(inputs_table)
}
### ... Function to retrieve results ----
pseudo_oncogenic_retrieve_results <- function(df) {
  if(is.null(df)){return(NULL)}
  mut_df <- filter(df,Mutation_filter==TRUE)
  to_return <- list()
  for (cna in c("Loss","Gain")){
    cnv_df <- filter(mut_df,CNA_type==cna)
    cnv_df <- filter(cnv_df,get(paste0(cna,"_filter"))==T)
    if (nrow(cnv_df)==0){
      next
    }else{
      regmod_df <- as.data.frame(t(sapply(cnv_df$Gene, function(gene){
        datarow <- cnv_df[cnv_df$Gene==gene,]
        if (cna=="Loss"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutLoss,datarow$MutLoss)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=N)
        }else if (cna=="Gain"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutGain,datarow$MutGain)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=N)
        }
        result <- reg_model_2w(dataf)
        result <- c(result,cna)
        return(result)
      })))
      colnames(regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                               "NoMutWT","MutWT","NoMutCNV","MutCNV","CNA_type")
      regmod_df$Gene <- rownames(regmod_df)
      to_return[[cna]] <- regmod_df
    }
  }
  return(to_return)
}


### ... Running 2way regression model ----
real_results <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_analysis_%s_%s.tsv", FREQ, SPLITMOD),
                           sep="\t", header=T)
if(grepl("Tissue",SPLITMOD)){
  real_results$model <- paste0(real_results$Tissue,".NA.",real_results$Stage)
}else if(grepl("Subtype",SPLITMOD)){
  real_results$model <- paste0(real_results$Tissue,".",real_results$Subtype,".",real_results$Stage)
}
permuted_outputs <- lapply(seq_along(permuted_matrixes),function(permuted_matrix_index){
  permuted_matrixes_list <- permuted_matrixes[[permuted_matrix_index]]
  permuted_matrix_name <- names(permuted_matrixes)[permuted_matrix_index]
  print(permuted_matrix_name)
  real_tested_pairs <- filter(real_results,model==permuted_matrix_name)[c("Gene","CNA_type")]
  permuted_inputs <- lapply(permuted_matrixes_list,pseudo_oncogenic_model_input,freqmut_threshold,freqcnv_threshold,real_tested_pairs)
  permuted_outputs <- lapply(permuted_inputs,pseudo_oncogenic_retrieve_results)
  permuted_outputs <- mapply(function(results,inputs){
    cna_list <- lapply(results,function(cna_results){
      cna_results <- merge(cna_results,inputs[c("Gene","MutFreq","LossFreq","GainFreq","Size","CNA_type")],by=c("Gene","CNA_type"))
      return(cna_results)
    })
    cna_results <- map_df(cna_list,~as.data.frame(.x))
    return(cna_results)
  },permuted_outputs,permuted_inputs,SIMPLIFY = F)
  return(permuted_outputs)
})
names(permuted_outputs) <- paste0(names(permuted_matrixes),".")
permuted_outputs <- unlist(permuted_outputs,recursive=F)
permuted_outputs <- map_df(permuted_outputs,~as.data.frame(.x),.id="Tissue")
permuted_outputs[c("Tissue","Subtype","Stage","Permutation")] <- str_split_fixed(permuted_outputs$Tissue,"\\.",4)
permuted_outputs <- split(permuted_outputs,permuted_outputs$Permutation)
saveRDS(permuted_outputs,
        sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_outputs_%s_%s.RDS",FREQ,SPLITMOD))



### 3. COUNTING REAL SIGNIFICANTS #####
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
subtypes <- setNames(unique(clinical_data[c("CANC_TYPE","CANC_SUBTYPE")]),
                     c("Tissue","Subtype"))
real_results <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_analysis_%s_%s.tsv", FREQ, SPLITMOD),
                           sep="\t", header=T)
if(grepl("Tissue",SPLITMOD)){
  if(grepl("PM",SPLITMOD)){
    tofill_counts_df <- expand.grid(Tissue = c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                             "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                    CNA_type = c("Loss","Gain"),
                                    Stage = c("Primary","Metastasis"))
  }else if(grepl("PWOWM",SPLITMOD)){
    tofill_counts_df <- expand.grid(Tissue = c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                             "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                    CNA_type = c("Loss","Gain"),
                                    Stage = c("Primary_WO_Metastasis","Primary_W_Metastasis","Metastasis"))
  }
}else if(grepl("Subtype",SPLITMOD)){
  tofill_counts_df <-
    merge(merge(subtypes,
                data.frame(CNA_type=c("Loss","Gain"))),
          data.frame(Stage=c("Primary","Metastasis")))
}
count_real_significants <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(cut){
  sig_results <- filter(real_results,P_value<=cut)
  sig_results$SIG <- ifelse(sig_results$CNA_type=="Gain"&sig_results$Estimate>0,T,
                            ifelse(sig_results$CNA_type=="Loss"&sig_results$Estimate<0,T,F))
  sig_results <- filter(sig_results,SIG)
  if(grepl("Tissue",SPLITMOD)){
    sig_counts <- sig_results %>%
      group_by(Tissue,Stage,CNA_type) %>%
      summarise(Real_Sigs_Sub=n())
    sig_counts <- sig_counts %>%
      group_by(Stage,CNA_type) %>%
      mutate(Real_Sigs=sum(Real_Sigs_Sub,na.rm=T))
  }else if(grepl("Subtype",SPLITMOD)){
    sig_counts <- sig_results %>%
      group_by(Tissue,Subtype,Stage,CNA_type) %>%
      summarise(Real_Sigs_Sub=n())
    sig_counts <- sig_counts %>%
      group_by(Tissue,Stage,CNA_type) %>%
      mutate(Real_Sigs=sum(Real_Sigs_Sub,na.rm=T))
    }
  sig_counts <- merge(tofill_counts_df,sig_counts,all.x=T)
  sig_counts[is.na(sig_counts)] <- 0
  sig_counts$P_cut <- cut
  if(cut==0.0502){
    print("25%")
  }else if(cut==0.1004){
    print("50%")
  }else if(cut==0.1506){
    print("75%")
  }else if(cut==1){
    print("100%")
  }
  return(sig_counts)
})
count_real_significants <- bind_rows(count_real_significants)
if(grepl("Tissue",SPLITMOD)){
  if(grepl("PM",SPLITMOD)){
    print("Number of rows should be: 8032000")
  }else if(grepl("PWOWM",SPLITMOD)){
    print("Number of rows should be: 12048000")
  }
}else if(grepl("Subtype",SPLITMOD)){
  print("Number of rows should be: 40160000")
}
print(nrow(count_real_significants))
write.table(count_real_significants, sprintf("./DATA/ANALYSIS_DATA/2wOG_permutation_real-counts_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)





### 4. COUNTING PERMUTATION SIGNIFICANTS #####
### This part is run in parallel 1:100 times #####
### ... Uploading files ----
permuted_results_all <- readRDS(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_outputs_%s_%s.RDS", FREQ, SPLITMOD))
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
subtypes <- setNames(unique(clinical_data[c("CANC_TYPE","CANC_SUBTYPE")]),
                     c("Tissue","Subtype"))
permutation <- args[2]
print(permutation)

if (!file.exists('./DATA/MODEL_OUTPUTS/')){
  dir.create('./DATA/MODEL_OUTPUTS/')
}
if (!file.exists('./DATA/MODEL_OUTPUTS/2way__OG/temp_permutation')){
  dir.create('./DATA/MODEL_OUTPUTS/2way__OG/temp_permutation')
}

### ... Counting number of significants at different P-value cutoffs, per cancer type, stage and CNA type ----
if(grepl("Tissue",SPLITMOD)){
  if(grepl("PM",SPLITMOD)){
    tofill_counts_df <- expand.grid(Tissue = c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                             "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                    CNA_type = c("Loss","Gain"),
                                    Stage = c("Primary","Metastasis"))
  }else if(grepl("PWOWM",SPLITMOD)){
    tofill_counts_df <- expand.grid(Tissue = c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                             "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                    CNA_type = c("Loss","Gain"),
                                    Stage = c("Primary_WO_Metastasis","Primary_W_Metastasis","Metastasis"))
  }
}else if(grepl("Subtype",SPLITMOD)){
  tofill_counts_df <-
    merge(merge(subtypes,
                data.frame(CNA_type=c("Loss","Gain"))),
          data.frame(Stage=c("Primary","Metastasis")))
}

permuted_results <- permuted_results_all[permutation][[1]]
permutation <- unique(permuted_results$Permutation)
significant_pairs_table <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(fdr_cut){
  significant_results <- filter(permuted_results,as.numeric(P_value)<=fdr_cut)
  significant_results$SIG <- ifelse(significant_results$CNA_type=="Gain"&significant_results$Estimate>0,T,
                                    ifelse(significant_results$CNA_type=="Loss"&significant_results$Estimate<0,T,F))
  significant_results <- filter(significant_results,SIG)
  if(grepl("Tissue",SPLITMOD)){
    sig_counts <- significant_results %>%
      group_by(Tissue,Stage,CNA_type) %>%
      summarise(Permut_Sigs_Sub=n())
    sig_counts <- sig_counts %>%
      group_by(Stage,CNA_type) %>%
      mutate(Permut_Sigs=sum(Permut_Sigs_Sub,na.rm=T))
  }else if(grepl("Subtype",SPLITMOD)){
    sig_counts <- significant_results %>%
      group_by(Tissue,Subtype,Stage,CNA_type) %>%
      summarise(Permut_Sigs_Sub=n())
    sig_counts <- sig_counts %>%
      group_by(Tissue,Stage,CNA_type) %>%
      mutate(Permut_Sigs=sum(Permut_Sigs_Sub,na.rm=T))
  }
  sig_counts <- merge(tofill_counts_df,sig_counts,all.x=T)
  sig_counts[is.na(sig_counts)] <- 0
  sig_counts$P_cut <- fdr_cut
  if(fdr_cut==0.0502){
    print("25%")
  }else if(fdr_cut==0.1004){
    print("50%")
  }else if(fdr_cut==0.1506){
    print("75%")
  }else if(fdr_cut==1){
    print("100%")
  }
  return(sig_counts)
})
significant_pairs_table <- bind_rows(significant_pairs_table)
significant_pairs_table$Permutation <- permutation
if(grepl("Tissue",SPLITMOD)){
  if(grepl("PM",SPLITMOD)){
    print("Number of rows should be: 8032000")
  }else if(grepl("PWOWM",SPLITMOD)){
    print("Number of rows should be: 12048000")
  }
}else if(grepl("Subtype",SPLITMOD)){
  print("Number of rows should be: 40160000")
}

write.table(significant_pairs_table, sprintf("./DATA/ANALYSIS_DATA/2way__OG/temp_permutation/2wOG_permutation_sig-counts_%s_%s_%s.tsv",
                                             FREQ, SPLITMOD, permutation),
            sep="\t", quote=F, row.names=F)



### ... When all 100 sig-counts file are created ----
count_real_significants <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2wOG_permutation_real-counts_%s_%s.tsv", FREQ, SPLITMOD),
                                        sep="\t", header=T)
files <- list.files(path = "./DATA/ANALYSIS_DATA/2way__OG/temp_permutation/",
                    pattern = "2wOG_permutation_sig-counts_*")
print("The number below should be 100")
print(length(files))

if(grepl("Tissue",SPLITMOD)){
  count_real_significants <- unique(count_real_significants[c("Stage","CNA_type","Real_Sigs","P_cut")])
}else if(grepl("Subtype",SPLITMOD)){
  count_real_significants <- unique(count_real_significants[c("Stage","CNA_type","Tissue","Real_Sigs","P_cut")])
}
count_file_list <- lapply(files,function(name){
  count_file <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__OG/temp_permutation/%s", name),
                           header=T, sep="\t")
  if(grepl("Tissue",SPLITMOD)){
    if(grepl("PM",SPLITMOD)){
      print("Number of rows should be: 8032000")
    }else if(grepl("PWOWM",SPLITMOD)){
      print("Number of rows should be: 12048000")
    }
  }else if(grepl("Subtype",SPLITMOD)){
    print("Number of rows should be: 40160000")
  }
  print(nrow(count_file))
  if(grepl("Tissue",SPLITMOD)){
    count_file <- count_file[c("Stage","CNA_type","Real_Sigs","P_cut","Permutation")]
    count_file <- unique(count_file)
    if(grepl("PM",SPLITMOD)){
      print("Number of rows should be: 803200")
    }else if(grepl("PWOWM",SPLITMOD)){
      print("Number of rows should be: 1204800")
    }
  }else if(grepl("Subtype",SPLITMOD)){
    count_file <- count_file[c("Stage","CNA_type","Tissue","Real_Sigs","P_cut","Permutation")]
    count_file <- unique(count_file)
    print("Number of rows should be: 8032000")
  }
  print(nrow(count_file))
  count_file <- merge(count_file,
                      count_real_significants)
  return(count_file)
})
print("all files read and merged with real counts")
count_files <- bind_rows(count_file_list)
print("all files bound into one")
if(grepl("Tissue",SPLITMOD)){
  if(grepl("PM",SPLITMOD)){
    print("Number of rows should be: 80320000")
  }else if(grepl("PWOWM",SPLITMOD)){
    print("Number of rows should be: 120480000")
  }
}else if(grepl("Subtype",SPLITMOD)){
  print("Number of rows should be: 803200000")
}
print(nrow(count_files))


### ... Calculating FDR ----
# FDR is the proportion of false positives by real positives
count_files$FDR <- count_files$Permut_Sigs / count_files$Real_Sigs
count_files$FDR[is.na(count_files$FDR)] <- NA
count_files$FDR[is.infinite(count_files$FDR)] <- NA
write.table(count_files, sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_all-counts-with-FDR_%s_%s.tsv", FREQ,SPLITMOD),
            sep="\t", quote=F, row.names=F)
saveRDS(count_files, sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_all-counts-with-FDR_%s_%s.RDS", FREQ, SPLITMOD))
print("saved permutation_all-counts-with-FDR")
#  Adding statistical information
if(grepl("Tissue",SPLITMOD)){
  count_files <- count_files %>%
    group_by(CNA_type,Stage,P_cut) %>%
    summarise(Mean_FDR=mean(FDR,na.rm=T),
           St_Err_FDR=plotrix::std.error(FDR,na.rm=T),
           Median_FDR=median(FDR,na.rm=T),
           St_Dev_FDR=sd(FDR,na.rm=T))
}else if(grepl("Subtype",SPLITMOD)){
  count_files <- count_files %>%
    group_by(CNA_type,Tissue,Stage,P_cut) %>%
    summarise(Mean_FDR=mean(FDR,na.rm=T),
           St_Err_FDR=plotrix::std.error(FDR,na.rm=T),
           Median_FDR=median(FDR,na.rm=T),
           St_Dev_FDR=sd(FDR,na.rm=T))
}
write.table(count_files, sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_FDR-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)



#### 5. Getting final FDR conversion table ----
permuted_significants_total <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_permutation_FDR-table_%s_%s.tsv", FREQ, SPLITMOD),
                                          sep="\t", header=T)
# # ### ... Comparing number of significants of BH-FDR vs Perm-FDR ----
FDR_conversion_table <- unique(permuted_significants_total[c("CNA_type","Stage","P_cut","Mean_FDR")])
write.table(FDR_conversion_table, sprintf("./DATA/ANALYSIS_DATA/2way__OG/2wOG_FDR-conversion-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)

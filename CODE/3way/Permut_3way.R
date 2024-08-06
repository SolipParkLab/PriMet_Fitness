
### ... PERMUTATION SCRIPT TO ASSESS SIGNIFICANCE OF INTERACTIONS ----
library(vegan)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
options(dplyr.summarise.inform = FALSE)



## Functions/Variables ----
### ... Changeable
n_permutations <- 100
SPLITMOD <- "Tissue-Stage-PM"
# SPLITMOD <- args[3]
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FDR_2way <- 10
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
###### WHICH FDR CORRECTION METHOD ARE WE USING?? ######
FDR_method <- "-PERM"
## Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



##### 1. RUNNING 3WAY MODEL WITH TWO-GENES #####
### ... Loading files ----
# Permuted matrices
permuted_matrixes <- readRDS(sprintf("./DATA/PROCESSED_DATA/PERMUTED_matrixes_%s.RDS",SPLITMOD))
# Glm_outputs from 3-way
real_results_twogenes <- readRDS(sprintf("./DATA/GLM_OUTPUTS/3way/3way_glm-outputs_%s_%s.RDS", FREQ, SPLITMOD))


### ... Creating function to generate two-genes binary matrix ONLY FOR SIGNIFICANT 2WAY GENES ----
more_pseudo_generating_twogenes_binary_mat <- function(tiss_mat, genes_a, genes_b){
  gene_a <- genes_a
  ##  Iterating over every Gene A
  twogenes_list <- bind_rows(lapply(genes_b,function(gene_b){
    ##  Matrix with only Gene A and Gene B information
    two_genes_bm <- tiss_mat[c(paste0(gene_a,"_mutation"),paste0(gene_a,"_Gain"),paste0(gene_a,"_Loss"),
                               paste0(gene_b,"_mutation"))]
    names(two_genes_bm) <- c("GA_mut","GA_Gain","GA_Loss","GB_mut")
    ##  Calculating alteration frequencies
    merged_gene_mat <- two_genes_bm %>% group_by(GB_mut) %>% summarise(Size = n(),
                                                                       FreqMutA = mean(GA_mut),
                                                                       FreqGainA = mean(GA_Gain),
                                                                       FreqLossA = mean(GA_Loss))
    muta <- mean(two_genes_bm$GA_mut)
    mutb <- mean(two_genes_bm$GB_mut)
    mutgain <- mean(two_genes_bm$GA_Gain)
    mutloss <- mean(two_genes_bm$GA_Loss)
    ##  Merging CNA information and turning loss into -1, wt into 0, gain into 1
    two_genes_bm$cnvs <- paste0(two_genes_bm$GA_Loss,two_genes_bm$GA_Gain)
    counts <- paste0(two_genes_bm$GA_mut,two_genes_bm$cnvs,two_genes_bm$GB_mut)
    counts <- c(counts,rep("1100",length(counts[counts=="1110"])))
    counts[counts=="1110"] <- "1010"
    counts <- c(counts,rep("0100",length(counts[counts=="0110"])))
    counts[counts=="0110"] <- "0010"
    counts <- c(counts,rep("0101",length(counts[counts=="0111"])))
    counts[counts=="0111"] <- "0011"
    counts <- c(counts,rep("1101",length(counts[counts=="1111"])))
    counts[counts=="1111"] <- "1011"
    ##  Calculating occurrences (mutA*cnv*mutB)
    possib_occur <- c("1011","0011","1101","0101","1001","0001",
                      "1010","0010","1100","0100","1000","0000")
    ##  Count of each occurrence
    N <- unlist(lapply(seq_along(possib_occur),function(x) sum(counts==possib_occur[x])))
    ##  Adding pseudocounts
    N <- N + 1
    N_df <- data.frame(matrix(N,
                              ncol=6,
                              nrow=2,
                              byrow=T))
    names(N_df) <-c("MutALoss","NoMutALoss","MutAGain","NoMutAGain","MutAWT","NoMutAWT")
    N_df$GB_mut <- c(1,0)
    input_row <- merge(merged_gene_mat,N_df,all=T)
    input_row <- bind_cols(input_row,
                           data.frame(GeneA = gene_a,
                                      GeneB = gene_b,
                                      Gene_Pair = paste0(gene_a,"_",gene_b),
                                      FreqMutA_all = muta,
                                      FreqMutB_all = mutb,
                                      FreqGainA_all = mutgain,
                                      FreqLossA_all = mutloss))
    return(input_row)}))
  return(twogenes_list)}



### ... Retrieving only significant 2way pairs from permuted results ----
real_results_twogenes$Model <- paste0(real_results_twogenes$Tissue,".NA.",real_results_twogenes$Stage)
real_results_twogenes_list <- split(real_results_twogenes,real_results_twogenes$Model)



### ... Generating two-genes inputs with sig 2way pairs, running 3way regression model ----
dir.create('./DATA/ANALYSIS_DATA/3way/perm_outputs/')
twogenes_results <- mapply(function(bm_list,model_sigs){
  mod <- as.numeric(args[1])
  model_sigs <- real_results_twogenes_list[[mod]]
  bm_list <- permuted_matrixes[[mod]]
  mod_name <- unique(model_sigs$Model)
  permuted_models_twogenes <- lapply(seq_along(bm_list),function(permut){
    bmat <- bm_list[[permut]]
    tested_genes <- split(model_sigs[c("GeneB")],list(model_sigs$GeneA,model_sigs$CNA_type))
    tested_genes <- tested_genes[lapply(tested_genes,nrow)!=0]
    twogenes_counts <- lapply(seq_along(tested_genes),function(index){
      set <- names(tested_genes)[index]
      gene_a <- str_split_fixed(set,"\\.",2)[,1]
      cna <- str_split_fixed(set,"\\.",2)[,2]
      genes_b <- tested_genes[[index]][[1]]
      twogenes_pair_counts <- more_pseudo_generating_twogenes_binary_mat(bmat, genes_b = genes_b, genes_a = gene_a)
      twogenes_pair_counts <- filter(twogenes_pair_counts, FreqMutB_all >= freqmut_threshold)
      results_2way <- twogenes_retrieve_2way_results(twogenes_pair_counts, cna)
      wide_twogenes_inputs <- pivot_wider(twogenes_pair_counts,names_from=GB_mut,values_from=2:11)
      names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),"_1","MutB")
      names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),"_0","NoMutB")
      if (cna=="Loss"){
        notcna <- "Gain"
        dataf <-  data.frame(mutA=c(1,0,1,0,1,0,1,0),
                             cnvA=c(-1,-1,0,0,-1,-1,0,0),
                             mutB=c(1,1,1,1,0,0,0,0))
      }else if (cna=="Gain"){
        notcna <- "Loss"
        dataf <-  data.frame(mutA=c(1,0,1,0,1,0,1,0),
                             cnvA=c(1,1,0,0,1,1,0,0),
                             mutB=c(1,1,1,1,0,0,0,0))
      }
      wide_twogenes_inputs <- wide_twogenes_inputs[-grep(paste0("A",notcna),names(wide_twogenes_inputs))]
      names(wide_twogenes_inputs) <- str_replace_all(names(wide_twogenes_inputs),paste0("A",cna),"ACNV")
      wide_twogenes_inputs$CNA_type <- cna
      reg_model <- as.data.frame(t(apply(wide_twogenes_inputs,1,function(datarow){
        N <- c(datarow["MutACNVMutB"],datarow["NoMutACNVMutB"],datarow["MutAWTMutB"],datarow["NoMutAWTMutB"],
               datarow["MutACNVNoMutB"],datarow["NoMutACNVNoMutB"],datarow["MutAWTNoMutB"],datarow["NoMutAWTNoMutB"])
        dataf$N <- as.numeric(N)
        result <- reg_model_3w_twog(dataf)
        return(result)
      },simplify=T)))
      colnames(reg_model) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance")
      wide_twogenes_inputs <- bind_cols(wide_twogenes_inputs,reg_model)
      wide_twogenes_inputs <- bind_cols(wide_twogenes_inputs,results_2way)
      return(wide_twogenes_inputs)
    })
    twogenes_counts <- bind_rows(twogenes_counts)
    twogenes_counts$Model <- unique(model_sigs$Model)
    twogenes_counts[c("Tissue","Subtype","Stage")] <- str_split_fixed(twogenes_counts$Model,"\\.",3)
    twogenes_counts$Permutation <- permut
    return(twogenes_counts)
  })
  permuted_models_twogenes <- map_df(permuted_models_twogenes,~as.data.frame(.x),.id="Permutation")
  saveRDS(permuted_models_twogenes,
          sprintf("./DATA/ANALYSIS_DATA/3way/perm_outputs/3way_permutation_outputs_%s_%s_%s.RDS", FREQ, SPLITMOD, mod_name))
  return(permuted_models_twogenes)
},
permuted_matrixes,
real_results_twogenes_list,
SIMPLIFY=F)



### ... Merging all files together ----
files <- list.files(path=sprintf("./DATA/ANALYSIS_DATA/3way/perm_outputs/",FREQ,SPLITMOD),
                    pattern="3way_permutation_outputs*")
names <- str_remove_all(str_remove_all(files,"3way_permutation_outputs"),".RDS")
twogenes_results <- lapply(files,function(names){return(readRDS(sprintf("./DATA/ANALYSIS_DATA/3way/%s", names)))})
names(twogenes_results) <- names



### ... Saving permutation outputs ----
saveRDS(twogenes_results,
        sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_outputs_%s_%s.RDS", FREQ, SPLITMOD))



##### 2. COUNTING SIGNIFICANTS #####
### ... Loading files ----
twogenes_results <- readRDS(sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_outputs_%s_%s.RDS", FREQ, SPLITMOD))
twogenes_results <- unlist(lapply(twogenes_results,function(tissue){
  permut_tissue <- split(tissue,tissue$Permutation)
  return(permut_tissue)
}),recursive=F)
twogenes_results <- lapply(1:n_permutations,function(perm){
  permutation_results <- twogenes_results[seq(perm,max(length(twogenes_results)),by=100)]
  permutation_results <- bind_rows(permutation_results)
  return(permutation_results)
})



### Counting significants ----
tofill_counts_df <- expand.grid(Tissue=c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary","Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                CNA_type=c("Loss","Gain"),
                                Stage=c("Primary","Metastasis"))
p_val_significant_twogenes_counts_list <- lapply(seq_along(twogenes_results)[as.numeric(args[1])],function(twog_mat_index){
  twog_mat <- twogenes_results[[twog_mat_index]]
  if(nrow(twog_mat)==0){
    twog_mat <- data.frame(P_value=NA)
  }
  counts_per_p_cut <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(fdr_cut){
    twog_sigs <- filter(twog_mat,P_value<=fdr_cut)
    count_table <- reshape2::melt(table(twog_sigs$Tissue,twog_sigs$Stage,twog_sigs$CNA_type),varnames=c("Tissue","Stage","CNA_type"),value.name="Permut_Sig_twog_pairs")
    sig_counts <- merge(tofill_counts_df,count_table,all.x=T)
    sig_counts[is.na(sig_counts)] <- 0
    sig_counts$P_cut <- fdr_cut
    return(sig_counts)
  })
  counts_per_p_cut <- bind_rows(counts_per_p_cut)
  return(counts_per_p_cut)})
names(p_val_significant_twogenes_counts_list) <- as.numeric(args[1]):as.numeric(args[2])
p_val_significant_twogenes_counts_list <- map_df(p_val_significant_twogenes_counts_list,~as.data.frame(.x), .id="Permutation")
dir.create("./DATA/ANALYSIS_DATA/3way/sig_counts/")
write.table(p_val_significant_twogenes_counts_list,
            sprintf("./DATA/ANALYSIS_DATA/3way/sig_counts/3way_permutation-fromPERM_sig-pairs_%s_%s_%s.tsv", FREQ, SPLITMOD, args[1]),
            sep="\t", quote=F, row.names=F)



files <- list.files(path="./DATA/ANALYSIS_DATA/3way/sig_counts/",
                    pattern="3way_permutation-fromPERM_sig-pairs")
res <- lapply(files,function(name){
  return(read.delim(sprintf("./DATA/ANALYSIS_DATA/sig_counts/3way/%s", name),
                    header=T, sep="\t"))
})
p_val_significant_twogenes_counts_list <- map_df(res,~as.data.frame(.x),.id="Permutation")
write.table(p_val_significant_twogenes_counts_list,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_bindrows-sig-pairs_%s_%s_FDR%s.tsv", FREQ, SPLITMOD, FDR_2way),
            sep="\t", quote=F, row.names=F)




p_val_significant_twogenes_counts_list <- read.delim(sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_bindrows-sig-pairs_%s_%s_FDR%s.tsv", FREQ, SPLITMOD, FDR_2way),
                                                     sep="\t", header=T)


## Counting real model significants ----
real_results_twogenes <- readRDS(sprintf("./DATA/GLM_OUTPUTS/3way/3way_glm-outputs_%s_%s.RDS", FREQ, SPLITMOD))
tofill_counts_df <- expand.grid(Tissue=c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary","Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                CNA_type=c("Loss","Gain"),
                                Stage=c("Primary","Metastasis"))
count_real_significants_twogenes <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(cut){
  sig_results <- filter(real_results_twogenes,P_value<=cut)
  sig_counts <- reshape2::melt(table(sig_results$Tissue,sig_results$Stage,sig_results$CNA_type),varnames=c("Tissue","Stage","CNA_type"),value.name="Real_Sig_twog_pairs")
  sig_counts <- merge(tofill_counts_df,sig_counts,all.x=T)
  sig_counts[is.na(sig_counts)] <- 0
  sig_counts$P_cut <- cut
  return(sig_counts)
})
count_real_significants_twogenes <- bind_rows(count_real_significants_twogenes)
write.table(count_real_significants_twogenes,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_real-sig-counts_%s_%s_FDR%s.tsv", FREQ, SPLITMOD, FDR_2way),
            sep="\t", quote=F, row.names=F)



## Merging count tables (permutation + real results) ----
count_real_significants_twogenes <- read.delim(sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM_real-sig-counts_%s_%s_FDR%s.tsv", FREQ, SPLITMOD, FDR_2way),
                                               sep="\t", header=T)
permuted_twogenes_counts_bytissue <- merge(count_real_significants_twogenes,
                                           p_val_significant_twogenes_counts_list,
                                           by=c("Tissue","CNA_type","Stage","P_cut"))



### ... FDR is the proportion of false positives by real positives
permuted_twogenes_counts_bytissue$FDR <- permuted_twogenes_counts_bytissue$Permut_Sig_twog_pairs/permuted_twogenes_counts_bytissue$Real_Sig_twog_pairs
permuted_twogenes_counts_bytissue$FDR[is.na(permuted_twogenes_counts_bytissue$FDR)] <- 0
permuted_twogenes_counts_bytissue$FDR[is.infinite(permuted_twogenes_counts_bytissue$FDR)] <- 0
### ... Adding statistical information
permuted_twogenes_counts_bytissue <- permuted_twogenes_counts_bytissue %>%
  group_by(CNA_type,Stage,P_cut) %>%
  mutate(Mean_FDR=mean(FDR,na.rm=T),
         St_Err_FDR=plotrix::std.error(FDR,na.rm=T),
         Median_FDR=median(FDR,na.rm=T))
### ... Overwriting previous count file with joined infos (permut + real) + FDR
write.table(permuted_twogenes_counts_bytissue,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM-0s_FDR-table_%s_%s.tsv",FREQ,SPLITMOD,FREQ,SPLITMOD),
            sep="\t", quote=F, row.names=F)



### ... Loading files ----
permuted_twogenes_counts <- read.delim(sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM-0s_FDR-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", header=T)
real_results_twogenes <- readRDS(sprintf("./DATA/ANALYSIS_DATA/3way/3way_PERM_analysis_%s_%s_FDR%s.RDS", FREQ, SPLITMOD, FDR_2way))

permuted_twogenes_counts_conflictive <- filter(permuted_twogenes_counts,Real_Sig_twog_pairs==0)
permuted_twogenes_counts_ok <- filter(permuted_twogenes_counts,Real_Sig_twog_pairs!=0)

permuted_twogenes_counts_conflictive$FDR <- permuted_twogenes_counts_conflictive$Permut_Sig_twog_pairs/permuted_twogenes_counts_conflictive$Real_Sig_twog_pairs
permuted_twogenes_counts_conflictive$FDR[is.na(permuted_twogenes_counts_conflictive$FDR)] <- 0
permuted_twogenes_counts_conflictive$FDR[is.infinite(permuted_twogenes_counts_conflictive$FDR)] <- 0
#  Adding statistical information
permuted_twogenes_counts_conflictive <- permuted_twogenes_counts_conflictive %>%
  group_by(CNA_type,Stage,P_cut) %>%
  mutate(Mean_FDR=mean(FDR,na.rm=T),
         St_Err_FDR=plotrix::std.error(FDR,na.rm=T),
         Median_FDR=median(FDR,na.rm=T))
permuted_twogenes_counts <- bind_rows(permuted_twogenes_counts_conflictive,permuted_twogenes_counts_ok)
write.table(permuted_twogenes_counts,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM-0s_FDR-table_%s_%s.tsv",FREQ,SPLITMOD,FREQ,SPLITMOD),
            sep="\t", quote=F, row.names=F)

permuted_twogenes_counts <- unique(permuted_twogenes_counts[c("CNA_type","Stage","P_cut","Mean_FDR","St_Err_FDR")])
write.table(permuted_twogenes_counts,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_permutation-fromPERM-0s_conversion-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)



### ... Making final conversion table ----
FDR_conversion_table <- read.delim(sprintf("./DATA/ANALYSIS_DATA/3way/3way_FDR-conversion-table_%s_%s.tsv", FREQ, SPLITMOD),
								   sep="\t", header=T)
FDR_conversion_table_all <- FDR_conversion_table %>% 
	group_by(CNA_type,Stage,P_cut) %>%
	summarise(Mean_FDR=mean(Mean_FDR,na.rm=T),
			  St_Err_FDR=plotrix::std.error(Mean_FDR,na.rm=T))
write.table(FDR_conversion_table_all,
            sprintf("./DATA/ANALYSIS_DATA/3way/3way_FDR-conversion-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)

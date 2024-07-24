
### ... PERMUTATION SCRIPT TO ASSESS SIGNIFICANCE OF INTERACTIONS ----
library(vegan)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(reshape2)
args <- commandArgs(trailingOnly = TRUE)
options(dplyr.summarise.inform = FALSE)




## Functions/Variables ----
### ... Changeable ---
n_permutations <- 100
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FDR_2way <- 10 
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)
### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Uploading files ----
binary_mats <- readRDS("./DATA/PROCESSED_DATA/001_oncokb_binary_matrixes.RDS")
clinical_data <- read.table("./DATA/PROCESSED_DATA/p_clinical-data.tsv",
                            sep = "\t", header=T)
cancgenedf <- read.csv("./DATA/PROCESSED_DATA/p_cancer-gene-list.tsv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)



### ... Subsetting clinical data according to model ----
clinical_data <- filter(clinical_data,N_ONCOGENIC_ALTERATIONS>0)
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



### ... Running permutations in parallel function ----
binary_mat_index <- as.numeric(args[1])
alter <- args[2]
binary_mat <- ready_bm[[binary_mat_index]]
binary_mat_name <- names(ready_bm)[binary_mat_index]
binary_alter <- binary_mat[grep(alter,names(binary_mat))]
permutation_matrixes_alter <- permatswap(binary_alter,times=n_permutations,method = "quasiswap")
dir.create('./DATA/PROCESSED_DATA/treatment_matrices/')
saveRDS(permutation_matrixes_alter$perm,
        sprintf("./DATA/PROCESSED_DATA/treatment_matrices/2wT-PERMUT_matrixes_%s%s.RDS",binary_mat_name,alter))


# When all of them have been done
files <- list.files(path="./DATA/PROCESSED_DATA/treatment_matrices/",
                    pattern="2wT-PERMUT_matrixes")

permutation_matrixes <- lapply(names(ready_bm),function(model){
  matrixes <- lapply(files[grep(model,files)],function(permutation_files){
    return(readRDS(sprintf("./002_processed_data/%s",permutation_files)))
  })
  merged_matrixes <- lapply(1:n_permutations,function(permut){
    bigmat <- merge(matrixes[[1]][permut],matrixes[[2]][permut],by="row.names")
    rownames(bigmat) <- bigmat$Row.names
    bigmat$Row.names <- NULL
    bigmat <- merge(bigmat,matrixes[[3]][permut],by="row.names")
    rownames(bigmat) <- bigmat$Row.names
    bigmat$Row.names <- NULL
    return(bigmat)
    })
  return(merged_matrixes)
})
names(permutation_matrixes) <- names(ready_bm)
saveRDS(permutation_matrixes,
        sprintf("./DATA/PROCESSED_DATA/PERMUTED_matrixes_Treatment-%s.RDS",SPLITMOD))



##### 2. RUNNING 2WAY REGRESSION MODEL OVER MATRICES #####
### ... Loading files ----
permuted_matrixes <- readRDS(sprintf("./DATA/PROCESSED_DATA/PERMUTED_matrixes_Treatment-%s.RDS",SPLITMOD))
### ... Function to create the input table from PERMUTED MATRIX ----
# In this matrix we don't check if gain / loss events in the same gene happen
# We simply count them both to maintain alteration frequencies
pseudo_oncogenic_model_input <- function(df,freqmut_threshold,freqcnv_threshold,pairs){
  if(is_tibble(df)){df <- as.data.frame(df)}
  size <- nrow(df)
  inputs_table <- mapply(function(gene,cna){
    gene_column_names <- paste0(gene,c("_mutation","_Loss","_Gain"))
    gene_bm <- df[gene_column_names]
    freqs <- apply(gene_bm,2,mean)
    reshape2::melt(table(gene_bm))
    gene_bm$cnvs <- paste0(gene_bm[,2],gene_bm[,3])
    counts <- paste0(gene_bm[,1],gene_bm$cnvs)
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
    pairs$Gene,pairs$CNA_type,SIMPLIFY=F)
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


# args <- c("1") # I think this is for a test
### ... Running 2way regression model ----
###### IF WE ARE RUNNING ONLY IN TESTED REAL PAIRS
real_results <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_PERM_analysis_%s_%s.tsv",FREQ,SPLITMOD),
                           sep="\t", header=T)
if(grepl("Tissue",SPLITMOD)){
  real_results$model <- paste0(real_results$Tissue,".NA.",real_results$Stage,".",real_results$Treatment)
}else if(grepl("Subtype",SPLITMOD)){
  real_results$model <- paste0(real_results$Tissue,".",real_results$Subtype,".",real_results$Stage,".",real_results$Treatment)
}
permuted_matrix_index <- as.numeric(args[1])
permuted_results <- lapply(seq_along(permuted_matrixes),function(permuted_matrix_index){
  permuted_matrixes_list <- permuted_matrixes[[permuted_matrix_index]]
  permuted_matrix_name <- names(permuted_matrixes)[permuted_matrix_index]
  real_tested_pairs <- filter(real_results,model==permuted_matrix_name)[c("Gene","CNA_type")]
  permuted_inputs <- lapply(permuted_matrixes_list, pseudo_oncogenic_model_input, freqmut_threshold, freqcnv_threshold, real_tested_pairs)
  permuted_outputs <- lapply(permuted_inputs, pseudo_oncogenic_retrieve_results)
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
names(permuted_results) <- names(permuted_matrixes)
saveRDS(permuted_outputs,
        sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT-PERMUT_matrix_outputs_%s.RDS", FREQ, SPLITMOD, permuted_matrix_name))
##  For loading 100 permuted files into a single file
files <- list.files(path="./DATA/ANALYSIS_DATA/2way__T/",
                    pattern="2wT-PERMUT_*")
model_names <- str_replace_all(str_replace_all(files,"2wT-PERMUT_matrix_outputs_",""),".RDS","")
permuted_results <- lapply(files,function(name){
  return(readRDS(sprintf("./DATA/ANALYSIS_DATA/2way__T/%s",name)))
})
names(permuted_results) <- model_names

permuted_results <- unlist(permuted_results,recursive=F)
names(permuted_results) <- unlist(lapply(model_names,function(name){return(rep(name,100))}))
## Joining cancer types across permutations
permuted_results <- lapply(1:n_permutations,function(permutation){
  model_outputs <- map_df(permuted_results[seq(permutation,max(length(permuted_results)),by=100)],~as.data.frame(.x),.id="Tissue")
  model_outputs[c("Tissue","Subtype","Stage","Treatment")] <- str_split_fixed(model_outputs$Tissue,"\\.",4)
  return(model_outputs)})
saveRDS(permuted_results,
        "./DATA/ANALYSIS_DATA/2way__T/2wT_permutation-outputs.RDS")



#### 3. COUNTING SIGNIFICANTS #####
### ... Uploading files ----
permuted_results <- readRDS("./DATA/ANALYSIS_DATA/2wT_permutation-outputs.RDS")


### ... Counting number of significants at different P-value cutoffs, per cancer type, stage and CNA type ----
tofill_counts_df <- expand.grid(Tissue=c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                         "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                CNA_type=c("Loss","Gain"),
                                Stage=c("Primary","Metastasis"),
                                Treatment=c("Treatment","Non_Treatment"))
significant_tables_counts <- lapply(seq_along(permuted_results)[as.numeric(args[1]):as.numeric(args[2])],function(m_ind){
  print(m_ind)
  m <- permuted_results[[m_ind]]
  new_significant_pairs_tables <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(fdr_cut){
      sig_m <- filter(m, as.numeric(P_value)<=fdr_cut)
    sig_m$SIG <- ifelse(sig_m$CNA_type=="Gain"&sig_m$Estimate>0,T,
                              ifelse(sig_m$CNA_type=="Loss"&sig_m$Estimate<0,T,F))
    sig_m <- filter(sig_m,SIG)
    sig_counts <- reshape2::melt(table(sig_m$Tissue,sig_m$Stage,sig_m$Treatment,sig_m$CNA_type),
                                 varnames=c("Tissue","Stage","Treatment","CNA_type"),
                                 value.name="Permut_Sigs")
    sig_counts <- merge(tofill_counts_df,sig_counts,all.x=T)
    sig_counts[is.na(sig_counts)] <- 0
    sig_counts$P_cut <- fdr_cut
    return(sig_counts)
  })
  new_significant_pairs_tables <- bind_rows(new_significant_pairs_tables)
  return(new_significant_pairs_tables)
})
names(significant_tables_counts) <- as.numeric(args[1]):as.numeric(args[2])
significant_tables_counts <- map_df(significant_tables_counts,~as.data.frame(.x),.id="Permutation")
significant_tables_counts <- significant_tables_counts %>%
  group_by(CNA_type,Stage,Treatment,Tissue,P_cut) %>%
  summarise(Permut_Sigs = mean(Permut_Sigs,na.rm=T))
write.table(significant_tables_counts,
            sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_permutation_sig-counts_%s.tsv",args[2]),
            sep="\t", quote=F, row.names=F)



##  For loading 100 permuted files into a single file
files <- list.files(path="./DATA/ANALYSIS_DATA/2way__T/",
                    pattern="2wT_permutation_sig-counts_counts*")
# files <- files[1:50] # I think this has to be commented for 1 to 100
res <- lapply(files,function(name){
  count_file <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__T/%s", name),
                           header=T, sep="\t")
  return(count_file)
})
resu <- map_df(res,~as.data.frame(.x),.id="Permutation")


# Depending on which grouping we want, comment one of these two options
# With Tissue
resu <- resu %>%
  group_by(CNA_type,Stage,Treatment,Tissue,P_cut) %>%
  summarise(Permut_Sigs = mean(Permut_Sigs,na.rm=T))
# Without Tissue
resu <- resu %>%
  group_by(CNA_type,Stage,Treatment,P_cut) %>%
  summarise(Permut_Sigs = sum(Permut_Sigs,na.rm=T))
saveRDS(resu,
        sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_bindrows-permutation_sig-counts_%s_%s.RDS", FREQ, SPLITMOD))



### ... Counting real model significants ----
real_results <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_PERM_analysis_%s_%s.tsv", FREQ, SPLITMOD),
                           sep="\t", header=T)
tofill_counts_df <- expand.grid(Tissue=c("Breast","Core-GI","Developmental-GI-Tract","Endocrine","Genitourinary",
                                         "Gynecologic","Head-and-Neck","Skin","Soft-Tissue","Thoracic"),
                                CNA_type=c("Loss","Gain"),
                                Stage=c("Primary","Metastasis"),
                                Treatment=c("Treatment","Non_Treatment"))
count_real_significants <- lapply(c(seq(0.000001,0.199999,0.000001),seq(0.2,1,0.001)),function(cut){
  sig_results <- filter(real_results,P_value<=cut)
  sig_results$SIG <- ifelse(sig_results$CNA_type=="Gain"&sig_results$Estimate>0,T,
                            ifelse(sig_results$CNA_type=="Loss"&sig_results$Estimate<0,T,F))
  sig_results <- filter(sig_results,SIG)
  sig_counts <- reshape2::melt(table(sig_results$Tissue,sig_results$Stage,sig_results$Treatment,sig_results$CNA_type),
                               varnames=c("Tissue","Stage","Treatment","CNA_type"),value.name="Real_Sigs")
  sig_counts <- merge(tofill_counts_df,sig_counts,all.x=T)
  sig_counts[is.na(sig_counts)] <- 0
  sig_counts$P_cut <- cut
  return(sig_counts)
})
count_real_significants <- bind_rows(count_real_significants)
write.table(count_real_significants, sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_permutation_real-counts_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)


### ... Merging count tables (permutation + real results) ----
## Uploading counts ----
significant_tables_counts <- readRDS(sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_bindrows-permutation_sig-counts_%s_%s.RDS", FREQ, SPLITMOD))
count_real_significants <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_permutation_real-counts_%s_%s.tsv", FREQ, SPLITMOD),
                                        sep="\t", header=T)

#### AQUÍ SE TENDRA QUE AGRUPAR IGUAL QUE CON LOS FALSOS SIGNIFICATIVOS
count_real_significants <- count_real_significants %>%
  group_by(CNA_type, Stage, Treatment, P_cut) %>%
  summarise(Real_Sigs = sum(Real_Sigs, na.rm=T))

write.table(count_real_significants,
            sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_permutation_summarised-real-counts_%s_%s.tsv",FREQ,SPLITMOD,FREQ,SPLITMOD),
            sep="\t", quote=F, row.names=F)

significant_tables_counts <- significant_tables_counts %>%
  group_by(P_cut, CNA_type, Stage, Treatment, Permutation) %>%
  summarise(Permut_Sigs=sum(Permut_Sigs, na.rm=T))

by_vect <- c("P_cut","CNA_type","Stage","Treatment")
permuted_significants_total <- merge(significant_tables_counts,
                                     count_real_significants,
                                     by = by_vect)

### ... FDR is the proportion of false positives by real positives
permuted_significants_total$FDR <- permuted_significants_total$Permut_Sigs/permuted_significants_total$Real_Sigs
### ... Adding statistical information
permuted_significants_total <- permuted_significants_total %>%
  group_by(CNA_type,Treatment,Stage,P_cut) %>%
  mutate(Mean_FDR = mean(FDR, na.rm=T),
         St_Err_FDR = plotrix::std.error(FDR, na.rm=T),
         Median_FDR = median(FDR, na.rm=T))
write.table(permuted_significants_total,
            sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_bindrows-permutation_FDR-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)



#### 4. PLOTTING THINGS #####
permuted_significants_total <- read.delim(sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_bindrows-permutation_FDR-table_%s_%s.tsv", FREQ, SPLITMOD),
                                          sep="\t", header=T)


### ... Comparing number of significants of BH-FDR vs Perm-FDR
FDR_conversion_table <- unique(permuted_significants_total[c("CNA_type","Stage","Treatment","P_cut","Mean_FDR")])
write.table(FDR_conversion_table,
            sprintf("./DATA/ANALYSIS_DATA/2way__T/2wT_FDR-conversion-table_%s_%s.tsv", FREQ, SPLITMOD),
            sep="\t", quote=F, row.names=F)

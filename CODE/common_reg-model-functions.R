## Commonly used functions for the regression models

## Function to negate "%in%"
"%!in%" <- Negate("%in%")

## Function to pad with NAs
na.pad <- function(x, len){
  x[1:len]
}


## Function to create the input table from binary matrix 
oncogenic_model_input <- function(df,freqmut_threshold,freqcnv_threshold){
  ### ... Changing tibble to df, if needed
  if(is_tibble(df)){df <- as.data.frame(df)}
  ### ... Row number indicates the size (number of pairs)
  size <- nrow(df)
  ### ... Retrieving all gene names (from column names) by eliminating the "_mutation" part
  genes <- str_replace_all(colnames(df)[grep("_mutation",colnames(df))],"_mutation","")
  ### ... Looping over gene names vector
  inputs_table <- lapply(genes,function(gene){
    ### ... Retrieving correspondent columns
    gene_column_names <- paste0(gene,c("_mutation","_Loss","_Gain"))
    gene_bm <- df[gene_column_names]
    ### ... Merging Loss and Gain counts into "cnvs"
    gene_bm$cnvs <- paste0(gene_bm[,2],gene_bm[,3])
    ### ... Transforming loss and gain counts to -1, 0 or 1
    gene_bm$cnvs <- ifelse(gene_bm$cnvs=="10",-1,
                           ifelse(gene_bm$cnvs=="01",1,0))
    ### ... Vector with possible occurrences
    possib_occur <- c("0-1","1-1","00","10","01","11")
    ### ... Count of each occurrence
    N <- unlist(lapply(1:6,function(x) sum(paste0(gene_bm[,1],gene_bm$cnvs)==possib_occur[x])))
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
                                      Gain_filter = gainfreq >= freqcnv_threshold),
                           N_df)
    return(input_row)})
  inputs_table <- bind_rows(inputs_table)
  return(inputs_table)
  }


## Function to generate regression 2way model
reg_model_2w <- function(dataf){
  form <- "N ~ mut*cnv"
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[4,]), 
              summary(mod)$null.deviance, 
              summary(mod)$deviance,
              dataf$N)
  return(result)}


## Function to create result table for the 2way model 
oncogenic_retrieve_results <- function(df) {
  if(is.null(df)){return(NULL)}
  ### ... Filtering by mutation frequency threshold
  mut_df <- filter(df,Mutation_filter==TRUE)
  to_return <- list()
  for (cna in c("Loss","Gain")){
    ### ... Filtering by cnv frequency threshold
    cnv_df <- filter(mut_df,get(paste0(cna,"_filter"))==T)
    if (nrow(cnv_df)==0){
      next
    }else{
      ### ... Running regression model over every gene in cnv_df
      regmod_df <- as.data.frame(t(sapply(cnv_df$Gene, function(gene){
        ### ... Retrieving the correspondent row
        datarow <- cnv_df[cnv_df$Gene==gene,]
        ### ... Creating input dataframe for the regression model
        if (cna=="Loss"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutLoss,datarow$MutLoss)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=N)
        }else if (cna=="Gain"){
          N <- c(datarow$NoMutWT,datarow$MutWT,datarow$NoMutGain,datarow$MutGain)
          dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=N)
        }
        ### ... Running the regression model
        result <- reg_model_2w(dataf)
        return(result)
      })))
      colnames(regmod_df) <- c("Estimate","Std_Error","z_value","P_value","Null_deviance","Residual_deviance",
                               "NoMutWT","MutWT","NoMutCNV","MutCNV")
      regmod_df$Gene <- rownames(regmod_df)
      to_return[[cna]] <- regmod_df
    }
  }
  return(to_return)
}


## Function to retrieve legend from ggplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


## Function to generate two - genes input binary matrix from cancer type binary matrix
generating_twogenes_binary_mat <- function(tiss_mat){
  ### ... List of Genes A
  genes_a <- str_replace_all(names(tiss_mat)[grep("_mutation",names(tiss_mat))],"_mutation","")
  ### ... Iterating over every Gene B
  big_twogenes_mat <- bind_rows(lapply(genes_a,function(gene_a){
    # print(paste0("Calculating matrix for gene ",grep(gene_a,genes_a)," out of ",length(genes_a)))
    ### ... List with Genes A (without Gene B)
    genes_b <- genes_a[-grep(gene_a,genes_a)]
    ### ... Iterating over every Gene A
    twogenes_list <- bind_rows(lapply(genes_b,function(gene_b){
      ### ... Matrix with only Gene A and Gene B information
      two_genes_bm <- tiss_mat[c(paste0(gene_a,"_mutation"),paste0(gene_a,"_Gain"),paste0(gene_a,"_Loss"),
                                 paste0(gene_b,"_mutation"))]
      names(two_genes_bm) <- c("GA_mut","GA_Gain","GA_Loss","GB_mut")
      ### ... Calculating alteration frequencies 
      merged_gene_mat <- two_genes_bm %>% group_by(GB_mut) %>% summarise(Size = n(),
                                                                         FreqMutA = mean(GA_mut),
                                                                         FreqGainA = mean(GA_Gain),
                                                                         FreqLossA = mean(GA_Loss))
      merged_gene_mat <- merge(data.frame(GB_mut=c(0,1)),merged_gene_mat,all.x=T)
      merged_gene_mat[is.na(merged_gene_mat)] <- 0
      muta <- mean(two_genes_bm$GA_mut)
      mutb <- mean(two_genes_bm$GB_mut)
      mutgain <- mean(two_genes_bm$GA_Gain)
      mutloss <- mean(two_genes_bm$GA_Loss)
      ### ... Merging CNA information and turning loss into -1, wt into 0, gain into 1
      two_genes_bm$cnvs <- paste0(two_genes_bm$GA_Loss,two_genes_bm$GA_Gain)
      two_genes_bm$cnvs <- ifelse(two_genes_bm$cnvs=="10",-1,
                                  ifelse(two_genes_bm$cnvs=="01",1,0))
      ### ... Calculating occurrences (mutA*cnv*mutB)
      possib_occur <- c("1-11","0-11","111","011","101","001",
                        "1-10","0-10","110","010","100","000")
      ### ... Count of each occurrence
      N <- unlist(lapply(seq_along(possib_occur),function(x) sum(paste0(two_genes_bm$GA_mut,two_genes_bm$cnvs,two_genes_bm$GB_mut)==possib_occur[x])))
      ### ... Adding pseudocounts
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
                                        FreqMutB = mutb))
      return(input_row)}))
    return(twogenes_list)}))
  return(big_twogenes_mat)}



## Function to generate regression 3way model between TWO GENES MODEL
reg_model_3w_twog <- function(dataf){
  form <- "N ~ mutA*cnvA*mutB"
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[8,]),
              summary(mod)$null.deviance,
              summary(mod)$deviance)
  return(result)}


## Function to return gene frequencies & counts in all tissues
calculate_genefreqs_and_counts <- function(bm,gen){
  gene_binary_cols <- bm[grep(paste0(gen,"_"),names(bm))]
  if(ncol(gene_binary_cols)==0){
    return(
      NULL
    )
  }
  loss <- gene_binary_cols[[grep("_Loss",names(gene_binary_cols))]]
  gain <- gene_binary_cols[[grep("_Gain",names(gene_binary_cols))]]
  mut <- gene_binary_cols[[grep("_mutation",names(gene_binary_cols))]]
  cnvs <- paste0(loss,gain)
  counts <- paste0(mut,cnvs)
  counts <- c(counts,rep("010",sum(counts=="011")))
  counts[which(counts=="011")] <- "010"
  counts <- c(counts,rep("110",sum(counts=="111")))
  counts[which(counts=="111")] <- "110"
  ### ... Calculating occurrences (mutA*cnv*mutB)
  possib_occur <- c("101","001","110","010","100","000")
  ### ... Count of each occurrence
  N <- unlist(lapply(seq_along(possib_occur),function(x) sum(counts==possib_occur[x])))
  ### ... We'll add pseudocounts later
  mutfreq <- sum(mut)/length(mut)
  lossfreq <- sum(loss)/length(cnvs)
  gainfreq <- sum(gain)/length(cnvs)
  return(c(mutfreq,lossfreq,gainfreq,N,length(mut)))
}
### ... Function to merge all significant non-target tissues information (freq,counts)
merge_genefreqs_and_counts <- function(df){
  df$MF <- df$MutFreq*df$Size
  df$LF <- df$LossFreq*df$Size
  df$GF <- df$GainFreq*df$Size
  rest_join <- data.frame(t(colSums(df[c(13:15,4:10)])))
  rest_join[1:3] <- rest_join[1:3]/rest_join$Size
  ### ... Including pseudocounts (+1 PER GENE, NOT PER TISSUE)
  rest_join[4:9] <- rest_join[4:9]+1
  names(rest_join)[1:3] <- c("MutFreq","LossFreq","GainFreq")
  rest_join$Tissue <- paste0(df$Tissue,collapse=", ")
  # rest_join$Subtype <- unique(df$Subtype)
  names(rest_join) <- paste0(names(rest_join),"Rest")
  return(rest_join)
}


## Function to generate regression 3way model for a given PAIR in the TISSUE-SPECIFIC MODEL
reg_model_tissue_3way <- function(dataf){
  form <- "N ~ mut*cnv*tiss"
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[8,]), 
              summary(mod)$null.deviance, 
              summary(mod)$deviance)
  return(result)
}


## Function to create result table for the 3way TISSUE-SPECIFIC MODEL 
tissue_retrieving_results <- function(gene_row) {
  ### ... Retrieving cna type (gain or loss) 
  cna <- gene_row[["CNA_type"]]
  ### ... Retrieving counts
  N_3way <- c(gene_row[["NoMutWTRest"]],gene_row[["MutWTRest"]],gene_row[[paste0("NoMut",cna,"Rest")]],gene_row[[paste0("Mut",cna,"Rest")]],
              gene_row[["NoMutWT"]],gene_row[["MutWT"]],gene_row[[paste0("NoMut",cna)]],gene_row[[paste0("Mut",cna)]])
  dataf_3way <- data.frame(
    mut=c(0,1,0,1,0,1,0,1),
    cnv=c(0,0,1,1,0,0,1,1),
    tiss=c(0,0,0,0,1,1,1,1),
    N=N_3way)
  ### ... Running 3way model
  result_3way <- reg_model_tissue_3way(dataf_3way)
  ### ... Retrieving counts for the target tissue
  N_2way_targ <- c(gene_row[["NoMutWT"]],gene_row[["MutWT"]],gene_row[[paste0("NoMut",cna)]],gene_row[[paste0("Mut",cna)]])
  if(cna=="Gain"){
    cnv <- c(0,0,1,1)
  }else if(cna=="Loss"){
    cnv <- c(0,0,-1,-1)
  }
  dataf_2way_targ <- data.frame(
    mut=c(0,1,0,1),
    cnv=cnv,
    N=N_2way_targ)
  ### ... Running 2way model for the target tissue
  result_2way_targ <- reg_model_2w(dataf_2way_targ)
  ### ... Retrieving counts for the rest of tissues
  N_2way_rest <- c(gene_row[["NoMutWTRest"]],gene_row[["MutWTRest"]],gene_row[[paste0("NoMut",cna,"Rest")]],gene_row[[paste0("Mut",cna,"Rest")]])
  dataf_2way_rest <- data.frame(
    mut=c(0,1,0,1),
    cnv=cnv,
    N=N_2way_rest)
  ### ... Running 2way model for the rest of tissues
  result_2way_rest <- reg_model_2w(dataf_2way_rest)
  result <- c(result_3way,result_2way_targ[1:(length(result_2way_targ)-4)],result_2way_rest[1:(length(result_2way_rest)-4)])
  return(result)
}


## Function to create result table for the TWO GENES MODEL, 2way 
## (2way model with the gene B mutated, 2way model with the gene B not mutated)
twogenes_retrieve_2way_results <- function(df,cna) {
  if(is.null(df)){return(NULL)}
  to_return <- list()
  regmod_df <- as.data.frame(t(sapply(unique(df$Gene_Pair), function(gene_pair){
    datarows <- df[df$Gene_Pair==gene_pair,]
    resultrows <- apply(datarows,1,function(datarow){
      if (cna=="Loss"){
        N <- c(datarow[["NoMutAWT"]],datarow[["MutAWT"]],datarow[["NoMutALoss"]],datarow[["MutALoss"]])
        dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,-1,-1), N=as.numeric(N))
      }else if (cna=="Gain"){
        N <- c(datarow[["NoMutAWT"]],datarow[["MutAWT"]],datarow[["NoMutAGain"]],datarow[["MutAGain"]])
        dataf <- data.frame(mut=c(0,1,0,1),cnv=c(0,0,1,1), N=as.numeric(N))
      }
      result <- reg_model_2w(dataf)
      return(result)
      })
    resultrow <- as.data.frame(t(c(resultrows[1:4],resultrows[1:4,2])))
    colnames(resultrow) <- c("Estimate_2w_NoMutB","Std_Error_2w_NoMutB","z_value_2w_NoMutB","P_value_2w_NoMutB","Estimate_2w_MutB","Std_Error_2w_MutB","z_value_2w_MutB","P_value_2w_MutB")
    return(resultrow)})))
  return(regmod_df)
}


## Function to generate 3WAY STAGE/TISSUE MODEL INPUT
model_input_3way <- function(bm0,bm1,gen){
  gene_binary_cols_0 <- bm0[grep(paste0(gen,"_"),names(bm0))]
  gene_binary_cols_1 <- bm1[grep(paste0(gen,"_"),names(bm1))]
  if(ncol(gene_binary_cols_0)==0|ncol(gene_binary_cols_1)==0){
    return(
      NULL
    )
  }
  counts <- lapply(list(gene_binary_cols_0,gene_binary_cols_1),function(gene_binary_cols){
    loss <- gene_binary_cols[[grep("Loss",names(gene_binary_cols))]]
    gain <- gene_binary_cols[[grep("Gain",names(gene_binary_cols))]]
    cnvs <- paste0(loss,gain)
    mut <- gene_binary_cols[[grep("mutation",names(gene_binary_cols))]]
    cnvs <- replace(cnvs, cnvs=="10",-1)
    cnvs <- replace(cnvs,cnvs=="00",0)
    cnvs <- replace(cnvs,cnvs=="01",1)
    oc_table <- data.frame(mut=c(0,1,0,1,0,1)
                           ,Var2=c(-1,-1,0,0,1,1))
    count_table <- merge(oc_table,
                         melt(table(mut
                                    ,as.numeric(cnvs)))
                         ,all.x=T)
    count_table$value[is.na(count_table$value)] <- 0
    N <- count_table$value
    return(c(N,length(mut)))
  })
  count_row <- setNames(unlist(counts),
                        c("NoMutLoss0","NoMutWT0","NoMutGain0",
                          "MutLoss0","MutWT0","MutGain0",
                          "Size0",
                          "NoMutLoss1","NoMutWT1","NoMutGain1",
                          "MutLoss1","MutWT1","MutGain1",
                          "Size1"))
  count_row["Size"] <- count_row["Size0"]+count_row["Size1"]
  count_row["MutFreq"] <- sum(count_row[grep("\\bMut",names(count_row))])/count_row["Size"]
  count_row["GainFreq"] <- sum(count_row[grep("Gain",names(count_row))])/count_row["Size"]
  count_row["LossFreq"] <- sum(count_row[grep("Loss",names(count_row))])/count_row["Size"]
  count_row[c(1:6,8:13)] <- count_row[c(1:6,8:13)] + 1
  return(count_row)
}


### ... Function to generate regression 3way model
reg_model_3w <- function(dataf){
  form <- "N ~ str*mut*cnv"
  mod <- glm(form, family = poisson(link = "log"), data = dataf)
  result <- c(unname(summary(mod)$coef[8,]), 
              summary(mod)$null.deviance, 
              summary(mod)$deviance)
  return(result)}


## Function to create result table for 3way model
retrieving_results_3w <- function(gene_row) {
  if(as.numeric(gene_row[["MutFreq"]])<freqmut_threshold){
    return(gene_row)
  }
  cna <- gene_row[["CNA_type"]]
  if(as.numeric(gene_row[[paste0(cna,"Freq")]])<freqcnv_threshold){
    return(gene_row)
  }
  if(cna=="Gain"){
    cnv <- c(0,0,0,0,1,1,1,1)
  }else if(cna=="Loss"){
    cnv <- c(0,0,0,0,-1,-1,-1,-1)
  }
  N <- as.numeric(c(gene_row[["NoMutWT0"]],gene_row[["NoMutWT1"]],gene_row[["MutWT0"]],gene_row[["MutWT1"]],
                    gene_row[[paste0("NoMut",cna,"0")]],gene_row[[paste0("NoMut",cna,"1")]],gene_row[[paste0("Mut",cna,"0")]],gene_row[[paste0("Mut",cna,"1")]]))
  dataf <- data.frame(str=c(0,1,0,1,0,1,0,1),mut=c(0,0,1,1,0,0,1,1),cnv=cnv,N=N)
  result <- setNames(reg_model_3w(dataf),
                     c("Estimate_3way","Std_Error_3way","z_value_3way","P_value_3way","Null_deviance_3way","Residual_deviance_3way"))
  return(c(gene_row,result))
}

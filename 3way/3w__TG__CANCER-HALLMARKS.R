
### ... Libraries ----
library(dplyr)
library(tidyr)
library(plotrix)
library(purrr)
library(ggplot2)
source("./CODE/common_reg-model-functions.R",local=T)



### ... Variables ----
FDR <- 'PERM'
FREQ <- 'mf1-cf10'
SPLITMOD <- 'Tissue-Stage-PM'
thr <- 10



### ... Load files ----
input_genes <- readRDS('./DATA/PROCESSED_DATA/492_input_genes.RDS') # Genes in the maf file
hallmarks <- read.table('./DATA/PROCESSED_DATA/CELL_CPTAC_Hallmark.txt',
                        sep = '\t', na.strings = '', header = T)
results <- filter(readRDS(sprintf('./DATA/ANALYSIS_DATA/3way__TG/3wTG_%s_analysis_%s_%s_FDR%s.RDS', FDR, FREQ, SPLITMOD, thr)),
                  SIG_FDR10_3way == T)



### ... Creating a list of FALSE pairs that will be used for permutations ----
false_pairs <- expand_grid(GeneA = input_genes, GeneB = input_genes) %>% 
  mutate(true_A = pmap_chr(., min), true_B = pmap_chr(., max)) %>% # Find the min and max in each row
  unite(check, c(true_A, true_B), remove = FALSE) %>% # Combine them in a "check" variable
  distinct(check, .keep_all = TRUE) %>% # Remove duplicates of the "check" variable
  dplyr::select(true_A, true_B)
false_pairs <- filter(false_pairs, true_A != true_B)
names(false_pairs) <- c('GeneA', 'GeneB')
false_pairs$Gene_Pair <- paste0(false_pairs$GeneA,'_',false_pairs$GeneB)



### ... Merging hallmarks information ----
# Need an extra step for mirror pairs in the same condition
res_corrected <- results[c('Stage','GeneA', 'GeneB')] %>% 
  mutate(true_A = pmap_chr(.[c('GeneA', 'GeneB')], min), true_B = pmap_chr(.[c('GeneA', 'GeneB')], max)) %>% 
  unite(check, c(Stage, true_A, true_B), remove = FALSE) %>% 
  distinct(check, .keep_all = TRUE) %>% 
  dplyr::select(Stage, true_A, true_B)
names(res_corrected) <- c('Stage', 'GeneA', 'GeneB')
res_corrected$Gene_Pair <- paste0(res_corrected$GeneA,'_',res_corrected$GeneB)

# Change to long format, so each gene can be merged to hallmarks table
results_long <- pivot_longer(unique(res_corrected),
                             cols = c('GeneA', 'GeneB'),
                             names_to = 'A_or_B',
                             values_to = 'Gene')
# Adding hallmarks info and then removing NA
a <- merge(results_long,
           hallmarks,
           by.x = 'Gene',
           by.y = 'Genes')
a_long <- na.omit(pivot_longer(a,
                               cols = colnames(a)[colnames(a) %!in% c('Gene', 'Stage', 'Gene_Pair', 'A_or_B')],
                               values_to = 'Presence',
                               names_to = 'Hallmark'))
# Define the possible hallmarks for each pair
a_wide <- na.omit(pivot_wider(a_long[colnames(a_long) != 'Gene'],
                              names_from = c('A_or_B'),
                              values_from = 'Presence'))
rm(a, a_long, results_long)



### ... Pairs per hallmark ----
hallmarks_pair_count_total <- as.data.frame(table(unique(a_wide[c('Stage', 'Gene_Pair', 'Hallmark')])['Hallmark']))
names(hallmarks_pair_count_total) <- c('Hallmark', 'Total_Count')
hallmarks_pair_count_stages <- as.data.frame(table(unique(a_wide[c('Gene_Pair', 'Hallmark', 'Stage')])[c('Hallmark', 'Stage')])) %>% 
  pivot_wider(values_from = 'Freq',
              names_from = 'Stage')
hallmarks_pair_count <- merge(hallmarks_pair_count_total,
                              hallmarks_pair_count_stages,
                              by = 'Hallmark')
rm(hallmarks_pair_count_total, hallmarks_pair_count_stages)
# Add hallmarks not matched with 0 counts
hallmarks_pair_count$Hallmark <- as.character(hallmarks_pair_count$Hallmark)
for (hm in colnames(hallmarks)[colnames(hallmarks) != 'Genes']){
  if (hm %!in% hallmarks_pair_count$Hallmark){
    hallmarks_pair_count <- rbind(hallmarks_pair_count, c(hm, 0, 0, 0))
  }}



### ... Compute N random permutations ----
N <- 10000 #Number of permutations
n_sig_pairs_prim <- as.numeric(nrow(filter(results, Stage == 'Primary'))) #Number of primary false pairs to sample
n_sig_pairs_meta <- as.numeric(nrow(filter(results, Stage == 'Metastasis'))) #Number of metastasis false pairs to sample

# Function for permutations
perm <- function(times, n_sig_pairs){
  j <- lapply(c(1:times), function(permutation){
    # Select randomly n_sig_pairs false pairs randomly
    fp <- false_pairs[sample(1:nrow(false_pairs), n_sig_pairs, replace = F),]
    # Change to long format
    fp_long <- pivot_longer(fp[c('Gene_Pair', 'GeneA', 'GeneB')],
                            cols = c('GeneA', 'GeneB'),
                            names_to = 'A_or_B',
                            values_to = 'Gene')
    # Merge to hallmarks and change to long format to remove NAs
    b <- merge(fp_long,
               hallmarks,
               by.x = 'Gene',
               by.y = 'Genes')
    b_long <- na.omit(pivot_longer(b,
                                   cols = colnames(b)[colnames(b) %!in% c('Gene', 'Gene_Pair', 'A_or_B')],
                                   values_to = 'Presence',
                                   names_to = 'Hallmark'))
    b_wide <- na.omit(pivot_wider(b_long[colnames(b_long) != 'Gene'],
                                  names_from = c('A_or_B'),
                                  values_from = 'Presence'))
    # Pairs per hallmark
    # length(unique(b_wide$Gene_Pair)) # UNIQUE pairs
    hallmarks_perm_count <- as.data.frame(table(unique(b_wide[c('Gene_Pair', 'Hallmark')])['Hallmark']))
    if (ncol(hallmarks_perm_count) > 1) {
      colnames(hallmarks_perm_count) <- c('Hallmark', paste0('Perm_',permutation,'_Count'))
      return(hallmarks_perm_count)
    } else {
      tempname <- paste0('Perm_',permutation,'_Count')
      return(data.frame('Hallmark' = 'Angiogenesis',
                        tempname = NA))
    }
  })
  return(j)
}



### ... Compute and join the permutation output in a single table ----
# Primary
set.seed(33) # NANOSEED
perm_table_prim <- perm(N, n_sig_pairs_prim) %>% 
  reduce(full_join, by = 'Hallmark')
row.names(perm_table_prim) <- perm_table_prim$Hallmark
perm_table_prim$Hallmark <- NULL
perm_table_prim <- as.data.frame(t(perm_table_prim))
# Add those hallmarks that do not appear in permutations
for (hm in colnames(hallmarks)[colnames(hallmarks) != 'Genes']){
  if (hm %!in% colnames(perm_table_prim)){
    perm_table_prim[,hm] <- NA
  }}

# Metastasis
perm_table_meta <- perm(N, n_sig_pairs_meta) %>% 
  reduce(full_join, by = 'Hallmark')
row.names(perm_table_meta) <- perm_table_meta$Hallmark
perm_table_meta$Hallmark <- NULL
perm_table_meta <- as.data.frame(t(perm_table_meta))
# Add those hallmarks that do not appear in permutations
for (hm in colnames(hallmarks)[colnames(hallmarks) != 'Genes']){
  if (hm %!in% colnames(perm_table_meta)){
    perm_table_meta[,hm] <- NA
  }}



### ... Change NA for 0s ----
perm_table_prim <- replace(perm_table_prim, is.na(perm_table_prim), 0)
perm_table_meta <- replace(perm_table_meta, is.na(perm_table_meta), 0)



### ... Compute empirical pvalues (one side only) ----
pairs_pval <- lapply(c('Primary', 'Metastasis'), function(Stage){
  # Select table
  ifelse(Stage == 'Primary',
         perm_table <- perm_table_prim,
         perm_table <- perm_table_meta)
  # Loop per stage
  p_vals<- lapply(hallmarks_pair_count$Hallmark, function(hm){
    # Loop per hallmark
    distribution <- as.data.frame(table(perm_table[, hm]))
    distribution$Var1 <- as.numeric(as.character(distribution$Var1))
    n_perm_obs <- sum(filter(distribution, Var1 >= as.numeric(filter(hallmarks_pair_count, Hallmark == hm)[,Stage]))$Freq)
    empirical_p_value <- n_perm_obs / sum(distribution$Freq)
    return(empirical_p_value)
    })
  return(p_vals)
})
hallmarks_pair_count$Pvalue_Metastasis <- as.numeric(pairs_pval[[2]])
hallmarks_pair_count$P_value_Primary <- as.numeric(pairs_pval[[1]])



### ... Compute averages and standard error ----
hallmarks_pair_count$Perm_Avg_Metastasis <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  avg <- mean(perm_table_meta[,hm], na.rm = T)
  return(avg)}))
hallmarks_pair_count$Perm_Avg_Primary <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  avg <- mean(perm_table_prim[,hm], na.rm = T)
  return(avg)}))
hallmarks_pair_count$Perm_Error_Metastasis <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  err <- std.error(perm_table_meta[,hm], na.rm = T)
  return(err)}))
hallmarks_pair_count$Perm_Error_Primary <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  err <- std.error(perm_table_prim[,hm], na.rm = T)
  return(err)}))



### ... Frequencies ----
hallmarks_pair_count$Freq_Metastasis <- as.numeric(hallmarks_pair_count$Metastasis) / n_sig_pairs_meta
hallmarks_pair_count$Freq_Primary <- as.numeric(hallmarks_pair_count$Primary) / n_sig_pairs_prim
hallmarks_pair_count$Avg_Perm_Freq_Metastasis <- hallmarks_pair_count$Perm_Avg_Metastasis / n_sig_pairs_meta
hallmarks_pair_count$Avg_Perm_Freq_Primary <- hallmarks_pair_count$Perm_Avg_Primary / n_sig_pairs_prim
hallmarks_pair_count$Error_Perm_Freq_Metastasis <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  err <- std.error(perm_table_prim[,hm] / n_sig_pairs_meta, na.rm = T)
  return(err)}))
hallmarks_pair_count$Error_Perm_Freq_Primary <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  err <- std.error(perm_table_prim[,hm] / n_sig_pairs_prim, na.rm = T)
  return(err)}))



### ... Differences ----
hallmarks_pair_count$Diff_Obs <- hallmarks_pair_count$Freq_Primary - hallmarks_pair_count$Freq_Metastasis
# Permuted differences
perm_freq_table_prim <- perm_table_prim / n_sig_pairs_prim
perm_freq_table_meta <- perm_table_meta / n_sig_pairs_meta
perm_freq_table_diff <- perm_freq_table_prim - perm_freq_table_meta[colnames(perm_freq_table_prim)]
# Permuted differences p_value
hallmarks_pair_count$Avg_Diff_Random <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  avg <- mean(perm_freq_table_diff[,hm], na.rm = T)
  return(avg)}))
# Empirical p-value (tail considered depending the sign of Diff_Obs)
hallmarks_pair_count$Diff_Pval <- unlist(lapply(hallmarks_pair_count$Hallmark, function(hm){
  if (filter(hallmarks_pair_count, Hallmark == hm)$Diff_Obs > 0) {
    p_val <- sum(perm_freq_table_diff[,hm] >= filter(hallmarks_pair_count, Hallmark == hm)$Diff_Obs) / 10000
  } else {
    p_val <- sum(perm_freq_table_diff[,hm] <= filter(hallmarks_pair_count, Hallmark == hm)$Diff_Obs) / 10000
  }
  return(p_val)
}))



### ... Saving file ----
write.table(hallmarks_pair_count,
            sprintf('./DATA/ANALYSIS_DATA/3way__TG/3wTG_%s_%s_%s_FDR%s_hallmarks.tsv', FDR, FREQ, SPLITMOD, thr),
            sep = '\t', row.names = F, quote = F)

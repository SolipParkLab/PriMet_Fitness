
### ... Loading libraries ----
library(dplyr)
library(tidyr)
library(plotrix)
library(purrr)
library(ggplot2)



### ... Variables ----
SPLITMOD <- "Tissue-Stage-PM"
freqmut_threshold <- "0.01"
freqcnv_threshold <- "0.10"
FREQ <- paste0("mf",
               as.numeric(freqmut_threshold)*100,
               "-cf",
               as.numeric(freqcnv_threshold)*100)



### ... Loading functions ----
source("./CODE/common_reg-model-functions.R",local=T)



### ... Load files ----
input_genes <- readRDS('./DATA/PROCESSED_DATA/492_input_genes.RDS')
hallmarks <- read.table('./DATA/PROCESSED_DATA/CELL_CPTAC_Hallmark.txt',
                        sep = '\t', na.strings = '', header = T)
results <- filter(read.delim('./DATA/ANALYSIS_DATA/2way__OG/2wOG_PERM_gene-perturbation_mf1-cf10_Tissue-Stage-PM.tsv',
                             sep = '\t'),
                  Behaviour == 'Perturbed')


### ... Merging hallmarks information ----
a <- merge(unique(results[,c('Gene','Function')]),
           hallmarks,
           by.x = 'Gene',
           by.y = 'Genes')
a_long <- na.omit(pivot_longer(a,
                               cols = colnames(a)[colnames(a) %!in% c('Gene', 'Function')],
                               values_to = 'Presence',
                               names_to = 'Hallmark'))



### ... Pairs per hallmark ----
hallmarks_gene_count <- as.data.frame(table(unique(a_long[c('Gene', 'Hallmark')])['Hallmark']))
names(hallmarks_gene_count) <- c('Hallmark', 'Total_Count')

# Add those hallmarks not matched with 0 counts
hallmarks_gene_count$Hallmark <- as.character(hallmarks_gene_count$Hallmark)
for (hm in colnames(hallmarks)[colnames(hallmarks) != 'Genes']){
  if (hm %!in% hallmarks_gene_count$Hallmark){
    hallmarks_gene_count <- rbind(hallmarks_gene_count, c(hm, 0))
  }}


### ... Compute observed frequencies ----
n_genes <- as.numeric(length(results$Gene))
hallmarks_gene_count$Total_Count <- as.numeric(hallmarks_gene_count$Total_Count)
hallmarks_gene_count$Obs_Freq <- hallmarks_gene_count$Total_Count / n_genes



### ... Compute N random permutations ----
N <- 10000 #Number of permutations
# Function for permutations
perm <- function(times, n_pert_genes){
  j <- lapply(c(1:times), function(permutation){
    # Select randomly n_sig_pairs false pairs randomly
    fg <- input_genes[sample(1:length(input_genes), n_pert_genes, replace = F)]
    # Merge to hallmarks and change to long format to remove NAs
    b <- filter(hallmarks, Genes %in% fg)
    b_long <- na.omit(pivot_longer(b,
                                   cols = colnames(b)[colnames(b) %!in% 'Genes'],
                                   values_to = 'Presence',
                                   names_to = 'Hallmark'))
    # Pairs per hallmark
    hallmarks_perm_count <- as.data.frame(table(unique(b_long[c('Genes', 'Hallmark')])['Hallmark']))
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
set.seed(33) # NANOSEED
perm_table_genes <- perm(N, n_genes) %>% 
  reduce(full_join, by = 'Hallmark')
row.names(perm_table_genes) <- perm_table_genes$Hallmark
perm_table_genes$Hallmark <- NULL
perm_table_genes <- as.data.frame(t(perm_table_genes))
# Add those hallmarks that do not appear in permutations
for (hm in colnames(hallmarks)[colnames(hallmarks) != 'Genes']){
  if (hm %!in% colnames(perm_table_genes)){
    perm_table_genes[,hm] <- NA
  }}
perm_table_genes <- replace(perm_table_genes, is.na(perm_table_genes), 0)



### ... Random frequencies ----
perm_freq_table_genes <- perm_table_genes / n_genes
hallmarks_gene_count$Avg_Random_Obs <- unlist(lapply(hallmarks_gene_count$Hallmark, function(hm){
  avg <- mean(perm_table_genes[,hm], na.rm = T)
  return(avg)}))
hallmarks_gene_count$Avg_Random_Freq <- unlist(lapply(hallmarks_gene_count$Hallmark, function(hm){
  avg <- mean(perm_freq_table_genes[,hm], na.rm = T)
  return(avg)}))



### Compute empirical P_values for frequencies
hallmarks_gene_count$Freq_Pvalue <- unlist(lapply(hallmarks_gene_count$Hallmark, function(hm){
  p_val <- sum(perm_freq_table_genes[,hm] >= filter(hallmarks_gene_count, Hallmark == hm)$Obs_Freq) / 10000
  return(p_val)
}))
# Save results
write.table(hallmarks_gene_count,
            sprintf('./DATA/ANALYSIS_DATA/2way__OG/2wOG_%s_perturbed-genes-hallmarks_%s_%s.tsv',
                    'PERM', FREQ, SPLITMOD),
            sep = '\t', row.names = F, quote = F)

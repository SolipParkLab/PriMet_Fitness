
### ... Import libraries ----
library(metafor)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)



### ... Loading functions ----
source('./CODE/common_reg-model-functions.R')



### ... Read and process input files ----
# Keep significant pairs only
t_res <- read.delim('./DATA/ANALYSIS_DATA/2way__T/2wT_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')
t_res_sig <- filter(t_res, SIG_FDR10 == T)
basic <- t_res_sig[c('Gene', 'Tissue', 'Stage', 'Treatment', 'CNA_type', 'P_value', 'Estimate_plot')]
# Model inputs
model_inputs_df <- read.delim('./DATA/GLM_INPUTS/2way__T/2w__T__inputs_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')



### ... Merge basic sig pairs information with input information ----
tarone_input <- merge(basic[colnames(basic) %!in% c('Treatment', 'P_value', 'Estimate_plot')],
                      model_inputs_df[c('Gene', 'Tissue', 'Stage', 'Treatment',
                                       'NoMutWT', 'MutWT', 'NoMutLoss', 'MutLoss', 'NoMutGain', 'MutGain')],
                      by = c('Gene', 'Tissue', 'Stage'))
tarone_input <- unique(tarone_input[colnames(tarone_input)])
# Vectors for 2 x 2 tables
# Upper left cell --> Coocurrence
tarone_input$ai <- ifelse(tarone_input$CNA_type == 'Gain', tarone_input$MutGain, tarone_input$MutLoss)
# Upper right cell --> Mut only
tarone_input$bi <- tarone_input$MutWT
# Lower left cell --> CNA only
tarone_input$ci <- ifelse(tarone_input$CNA_type == 'Gain', tarone_input$NoMutGain, tarone_input$NoMutLoss)
# Lower right cell --> Both wild types
tarone_input$di <- tarone_input$NoMutWT
# Top row sum
tarone_input$n1i <- tarone_input$ai + tarone_input$bi
# Bottom row sum
tarone_input$n2i <- tarone_input$ci + tarone_input$di


### Odds Ratios ----
tarone_input$Odds_Ratio <- (tarone_input$ai * tarone_input$di) / (tarone_input$bi * tarone_input$ci)



### Tarone's Test ----
# For each significant result, I have to provide two tables: Treat and NonTreat
# So, the length of each input vector will be 2
tarone_results <- apply(tarone_input, 1, function(row){
  treat <- filter(tarone_input, Gene == row['Gene'] & Tissue == row['Tissue'] &
                    Stage == row['Stage'] & CNA_type == row['CNA_type'] & Treatment == 'Treatment')
  non_treat <- filter(tarone_input, Gene == row['Gene'] & Tissue == row['Tissue'] &
                        Stage == row['Stage'] & CNA_type == row['CNA_type'] & Treatment == 'Non_Treatment')
  if (row['Stage'] == 'Treatment'){
    tarone <- rma.mh(ai = c(treat$ai, non_treat$ai),
                     bi = c(treat$bi, non_treat$bi),
                     ci = c(treat$ci, non_treat$ci),
                     di = c(treat$di, non_treat$di),
                     n1i = c(treat$n1i, non_treat$n1i),
                     n2i = c(treat$n2i, non_treat$n2i),
                     measure = 'OR',
                     add = 0.5)
  } else {
    tarone <- rma.mh(ai = c(non_treat$ai, treat$ai),
                     bi = c(non_treat$bi, treat$bi),
                     ci = c(non_treat$ci, treat$ci),
                     di = c(non_treat$di, treat$di),
                     n1i = c(non_treat$n1i, treat$n1i),
                     n2i = c(non_treat$n2i, treat$n2i),
                     measure = 'OR',
                     add = 0.5)
  }
  return(c('Tarone_Stat' = tarone$TA, 'Tarone_Pval' = tarone$TAp))
})
tarone_results <- as.data.frame(t(tarone_results))
tarone_results <- cbind(tarone_input, tarone_results)



### ... Creating the output files ----
output <- pivot_wider(tarone_results[c('Gene','Tissue','Stage','CNA_type','Treatment','Odds_Ratio','Tarone_Stat','Tarone_Pval')],
                      names_from = 'Treatment',
                      values_from = c('Odds_Ratio'))
output <- output %>% 
  rename('Odds_Treat' = 'Treatment',
       'Odds_NonTreat' = 'Non_Treatment')
output[,'-log10(Pval)'] <- -log10(output$Tarone_Pval)
output[,'log2(Treat) - log2(NonTreat)'] <- log2(output$Odds_Treat) - log2(output$Odds_NonTreat)
full_output <- merge(output,
                     pivot_wider(t_res[c('Gene', 'Tissue', 'Stage', 'Treatment', 'CNA_type', 'SIG_FDR10')],
                                 names_from = 'Treatment',
                                 values_from = c('SIG_FDR10'),
                                 names_prefix = 'SIG_FDR10_'),
                     by = c('Gene', 'Tissue', 'Stage', 'CNA_type'))
# Saving
write.table(full_output, './DATA/ANALYSIS_DATA/2way__T/2wT_TvsNT_Tarone-output.tsv',
            sep = '\t', row.names = F, quote = F)



### ... Load permutation output and running Tarone in permutation ----
perm_input <- readRDS('./DATA/ANALYSIS_DATA/2way__T/2wT_permutation-outputs.RDS')
tarone_res <- full_output
tarone_res$Sig_State <- ifelse(tarone_res$SIG_FDR10_Treatment == T & tarone_res$SIG_FDR10_Non_Treatment == T, 'Both',
                               ifelse(tarone_res$SIG_FDR10_Treatment == T & tarone_res$SIG_FDR10_Non_Treatment == F, 'Treatment',
                                      ifelse(tarone_res$SIG_FDR10_Treatment == F & tarone_res$SIG_FDR10_Non_Treatment == T,
                                             'Non_Treatment', NA)))
tarone_res <- tarone_res[!is.na(tarone_res$Sig_State),]
n_pairs <- as.numeric(nrow(tarone_res))



### ... Run Tarone's test for the same 97 pairs in each randomization
false_output <- lapply(c(1:length(perm_input)), function(i){
  false_output <- perm_input[[i]]
  # Sample n_pairs false pairs
  to_tarone <- merge(false_output,
                     tarone_res[c('Gene','Tissue','Stage','CNA_type')],
                     by = c('Gene','Tissue','Stage','CNA_type'))
  pairs <- unique(to_tarone[c('Gene','Tissue', 'Stage', 'CNA_type')])
  for (column in c('MutCNV', 'MutWT', 'NoMutCNV', 'NoMutWT')) {to_tarone[,column] <- as.numeric(to_tarone[,column])}
  # Run Tarone test
  tarone_results <- apply(pairs, 1, function(row){
    treat <- filter(to_tarone, Gene == row['Gene'] & Tissue == row['Tissue'] &
                      Stage == row['Stage'] & CNA_type == row['CNA_type'] & Treatment == 'Treatment')
    non_treat <- filter(to_tarone, Gene == row['Gene'] & Tissue == row['Tissue'] &
                          Stage == row['Stage'] & CNA_type == row['CNA_type'] & Treatment == 'Non_Treatment')
    tarone <- rma.mh(ai = c(treat$MutCNV, non_treat$MutCNV),
                     bi = c(treat$MutWT, non_treat$MutWT),
                     ci = c(treat$NoMutCNV, non_treat$NoMutCNV),
                     di = c(treat$NoMutWT, non_treat$NoMutWT),
                     n1i = c(treat$MutCNV + treat$MutWT, non_treat$MutCNV + non_treat$MutWT),
                     n2i = c(treat$NoMutCNV + treat$NoMutWT, non_treat$NoMutCNV + non_treat$NoMutWT),
                     measure = 'OR',
                     add = 0.5)
    return(c('Tarone_Stat' = tarone$TA, 'Tarone_Pval' = tarone$TAp))
  })
  tarone_results <- as.data.frame(t(tarone_results))
  pairs[,c(paste0('Tarone_Stat','_',i))] <- tarone_results$Tarone_Stat
  pairs[,c(paste0('Tarone_Pval','_',i))] <- tarone_results$Tarone_Pval
  return(pairs)
})
false_output_concat <- false_output %>% 
  reduce(full_join, by = c('Gene', 'Tissue', 'Stage', 'CNA_type'))
saveRDS(false_output_concat, './DATA/ANALYSIS_DATA/2way__T/2wT_TvsNT_Tarone-randomization-output.RDS')



### ... Selecting only Tarone's statistic ----
t_random_stats <- false_output_concat[,c('Gene', 'Tissue', 'Stage', 'CNA_type',
                                         colnames(false_output_concat)[grepl('Stat', colnames(false_output_concat))])]



### ... Compute empirical P-value for each pair ----
tarone_res$Tarone_Empirical_Pval <- apply(tarone_res, 1, function(row){
  filt <- filter(t_random_stats, Gene == row['Gene'] & Tissue == row['Tissue'] & Stage == row['Stage'] & CNA_type == row['CNA_type'])
  comp <- filt[,c(5:ncol(filt))] >= as.numeric(row['Tarone_Stat'])
  return('Tarone_Emp_Pval' = rowSums(comp)/100)
})


### ... Adjust P-value using Bonferroni and ALL pairs ----
tarone_res$Adj_TarEmpPval_AllGenes_BH <- p.adjust(tarone_res$Tarone_Empirical_Pval, method = 'fdr')
tarone_res$SIG_BH10_AllGenes_Emp <- ifelse(tarone_res$Adj_TarEmpPval_AllGenes_BH <= 0.1, T, F)
tarone_res$SIG_BH20_AllGenes_Emp <- ifelse(tarone_res$Adj_TarEmpPval_AllGenes_BH <= 0.2, T, F)
# Saving file
write.table(tarone_res, './DATA/ANALYSIS_DATA/2way__T/2wT_TvsNT_Tarone-output-adj-pval.tsv',
            sep = '\t', row.names = F, quote = F)

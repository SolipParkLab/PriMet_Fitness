
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
# Keep significant pairs only from treatment results
t_res <- read.delim('./DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_PERM_analysis_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')
t_res_sig <- filter(t_res, SIG_FDR10 == T)
basic <- t_res_sig[c('Gene', 'Tissue', 'Stage', 'Treatment', 'CNA_type', 'P_value', 'Estimate_plot')]
# Glm inputs
model_inputs_df <- read.delim('./DATA/GLM_INPUTS/2way_Treatment/2way_Treatment_glm-inputs_mf1-cf10_Tissue-Stage-PM.tsv', sep = '\t')



### ... Merge basic sig pairs information with input information ----
tarone_input <- merge(basic[colnames(basic) %!in% c('Stage', 'P_value', 'Estimate_plot')],
                      model_inputs_df[c('Gene', 'Tissue', 'Stage', 'Treatment',
                                        'NoMutWT', 'MutWT', 'NoMutLoss', 'MutLoss', 'NoMutGain', 'MutGain')],
                      by = c('Gene', 'Tissue', 'Treatment'))
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



### ... Odds Ratios ----
tarone_input$Odds_Ratio <- (tarone_input$ai * tarone_input$di) / (tarone_input$bi * tarone_input$ci)



### ... Tarone's Test ----
# For each significant result, I have to provide two tables: Treat and NonTreat
# So, the length of each input vector will be 2
tarone_results <- apply(tarone_input, 1, function(row){
  pri <- filter(tarone_input, Gene == row['Gene'] & Tissue == row['Tissue'] &
                  Treatment == row['Treatment'] & CNA_type == row['CNA_type'] & Stage == 'Primary')
  meta <- filter(tarone_input, Gene == row['Gene'] & Tissue == row['Tissue'] &
                   Treatment == row['Treatment'] & CNA_type == row['CNA_type'] & Stage == 'Metastasis')
  if (row['Stage'] == 'Primary'){
    tarone <- rma.mh(ai = c(pri$ai, meta$ai),
                     bi = c(pri$bi, meta$bi),
                     ci = c(pri$ci, meta$ci),
                     di = c(pri$di, meta$di),
                     n1i = c(pri$n1i, meta$n1i),
                     n2i = c(pri$n2i, meta$n2i),
                     measure = 'OR',
                     add = 0.5)
  } else {
    tarone <- rma.mh(ai = c(meta$ai, pri$ai),
                     bi = c(meta$bi, pri$bi),
                     ci = c(meta$ci, pri$ci),
                     di = c(meta$di, pri$di),
                     n1i = c(meta$n1i, pri$n1i),
                     n2i = c(meta$n2i, pri$n2i),
                     measure = 'OR',
                     add = 0.5)
  }
  return(c('Tarone_Stat' = tarone$TA, 'Tarone_Pval' = tarone$TAp))
})
tarone_results <- as.data.frame(t(tarone_results))
tarone_results <- cbind(tarone_input, tarone_results)



### ... Creating the output files ----
output <- pivot_wider(tarone_results[c('Gene','Tissue','Stage','CNA_type','Treatment','Odds_Ratio','Tarone_Stat','Tarone_Pval')],
                      names_from = 'Stage',
                      values_from = c('Odds_Ratio'))
output <- output %>% 
  rename('Odds_Pri' = 'Primary',
         'Odds_Meta' = 'Metastasis')
output[,'-log10(Pval)'] <- -log10(output$Tarone_Pval)
output[,'log2(Treat) - log2(NonTreat)'] <- log2(output$Odds_Pri) - log2(output$Odds_Meta)
full_output <- merge(output,
                     pivot_wider(t_res[c('Gene', 'Tissue', 'Stage', 'Treatment', 'CNA_type', 'SIG_FDR10')],
                                 names_from = 'Stage',
                                 values_from = 'SIG_FDR10',
                                 names_prefix = 'SIG_FDR10_'),
                     by = c('Gene', 'Tissue', 'Treatment', 'CNA_type'))
# Saving
write.table(full_output, './DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_PvsM_Tarone-output.tsv',
            sep = '\t', row.names = F, quote = F)



### ... Load permutation output ----
perm_input <- readRDS('./DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_permutation-outputs.RDS')
tarone_res <- full_output
tarone_res$Sig_State <- ifelse(tarone_res$SIG_FDR10_Primary == T & tarone_res$SIG_FDR10_Metastasis == T, 'Both',
                               ifelse(tarone_res$SIG_FDR10_Primary == T & tarone_res$SIG_FDR10_Metastasis == F, 'Primary',
                                      ifelse(tarone_res$SIG_FDR10_Primary == F & tarone_res$SIG_FDR10_Metastasis == T,
                                             'Metastasis', NA)))
tarone_res <- tarone_res[!is.na(tarone_res$Sig_State),]
n_pairs <- as.numeric(nrow(tarone_res))



### ... Run Tarone's test for the same 97 pairs in each randomization ----
false_output <- lapply(c(1:length(perm_input)), function(i){
  false_output <- perm_input[[i]]
  # Sample n_pairs false pairs
  to_tarone <- merge(false_output,
                     tarone_res[c('Gene','Tissue','Treatment','CNA_type')],
                     by = c('Gene','Tissue','Treatment','CNA_type'))
  pairs <- unique(to_tarone[c('Gene','Tissue', 'Treatment', 'CNA_type')])
  for (column in c('MutCNV', 'MutWT', 'NoMutCNV', 'NoMutWT')) {to_tarone[,column] <- as.numeric(to_tarone[,column])}
  # Run Tarone test
  tarone_results <- apply(pairs, 1, function(row){
    pri <- filter(to_tarone, Gene == row['Gene'] & Tissue == row['Tissue'] &
                    Treatment == row['Treatment'] & CNA_type == row['CNA_type'] & Stage == 'Primary')
    meta <- filter(to_tarone, Gene == row['Gene'] & Tissue == row['Tissue'] &
                     Treatment == row['Treatment'] & CNA_type == row['CNA_type'] & Stage == 'Metastasis')
    tarone <- rma.mh(ai = c(pri$MutCNV, meta$MutCNV),
                     bi = c(pri$MutWT, meta$MutWT),
                     ci = c(pri$NoMutCNV, meta$NoMutCNV),
                     di = c(pri$NoMutWT, meta$NoMutWT),
                     n1i = c(pri$MutCNV + pri$MutWT, meta$MutCNV + meta$MutWT),
                     n2i = c(pri$NoMutCNV + pri$NoMutWT, meta$NoMutCNV + meta$NoMutWT),
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
  reduce(full_join, by = c('Gene', 'Tissue', 'Treatment', 'CNA_type'))
saveRDS(false_output_concat, './DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_PvsM_Tarone-randomization-output.RDS')


### ... Load permutation output ----
t_random_stats <- false_output_concat[,c('Gene', 'Tissue', 'Treatment', 'CNA_type',
                                         colnames(false_output_concat)[grepl('Stat', colnames(false_output_concat))])]



### ... Compute empirical P-value for each pair ----
tarone_res$Tarone_Empirical_Pval <- apply(tarone_res, 1, function(row){
  filt <- filter(t_random_stats, Gene == row['Gene'] & Tissue == row['Tissue'] & Treatment == row['Treatment'] & CNA_type == row['CNA_type'])
  comp <- filt[,c(5:ncol(filt))] >= as.numeric(row['Tarone_Stat'])
  return('Tarone_Emp_Pval' = rowSums(comp)/100)
})



### ... Adjust P-value using Bonferroni and ALL pairs ----
tarone_res$Adj_TarEmpPval_AllGenes_BH <- p.adjust(tarone_res$Tarone_Empirical_Pval, method = 'fdr')
tarone_res$SIG_BH10_AllGenes_Emp <- ifelse(tarone_res$Adj_TarEmpPval_AllGenes_BH <= 0.1, T, F)
tarone_res$SIG_BH20_AllGenes_Emp <- ifelse(tarone_res$Adj_TarEmpPval_AllGenes_BH <= 0.2, T, F)
# Saving file
write.table(tarone_res, './DATA/ANALYSIS_DATA/2way_Treatment/2way_Treatment_PvsM_Tarone-output-adj-pval.tsv',
            sep = '\t', row.names = F, quote = F)

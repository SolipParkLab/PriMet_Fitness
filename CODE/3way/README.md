## 3 way code

The code here allows to run the 3-way log-linear regression models. 

1.  glm-inputs: code to compute the inputs for the log-linear regression models from the binary matrix. Eight classes in total are required depending on whether geneB has mutations (NoMutA_WT_MutB, MutA_WT_MutB, NoMutA_CNA_MutB and MutA_CNA_MutB) or not (NoMutA_WT_NoMutB, MutA_WT_NoMutB, NoMutA_CNA_NoMutB and MutA_CNA_NoMutB).

2. glm-outputs: code to run the log-linear regression models.

3. Permut: code to obtain the number of false significant pairs and compute the average FDR for 3-way models.

4. fdr-correction: false discovery rate correction of the log-linear regression's *P*-value based on the permuted counts.

5. Hallmarks: code to perform the cancer hallmarks analysis for the significant 3-way pairs.
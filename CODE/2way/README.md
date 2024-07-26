## 2 way code

This code allows to run the 2-way models, including oncogenic (across cancer types and subtypes), treatment and confounding factors models.

- glm-inputs: code to compute the inputs for the log-linear regression models from the binary matrices. Four types of samples are required: (i) NoMut_WT, the number of samples that have neither mutation nor CNAs; (ii) Mut_WT, the number of samples with only mutations and no CNAs; (iii) NoMut_CNA, the number of samples with only CNAs and no mutations, and (iv) Mut_CNA, the number of samples with both mutation and CNAs.

- glm-outputs: code to run the log-linear regression models.

- Permut: code to compute the permutated matrices, run the log-linear regression models using the permuted matrices, obtain the number of false significant pairs and compute the average FDR as described in the Methods section.

- fdr-class (if applied): false discovery rate correction of the log-linear regression's *P*-value based on the permuted counts and classification of tested genes as described by Park *et al* in 2021.

- Additional analysis: code for the specific analysis we have done, such as perturbation status in the oncogenic models or odds ratio heterogeneity test in treatment models.

**Side note**: in the case of position specific analysis in treatment, all the steps are included within the same script.
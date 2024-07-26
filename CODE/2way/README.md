## 2 way code

This code allows to run the 2-way models, including oncogenic (across cancer types and subtypes), treatment and confounding factors models.

- glm-inputs: code to compute the inputs for the log-linear regression models from the binary matrices.

- glm-outputs: code to run the log-linear regression models.

- Permut: code to run the permutation matrices, run the log-linear regression models using the permuted matrices, obtain the number of false significant pairs and compute the average FDR as described in the Methods section.

- fdr-class (if applied): false discovery rate correction of the log-linear regression *P*-value based on the permuted counts and classification of tested genes as described in Park *et al* (2019).

- Additional analysis: code for the additional analysis we have done, such as perturbation status in the oncogenic models or odds ratio heterogeneity in treatment.

Side note: in the case of position specific analysis in treatment, all the steps are included within the same script.
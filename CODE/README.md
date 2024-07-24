### In this folder you can find all the code necesary to run all the analysis in our work.

- *2way* folder has the code for all the 2-way log-linear regression models, including oncogenic (across cancer types and subtypes), treatment, aggressiveness (number of metastasis) and confounding factors models.

- *3way* folder has the code for the 3-way log-linear regression models.

- *Permutations* folder has the code to generate the permuted matrices for 2-way and 3-way models and compute the False Discovery Rate (FDR).

- *Preprocessing* folder has the code to process the raw data (downloaded as specified in the Methods section)

- *Survival* folder has the code to classify the samples according their mutation status based on the significant and non-sifnificant 2-way and 3-way interactions and run the Cox regression models.

The additional script includes functions used by other scripts.
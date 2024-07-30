# PriMet_Fitness

**CODE** includes the methods for analyzing co-occurrences between mutations and CNAs within a cancer gene using log-linear regression models. We used log-linear regression models with Poisson functions to determine the significance of these co-occurrences.

In **DATA**, we provide the input files for the log-linear regression models.

For the 2-way models, four types of sample counts are required:
1. NoMut_WT: Samples with neither mutation nor CNAs.
2. Mut_WT: Samples with only mutations and no CNAs.
3. NoMut_CNA: Samples with only CNAs and no mutations.
4. Mut_CNA: Samples with both mutations and CNAs. Copy number gain and loss events have been analyzed independently.

For 3-way interactions, these four classes need to be counted twice, depending on the mutation status of the second gene (geneB).
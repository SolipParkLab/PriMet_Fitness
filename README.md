# PriMet_Fitness

**CODE** includes the methods for analysing coocurrences between mutations and CNAs within a cancer gene using log-linear regression models. We used log-linear regression models with Poisson functions to determine the significance of the coocurrences within a gene. Several models have been developed:
- 2-way interaction between mutation and CNAs within a gene - tissue - stage pair.
- 2-way interaction between mutation and CNAs within a gene - tissue - stage - treatment_status pair.
- 3-way interaction between two genes with three possible genomic alterations (geneA-mutation:geneA-CNA:geneB-mutation) 

In **DATA** we provide the input files for the log-linear regression models. Four types of samples are required in the case of 2-way models: (i) NoMut_WT, the number of samples that have neither mutation nor CNAs; (ii) Mut_WT, the number of samples with only mutations and no CNAs; (iii) NoMut_CNA, the number of samples with only CNAs and no mutations, and (iv) Mut_CNA, the number of samples with both mutation and CNAs. Copy number gain and loss events have been analyzed independently.

In the case of 3-way interactions, these four classes have to be counted twice, depending on the mutation status of the second gene of the 3-way interaction (i.e geneB). Eight classes in total are required:
- If geneB has no mutations: NoMutA_WT_NoMutB, MutA_WT_NoMutB, NoMutA_CNA_NoMutB and MutA_CNA_NoMutB.
- If geneB has mutations: NoMutA_WT_MutB, MutA_WT_MutB, NoMutA_CNA_MutB and MutA_CNA_MutB.
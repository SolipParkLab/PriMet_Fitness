## DATA

Files that are necessary to run the scripts are described here. Some of them were dowloaded directly from the source, so we describe where to download them.

* Clinical information, somatic mutation data and copy-number alteration data can be downloaded from [cBioPortal](https://www.cbioportal.org/study/summary?id=msk_met_2021)
* Treatment data can be collected from Supplementary table 1B and metastatic biopsy location data from Supplementary table 1C in the [MSK-MET publication](https://www.sciencedirect.com/science/article/pii/S0092867422000034?via%3Dihub) (DOI: 10.1016/j.cell.2022.01.003)
* Maf with clonality information is created by merging the maf file annotated with [facets-suite](https://github.com/mskcc/facets-suite) and the clinical data (check preprocessing the script for exact procedure).
* Hallmarks file can be downloaded from [CPTAC study](https://www.cell.com/cell/fulltext/S0092-8674%2823%2900780-8?dgcid=raven_jbs_aip_email) (DOI: 10.1016/j.cell.2023.07.014)
* The binary matrix (**Binary_Matrix.tsv**) was created in the preprocessing script and merged to clinical data. Confounding factor levels (HIGH or LOW) were assigned using the median.
* **492_input_genes.tsv** is the list of genes in the maf file.
* **p_cancer-gene-list.tsv** includes a list of cancer genes and their function as tumor drivers (oncogene, tumor supressor or both). Can be downloaded from [oncoKB](https://www.oncokb.org/cancer-genes).
* **p_gene-names_conversion_table.tsv** includes a list of genes from the maf file that don't have the current HUGO symbol. This table also provides the current HUGO symbols.
* **clonality_results.tsv** has sample clonality information, obtained from the functional clonal maf (created in the preprocessing script).
* **Surv-clinical-data.tsv** has information used to create the sample classification based on 2-way and 3-way pairs for Cox regression models.

Additionaly, we provide sample counts for 2-way and 3-way log-linear regression models.
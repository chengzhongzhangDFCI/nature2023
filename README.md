# Project Repository: 
## Heritable transcriptional defects from aberrations of nuclear architecture

### Authors: 
- Nikos Mynhier
- Gregory Brunette
- Lanting Li
- Etai Jacob
- Cheng-Zhong Zhang

### README:

This repository holds the code required to rerun the analysis performed on the article's sequencing data. 

### Preprocessing 

This directory holds the source (src) code for our snakemake pipeline that is used to align and preprocess sequencing data. 

This code serves as a reference for we processed our sequencing data. To run this code one must adjust the snakemake configuration appropriately and download raw sequencing data from the public repository. 

### scRNA (single cell RNA sequencing analysis)

This directory holds the source code and data needed to perform the single cell sequencing analysis articulated in the methods of the article. 

This code can be used to rerun our analysis from provided summary data files. Our summary data files include aggregate TPM and allele specific expression data from each cell in the study.

Summary data files can  be found at the following location:
/nature_2023_submission_code_2/scRNA/data/analysis_data

To rerun our analysis: 
1) Adjust the paths in /nature_2023_submission_code_2/scRNA/src/analysis_scripts/analysis_macro.R 
2) Download necessary packages based on those required in analysis_macro.R or using /nature_2023_submission_code_2/preprocessing/src/environment/R-packages-to-download.R
3) Run /nature_2023_submission_code_2/scRNA/src/analysis_scripts/analysis_macro.R 

#Note: This will create many data file intermediates and plots.

### ATACSeq

### BulkRNA

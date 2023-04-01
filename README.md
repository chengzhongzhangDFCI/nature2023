# Project Repository: 
## Heritable transcriptional defects from aberrations of nuclear architecture

### Authors: 
- Nikos Mynhier (single-cell RNA-Seq analysis)
- Lanting Li (bulk RNA-Seq analysis)
- Gregory Brunette (bulk ATAC-Seq analysis and DNA copy-number analysis)
- Shiwei Liu and Stamatis Papathanasiou (imaging data analysis)

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

For any questions about the scripts or generated plots please feel free to reach out to us for clarification. 

### bulkATAC

Source code and results for bulk ATAC sequencing of RPE-1 chr4 bridge clones.

./bulkATAC/data/BigWigs --> normalized ATAC peak signal after normalization for library size and copy-number
./bulkATAC/data/FC_tables --> genome-wide normalized ATAC signal for 1 Mb intervals
./bulkATAC/data/ATAC_fragcounts_raw.txt --> unprocessed ATAC seq fragment counts in peaks
./bulkATAC/data/bg_gc.dat --> ATAC seq background peak set for permutation analysis

### bulkRNA

See `./bulkRNA/src/workbook_bulkRNAanalysis.ipynb` for details.

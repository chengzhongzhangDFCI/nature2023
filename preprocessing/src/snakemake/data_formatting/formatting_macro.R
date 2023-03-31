##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

########################################
### Input Desired script Directories ###
########################################
#args <- commandArgs(trailingOnly = TRUE)
args <- c("/Users/nikos_mynhier/Documents/GitHub/nature_2023_submission_code_2/preprocessing/src/snakemake", "/Users/nikos_mynhier/Documents/GitHub/nature_2023_submission_code_2/scRNA/data",  "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
experiments <- args[3:length(args)]

#######################
### Source Packages ###
#######################

require(tidyverse)  # data manipulation
require(cluster)    # clustering algorithms
require(factoextra) # clustering algorithms & visualization
require(reshape)
require(data.table)
require(gtools)
require(readxl)
require(ggplot2)
require(matrixStats)
require(reshape2)
require(gridExtra)
require(dplyr)
require(ggthemes)
require(rlang)
require(stringr)
require(readr)

###################
### Run scripts ###
###################

#Collect and process ASE data
#WARNING: LONG RUN TIME
message("Running QC")
source( sprintf('%s/data_formatting/QC_calculations.R', scriptsdir) )

#Collect and process ASE data
#WARNING: LONG RUN TIME
message("Aggregating ASE")
source( sprintf('%s/data_formatting/ASE_calculations_v2.R', scriptsdir) )

#Collect and process mRNA transcript data
#WARNING: LONG RUN TIME
message("Aggregating TPM")
source( sprintf('%s/data_formatting/TPM_calculations.R', scriptsdir) )














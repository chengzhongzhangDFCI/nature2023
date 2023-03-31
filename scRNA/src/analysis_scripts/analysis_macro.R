##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
# This script runs, in order, all the scripts used in the analysis.
# This script was written by Nikos Mynhier 

########################################
### Input Desired script Directories ###
########################################

#args <- commandArgs(trailingOnly = TRUE)
args <- c("/Users/nikos_mynhier/Documents/GitHub/nature_2023_submission_code_2/scRNA/src", "/Users/nikos_mynhier/Documents/GitHub/nature_2023_submission_code_2/scRNA/data/analysis_data")
print(args)
scriptsdir <- args[1]
dirpath <- args[2]

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

#Calculation TPM ratios and adjust for cell specific effects
message("Normalizing data")
source( sprintf('%s/analysis_scripts/normalize_data.R', scriptsdir) )

#Aggregate by inverse variance at bin, arm, and chr levels
message("Binning data")
source( sprintf('%s/analysis_scripts/agg_by_inv_var.R', scriptsdir) )

#Curate reference distributions for each chr and arm 
message("Creating reference distributions")
source( sprintf('%s/analysis_scripts/create_ref_distr.R', scriptsdir) )

#Calculate p-values at chr and arm level
message("Calculating p-values")
source( sprintf('%s/analysis_scripts/pval_calculations_2.R', scriptsdir) )

#Automatically call abnormal CN families and global summary 
message("Analyzing results")
source( sprintf('%s/analysis_scripts/automated_analysis_4.R', scriptsdir) )

#Generate Visuals
message("Creating individual visuals")
source( sprintf('%s/visuals_scripts/visual_generator_v2.R', scriptsdir) )

#Generate Summary Visuals
message("Creating summary visuals")
source( sprintf('%s/visuals_scripts/summary_visual_generator_2.R', scriptsdir) )














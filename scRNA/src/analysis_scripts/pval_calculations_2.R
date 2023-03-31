##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

###############
### Imports ###
###############

#Annotation list with information on every cell and IDs of control samples
anno <- data.table(read.csv( sprintf("%s/annotation_list.csv", dirpath)))
controlSampleIDs <- readRDS( sprintf("%s/controlSampleIDs.rds", dirpath))

#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#Get MN cells
MN_Cells <- anno_gen_1_2$WTA.plate[which(anno_gen_1_2$Sister1 == 1)]

#import samples that were manually determined to be aneuploid
hand_samples <- readRDS(sprintf("%s/grouped_control_aneuploidies.rds", dirpath))

#import gene and arm annotations
geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos_2.rds", dirpath))
armRanges <- readRDS(sprintf("%s/armRanges.rds", dirpath))

#Read in raw TPM matrix
TPM <- readRDS(file=sprintf("%s/TPM.bygene.rds", dirpath))

#read in TPM by chr
TPM_bychr <- readRDS(file=sprintf("%s/TPM.inv_var.bychr.rds", dirpath))

#read in AS_TPM by chr
AS_TPM_bychr <- readRDS(file=sprintf("%s/AS-TPM.inv_var.bychr.rds", dirpath))

#read in TPM by arm
TPM_byarm <- readRDS(file=sprintf("%s/TPM.inv_var.byarm.rds", dirpath))

#read in AS_TPM by chr
AS_TPM_byarm <- readRDS(file=sprintf("%s/AS-TPM.inv_var.byarm.rds", dirpath))

#################
### Functions ###
#################

#function to get pvals through z-test
pval_ztest <- function(x, pop_distr){
  #x is a vector of points to check against population statistics
  #pop_distr is a reference vector from which to pull statistics
  mu <- mean(pop_distr)
  sd <- sd(pop_distr)
  #calc z-score for all values in x against stats from pop_distr
  z_score <- ((x - mu) / sd)
  #calculate a p-value for each element in the vector of z-scores based on normal distr
  pvals <- 2*pnorm( -abs(z_score) ) #*2 for two tailed test
  #return a p-value vector
  return(pvals)
}

#function to get pvals for a matrix 
run_ztest_for_matrix <- function(mat, pop_distr, anno_cols = 1){
  #run pval_ztest function on non-annotation columns of input matrix
  changeCols <- colnames(mat)[-c(1:anno_cols)]
  pval_mat <- copy(mat)
  pval_mat[,(changeCols):= lapply(.SD, pval_ztest, pop_distr = pop_distr), .SDcols = changeCols] 
  #return pval matrix 
  return(pval_mat)
}

####################################################
### Use ref distributions to calculate TPM pvals ###
####################################################

#get global norm factor as QC
TPM_global_factors <- global_factors(TPM)
saveRDS(TPM_global_factors, file = sprintf("%s/TPM_global_norm_factor.rds", dirpath))

#Pull norm factors for ref samples
mono_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% mono_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% mono_TPM$ID)]))
ctrl_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% ctrl_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% ctrl_TPM$ID)]))
tri_factors <- data.table(cbind(ID = names(TPM_global_factors[which(names(TPM_global_factors) %in% tri_TPM$ID)]), factor = TPM_global_factors[which(names(TPM_global_factors) %in% tri_TPM$ID)]))

#add norm factors to ref data tables
setkey(mono_TPM, ID); setkey(mono_factors, ID); mono_TPM <- merge(mono_TPM, mono_factors);
setkey(ctrl_TPM, ID); setkey(ctrl_factors, ID); ctrl_TPM <- merge(ctrl_TPM, ctrl_factors);
setkey(tri_TPM, ID); setkey(tri_factors, ID); tri_TPM <- merge(tri_TPM, tri_factors);

#run pval calculations on abnormal CN states tables (Not chr specific)
TPM_bychr_mono_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = mono_TPM$vals, anno_cols = 1)
TPM_bychr_tri_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = tri_TPM$vals, anno_cols = 1)
TPM_bychr_ctrl_pvals <- run_ztest_for_matrix(mat=TPM_bychr, pop_distr = ctrl_TPM_noX_no10$vals, anno_cols = 1)

#save pvals for each state
saveRDS(TPM_bychr_mono_pvals, file = sprintf("%s/pval_matrix_loss_bychr.rds", dirpath))
saveRDS(TPM_bychr_ctrl_pvals, file = sprintf("%s/pval_matrix_control_bychr.rds", dirpath))
saveRDS(TPM_bychr_tri_pvals, file = sprintf("%s/pval_matrix_gain_bychr.rds", dirpath))

###################################################################
### Run pval calculations on ref AS_TPM chr level distributions ###
###################################################################

#read in hap spec CN ref distributions
ref_hap_spec_CN <- readRDS(file = sprintf("%s/ref_hap_spec_CN.rds", dirpath))

AS_TPM_bychr_ctrl_pvals_chr_spec <- list()
AS_TPM_bychr_ctrl_pvals_chr_spec$A <- data.table()
AS_TPM_bychr_ctrl_pvals_chr_spec$B <- data.table()
AS_TPM_bychr_mono_pvals_chr_spec <- list()
AS_TPM_bychr_mono_pvals_chr_spec$A <- data.table()
AS_TPM_bychr_mono_pvals_chr_spec$B <- data.table()
for(i in 1:length(unique(ctrl_AS_TPM$chr))){
  chr <- unique(ctrl_AS_TPM$chr)[i]
  chr_str <- paste0("chr", unique(ctrl_AS_TPM$chr)[i])
  if(!chr_str %in% c("chr10b", "chr23")){
    #calculate disomy pvals for allele A
    tmp_pvals_A <- run_ztest_for_matrix(mat=AS_TPM_bychr$A[which(AS_TPM_bychr$A$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "A&B"))], anno_cols = 1)
    AS_TPM_bychr_ctrl_pvals_chr_spec$A <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$A, tmp_pvals_A)
    #calculate monosomy pvals for allele A
    tmp_pvals_A_mono <- run_ztest_for_matrix(mat=AS_TPM_bychr$A[which(AS_TPM_bychr$A$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "A"))], anno_cols = 1)
    AS_TPM_bychr_mono_pvals_chr_spec$A <- rbind(AS_TPM_bychr_mono_pvals_chr_spec$A, tmp_pvals_A_mono)
    #calculate disomy pvals for allele B
    tmp_pvals_B <- run_ztest_for_matrix(mat=AS_TPM_bychr$B[which(AS_TPM_bychr$B$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "A&B"))], anno_cols = 1)
    AS_TPM_bychr_ctrl_pvals_chr_spec$B <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$B, tmp_pvals_B)
    #calculate monosomy pvals for allele B
    tmp_pvals_B_mono <- run_ztest_for_matrix(mat=AS_TPM_bychr$B[which(AS_TPM_bychr$B$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "B"))], anno_cols = 1)
    AS_TPM_bychr_mono_pvals_chr_spec$B <- rbind(AS_TPM_bychr_mono_pvals_chr_spec$B, tmp_pvals_B_mono)
  } else { #only run monosomy pvals because we dont have a disomy reference
    if(chr_str == "chr23") {chr_str <- "chrX"}
    #calculate monosomy pvals for allele A
    tmp_pvals_A_mono <- run_ztest_for_matrix(mat=AS_TPM_bychr$A[which(AS_TPM_bychr$A$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "A"))], anno_cols = 1)
    AS_TPM_bychr_mono_pvals_chr_spec$A <- rbind(AS_TPM_bychr_mono_pvals_chr_spec$A, tmp_pvals_A_mono)
    #calculate monosomy pvals for allele B
    tmp_pvals_B_mono <- run_ztest_for_matrix(mat=AS_TPM_bychr$B[which(AS_TPM_bychr$B$chr == chr_str),], pop_distr = ref_hap_spec_CN$vals[intersect(which(ref_hap_spec_CN$chr == chr), which(ref_hap_spec_CN$hap == "B"))], anno_cols = 1)
    AS_TPM_bychr_mono_pvals_chr_spec$B <- rbind(AS_TPM_bychr_mono_pvals_chr_spec$B, tmp_pvals_B_mono)
  }
}

#save pvals for each state
saveRDS(AS_TPM_bychr_mono_pvals_chr_spec, file = sprintf("%s/pval_AS_TPM_bychr_mono_pvals_chr_spec.rds", dirpath))
saveRDS(AS_TPM_bychr_ctrl_pvals_chr_spec, file = sprintf("%s/pval_AS_TPM_bychr_ctrl_pvals_chr_spec.rds", dirpath))

###################################################################
### Run pval calculations on ref AS_TPM arm level distributions ###
###################################################################

AS_TPM_byarm_ctrl_pvals_arm_spec <- list()
AS_TPM_byarm_ctrl_pvals_arm_spec$A <- data.table()
AS_TPM_byarm_ctrl_pvals_arm_spec$B <- data.table()
for(i in 1:length(unique(ctrl_AS_TPM_byarm_A$chr))){
  chr <- unique(ctrl_AS_TPM_byarm_A$chr)[i]
  row <- which(AS_TPM_byarm$A$arm == chr)
  tmp_pvals_A <- run_ztest_for_matrix(mat=AS_TPM_byarm$A[row,], pop_distr = ctrl_AS_TPM_byarm_A$vals[ctrl_AS_TPM_byarm_A$chr == chr], anno_cols = 1)
  AS_TPM_byarm_ctrl_pvals_arm_spec$A <- rbind(AS_TPM_byarm_ctrl_pvals_arm_spec$A, tmp_pvals_A)
  tmp_pvals_B <- run_ztest_for_matrix(mat=AS_TPM_byarm$B[row,], pop_distr = ctrl_AS_TPM_byarm_B$vals[ctrl_AS_TPM_byarm_B$chr == chr], anno_cols = 1)
  AS_TPM_byarm_ctrl_pvals_arm_spec$B <- rbind(AS_TPM_byarm_ctrl_pvals_arm_spec$B, tmp_pvals_B)
}

saveRDS(AS_TPM_byarm_ctrl_pvals_arm_spec, file = sprintf("%s/AS_TPM_pval_matrix_control_byarm_armspec.rds", dirpath))


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

get_vals_from_mat <- function(samples_df, mat, add_chr = T) {
  vals <- c()
  for(i in 1:nrow(samples_df)) {
    cell <- samples_df$ID[i]
    if(add_chr) {
      chr <- paste0("chr", samples_df$chr[i]); if(length(grep(chr, pattern = 23))>0) {chr <- str_replace(chr, pattern = "23", replacement = "X")};
    } else {chr <- samples_df$chr[i]; if(length(grep(chr, pattern = 23))>0) {chr <- str_replace(chr, pattern = "23", replacement = "X")};}
    row <- which(mat$chr == chr)
    unname(unlist(mat[row, ..cell]))
    vals <- append(vals, unlist(mat[row, ..cell]))
  }
  samples_vals <- cbind(samples_df, vals)
  return(samples_vals)
}

#Get the factors used to do global normalization 
global_factors <- function(mat, anno_cols = 6) {
  #Copy dataframe and remove anno columns
  cols <- c(1:anno_cols)
  dt2go <- copy(data.table(mat)[,-..cols])
  log2go <- log2(dt2go)
  log2go[log2go == -Inf] <- NA
  
  #<(log(TPM_i)>_genes - <<log(TPM_i)>_ctrl>_genes
  cell_means <- colMeans(log2go, na.rm = T)
  ctrl_means <- mean(colMeans(log2go[,..controlSampleIDs], na.rm = T), na.rm = T)
  cell_TPM_diff <- 2^(cell_means - ctrl_means)
  
  #return normed ratios
  return(cell_TPM_diff)
}

##################################################
### Create ref distributions for each CN state ###
##################################################

#split ref samples by CN state
mono_samples <- hand_samples[which(hand_samples$CN == 1),]
ctrl_samples <- data.table(cbind(ID = rep(controlSampleIDs, each=24), CN = 2, chr = rep(c(1:9, "10a", "10b", 11:23), length(controlSampleIDs))))
tri_samples <- hand_samples[which(hand_samples$CN == 3),]

#run get_vals_from_mat function on for each ref group
mono_TPM <- get_vals_from_mat(samples_df=mono_samples, mat=TPM_bychr)
ctrl_TPM <- get_vals_from_mat(samples_df=ctrl_samples, mat=TPM_bychr)
ctrl_TPM_noX_no10 <- ctrl_TPM[-which(ctrl_TPM$chr %in% c("23", "10b")),]
tri_TPM <- get_vals_from_mat(samples_df=tri_samples, mat=TPM_bychr)

#combine ref TPM values
ref_TPM <- rbind(mono_TPM, ctrl_TPM_noX_no10, tri_TPM)

#save ref TPM values
saveRDS(ref_TPM, file = sprintf("%s/ref_TPM.rds", dirpath))

#####################################################
### Create ref distributions for each AS-CN state ###
#####################################################

#run get_vals_from_mat function to get AS_TPM (Hap A) values for each ref cell
ctrl_AS_TPM_A <- get_vals_from_mat(samples_df=ctrl_samples, mat=AS_TPM_bychr$A)

#run get_vals_from_mat function to get AS_TPM (Hap B) values for each ref cell
ctrl_AS_TPM_B <- get_vals_from_mat(samples_df=ctrl_samples, mat=AS_TPM_bychr$B)

#Get AS_TPM from correct haplotypes
ctrl_AS_TPM_AB <- cbind(ctrl_AS_TPM_A[,1:3], A = as.numeric(unlist(ctrl_AS_TPM_A[,4])), B = as.numeric(unlist(ctrl_AS_TPM_B[,4])))
ctrl_AS_TPM <- rbind(cbind(ctrl_AS_TPM_AB[,1:3], hap = "A", vals = unlist(ctrl_AS_TPM_AB[,4])), cbind(ctrl_AS_TPM_AB[,1:3], hap = "B", vals = unlist(ctrl_AS_TPM_AB[,5])))
ctrl_AS_TPM$CN <- as.numeric(ctrl_AS_TPM$CN) - 1
ctrl_AS_TPM_noX_no10 <- ctrl_AS_TPM[-which(ctrl_AS_TPM$chr %in% c("23", "10b")),] #remove these as they are not diploid CN state

#Adjust ctrl_TPM for binding with AS value
ctrl_TPM_noX_no10_2 <- copy(ctrl_TPM_noX_no10) #copy the control (allele NON-specific) values from other matrix
ctrl_TPM_noX_no10_2$hap <- "A&B" #add hap annotation for diploid control cells
ctrl_TPM_noX_no10_2 <- ctrl_TPM_noX_no10_2[,c("ID", "CN", "chr", "hap", "vals")] #keep column order the same as AS_TPM
ctrl_TPM_noX_no10_2$vals <- 2*ctrl_TPM_noX_no10_2$vals #to make the base CN = 2

#bind all together
ref_hap_spec_CN <- rbind(ctrl_AS_TPM, ctrl_TPM_noX_no10_2)

#manually identified trisomies in control cells
ctrl_tri12_cells <- c("181012_5E", "181012_5F", "181012_5H", "170208_A1", "170208_A2", "210720_7A", "170208_D1", "170208_D2")

#remove trisomy 12s
ref_hap_spec_CN <- ref_hap_spec_CN[-intersect(which(ref_hap_spec_CN$chr == "12"), which(ref_hap_spec_CN$ID %in% ctrl_tri12_cells)),]

#save ref TPM values
saveRDS(ref_hap_spec_CN, file = sprintf("%s/ref_hap_spec_CN.rds", dirpath))

################################
### Plot AS_CN distributions ###
################################

#plot the reference distributions for hap specific CN states
hap_spec_CN_distrs <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(data = ref_hap_spec_CN, mapping = aes(x = chr, y = vals, color = hap), outlier.alpha = 0)  +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

###################################################################################################
### Create ref distributions for disomy reference distribution split by arm, haplotype, per chr ###
###################################################################################################

#Copy and adjust ctrl samples
ctrl_samples_arm <- copy(ctrl_samples)
ctrl_samples_arm$chr[ctrl_samples_arm$chr == c("10a")] <- "10"
ctrl_samples_arm <- ctrl_samples_arm[!which(ctrl_samples_arm$chr == "10b"),]
ctrl_samples_arm$chr[which(ctrl_samples_arm$chr == "23")] <- "X"

#stack ctrl samples matrix and add column with arms
ctrl_samples_arms <- append(paste(ctrl_samples_arm$chr, "p", sep = ""), paste(ctrl_samples_arm$chr, "q", sep = ""))
ctrl_samples_byarm <- rbind(ctrl_samples_arm, ctrl_samples_arm)
ctrl_samples_byarm$chr <- ctrl_samples_arms

#change colnames for arms to match format of function
colnames(AS_TPM_byarm$A)[1] <- "chr"
colnames(AS_TPM_byarm$B)[1] <- "chr"

#reduce the matrix to arms in the AS_TPM data
ctrl_samples_byarm <- ctrl_samples_byarm[which(ctrl_samples_byarm$chr %in% unique(AS_TPM_byarm$A$chr)),]

#run get_vals_from_mat function to get AS_TPM (Hap A) values for each ref cell
ctrl_AS_TPM_byarm_A <- get_vals_from_mat(samples_df=ctrl_samples_byarm, mat=AS_TPM_byarm$A, add_chr = F)

#run get_vals_from_mat function to get AS_TPM (Hap B) values for each ref cell
ctrl_AS_TPM_byarm_B <- get_vals_from_mat(samples_df=ctrl_samples_byarm, mat=AS_TPM_byarm$B, add_chr = F)

#Get AS_TPM from correct haplotypes
ctrl_AS_TPM_byarm_AB <- cbind(ctrl_AS_TPM_byarm_A[,1:3], A = as.numeric(unlist(ctrl_AS_TPM_byarm_A[,4])), B = as.numeric(unlist(ctrl_AS_TPM_byarm_B[,4])))
ctrl_AS_TPM_byarm <- cbind(ctrl_AS_TPM_byarm_AB[,1:3], vals = rowMeans(as.matrix(ctrl_AS_TPM_byarm_AB[,4:5])))

#bind all together
ref_AS_TPM_byarm <- ctrl_AS_TPM_byarm
ref_AS_TPM_byarm$CN <- as.numeric(ref_AS_TPM_byarm$CN) - 1

#save ref TPM values
saveRDS(ref_AS_TPM_byarm, file = sprintf("%s/ref_AS_TPM_byarm.rds", dirpath))

####################################################
### Visualize distribution of AS_TPM byarm byhap ###
####################################################

ctrl_AS_TPM_byarm_AB_melted <- melt(ctrl_AS_TPM_byarm_AB, measure.vars = c("A", "B"))
ctrl_AS_TPM_byarm_AB_melted$value[which(ctrl_AS_TPM_byarm_AB_melted$chr == "10q")] <- ctrl_AS_TPM_byarm_AB_melted$value[which(ctrl_AS_TPM_byarm_AB_melted$chr == "10q")]*1.5
ctrl_AS_TPM_byarm_AB_melted$value[which(ctrl_AS_TPM_byarm_AB_melted$chr %in% c("Xp", "Xq"))] <- ctrl_AS_TPM_byarm_AB_melted$value[which(ctrl_AS_TPM_byarm_AB_melted$chr %in% c("Xp", "Xq"))]/(mean(ctrl_AS_TPM_byarm_AB_melted[intersect(which(ctrl_AS_TPM_byarm_AB_melted$chr %in% c("Xp", "Xq")), which(ctrl_AS_TPM_byarm_AB_melted$variable == "B")),]$value))

arm_hap_spec_distrs <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_boxplot(data = ctrl_AS_TPM_byarm_AB_melted, mapping = aes(x = chr, y = value, color = variable), outlier.alpha = 0)  +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

#Archived version of AS-CN ref distribuions
#####################################################
### Create ref distributions for each AS-CN state ###
#####################################################

#run get_vals_from_mat function to get AS_TPM (Hap A) values for each ref cell
mono_AS_TPM_A <- get_vals_from_mat(samples_df=mono_samples, mat=AS_TPM_bychr$A)
ctrl_AS_TPM_A <- get_vals_from_mat(samples_df=ctrl_samples, mat=AS_TPM_bychr$A)
tri_AS_TPM_A <- get_vals_from_mat(samples_df=tri_samples, mat=AS_TPM_bychr$A)

#run get_vals_from_mat function to get AS_TPM (Hap B) values for each ref cell
mono_AS_TPM_B <- get_vals_from_mat(samples_df=mono_samples, mat=AS_TPM_bychr$B)
ctrl_AS_TPM_B <- get_vals_from_mat(samples_df=ctrl_samples, mat=AS_TPM_bychr$B)
tri_AS_TPM_B <- get_vals_from_mat(samples_df=tri_samples, mat=AS_TPM_bychr$B)

#Get AS_TPM from correct haplotypes
mono_AS_TPM_AB <- cbind(mono_AS_TPM_A[,1:3], A = as.numeric(unlist(mono_AS_TPM_A[,4])), B = as.numeric(unlist(mono_AS_TPM_B[,4])))
mono_AS_TPM <- cbind(mono_AS_TPM_AB[,1:3], aneuploid = rowMins(as.matrix(mono_AS_TPM_AB[,4:5])), normal = rowMaxs(as.matrix(mono_AS_TPM_AB[,4:5])))
ctrl_AS_TPM_AB <- cbind(ctrl_AS_TPM_A[,1:3], aneuploid = as.numeric(unlist(ctrl_AS_TPM_A[,4])), normal = as.numeric(unlist(ctrl_AS_TPM_B[,4])))
ctrl_AS_TPM_AB_noX_no10 <- ctrl_AS_TPM_AB[-which(ctrl_AS_TPM_AB$chr %in% c("23", "10b")),] #remove these as they are not diploid CN state
tri_AS_TPM_AB <- cbind(tri_AS_TPM_A[,1:3], A = as.numeric(unlist(tri_AS_TPM_A[,4])), B = as.numeric(unlist(tri_AS_TPM_B[,4])))
tri_AS_TPM <- cbind(tri_AS_TPM_AB[,1:3], aneuploid = rowMaxs(as.matrix(tri_AS_TPM_AB[,4:5])), normal = rowMins(as.matrix(tri_AS_TPM_AB[,4:5])))

#bind all together
ref_AS_TPM <- rbind(mono_AS_TPM, ctrl_AS_TPM_AB_noX_no10, tri_AS_TPM)
ref_AS_TPM$CN <- as.numeric(ref_AS_TPM$CN) - 1
ref_AS_TPM_melt <- melt(ref_AS_TPM, measure.vars = c(colnames(ref_AS_TPM)[-c(1:3)]))

#save ref TPM values
saveRDS(ref_AS_TPM_melt, file = sprintf("%s/ref_AS_TPM.rds", dirpath))
saveRDS(ref_AS_TPM, file = sprintf("%s/ref_AS_TPM_unmelt.rds", dirpath))

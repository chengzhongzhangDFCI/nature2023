##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

#################################
### Read in ASE and TPM files ###
#################################

#Read in data binned by genomic coordinates
TPM_nolim <- readRDS(file=sprintf("%s/TPM.nolim.rds", dirpath))
TPM <- readRDS(file=sprintf("%s/TPM.bygene.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/ASE.bygene.rds", dirpath))

#Control samples
controlSampleIDs <- readRDS(sprintf("%s/controlSampleIDs.rds", dirpath))

#######################################
### Function to normalize a matrix ####
#######################################

norm_data_table <- function(mat, anno_cols = 6) {
   #Copy dataframe and remove anno columns
   cols <- c(1:anno_cols)
   dt2go <- copy(data.table(mat)[,-..cols])
   
   #Divide by row means to get ratio to normal 
   rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
   rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), rowmeans)))
   dt3go <- (dt2go / rowmeans_mat_reciprocal)
   dt3go[dt3go == "NaN"] <- NA
   dt3go[dt3go == Inf] <- NA
   
   #add annotations to ratios
   ratios <- cbind(data.table(mat)[,..cols], dt3go)
   
   #return normed ratios
   return(ratios)
}

##################################
### aggregate TPM whole genome ###
##################################

adjust_by_inv_var <- function(mat, anno_cols = 6, w_mat = mat) {
   
   #read in geneRanges
   geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", dirpath))
   
   #Ensure that chr columns in character and sort mat
   mat <- data.table(mat)
   mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
   setkey(mat, chr, start, end)
   
   #Ensure that chr columns in character and sort weighting mat
   w_mat <- data.table(w_mat)
   w_mat$chr <- as.character(w_mat$chr) #turn factor chrs to character chrs
   setkey(w_mat, chr, start, end)
   
   #get the inv variance of each row #TODO: only get var over control cells?
   cols <- c(1:anno_cols)
   gene_inv_variance <- 1/rowVars(as.matrix(w_mat[,..controlSampleIDs]), na.rm = T)
   gene_inv_variance[gene_inv_variance == Inf] <- 0 #In this version points with zero variance are not considered.
   gene_inv_variance[is.na(gene_inv_variance)] <- 0
   
   #merge weights with mat for aggregation
   mat_w <- data.table(cbind(weight = gene_inv_variance, w_mat))
   
   #aggregate with mean weighted by inv variance
   cols <- c(1:(anno_cols+1))
   sum_cols <- colnames(mat_w[,-..cols])
   mat_w_cell <- data.table(mat_w %>% summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
   
   #adjust the input matrix by the adjustment factor determined from weighting mat
   cols <- c(1:anno_cols)
   adjustment_mat <- data.table(t((replicate(nrow(mat), unlist(mat_w_cell)))))
   mat_adjusted <- (mat[,-..cols] / adjustment_mat)
   mat_out <- cbind(mat[,..cols], mat_adjusted)
   
   return(mat_out)
}


adjust_by_inv_var_2 <- function(mat, anno_cols = 6) {
   
   #read in geneRanges
   geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", dirpath))
   
   #Ensure that chr columns in character and sort mat
   mat <- data.table(mat)
   mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
   setkey(mat, chr, start, end)
   
   #get the inv variance of each row #TODO: only get var over control cells?
   cols <- c(1:anno_cols)
   gene_inv_variance <- 1/rowVars(as.matrix(mat[,..controlSampleIDs]), na.rm = T)
   
   #merge weights with mat for aggregation
   mat_w <- data.table(cbind(weight = gene_inv_variance, mat))
   mat_w$weight[is.na(mat_w$weight)] <- 0
   
   #remove chrX
   mat_w_noX <- mat_w %>% filter(chr != "chrX")
   
   #address boundaries of variance calculations
   gene_inv_variance_noinf <- mat_w_noX$weight[which(mat_w_noX$weight != Inf)]
   w_th <- quantile(probs = c(.95), gene_inv_variance_noinf, na.rm = T)
   mat_w$weight[mat_w$weight > w_th] <- w_th #cap weights at 95% quantile excluding chrX
   
   #aggregate with mean weighted by inv variance
   cols <- c(1:(anno_cols+1))
   sum_cols <- colnames(mat_w[,-..cols])
   mat_w_cell <- data.table(mat_w %>% summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
   
   #adjust the input matrix by the adjustment factor determined from weighting mat
   cols <- c(1:anno_cols)
   adjustment_mat <- data.table(t((replicate(nrow(mat), unlist(mat_w_cell)))))
   mat_adjusted <- (mat[,-..cols] / adjustment_mat)
   mat_out <- cbind(mat[,..cols], mat_adjusted)
   
   return(mat_out)
}

########################################################################
### Split chr10 into two parts by adding additional chr factor level ###
########################################################################

split_chr10_anno <- function(mat){
   
   #change chr column to characters instead of factors
   mat$chr <- as.character(unlist(mat$chr))

   #get genes before and after 61Mb
   mat_10 <- mat[which(mat$chr == "chr10")]
   mat_10a_genes <- mat_10$gene_id[which(mat_10$start < 61000000)]
   mat_10b_genes <- mat_10$gene_id[which(mat_10$start >= 61000000)]

   #change regions before 61Mb to a and after 61Mb to b
   mat$chr[which(mat$gene_id %in% mat_10a_genes)] <- "chr10a"
   mat$chr[which(mat$gene_id %in% mat_10b_genes)] <- "chr10b"
   
   return(mat)
   
}

################################
### Norm all data structures ###
################################

#read in geneRanges to change chr entries
geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos.rds", dirpath))

#split chr10 into two chr entries
geneRanges2 <- split_chr10_anno(geneRanges)
saveRDS(geneRanges2, sprintf("%s/geneRanges_Nikos_2.rds", dirpath))

#split chr10 into two chr entries
TPM_nolim2 <- split_chr10_anno(TPM_nolim)
saveRDS(TPM_nolim2, file=sprintf("%s/TPM.nolim.rds", dirpath))

#split chr10 into two chr entries
ASE2 <- list()
ASE2$A <- split_chr10_anno(ASE$A)
ASE2$B <- split_chr10_anno(ASE$B)
ASE2$AF <- split_chr10_anno(ASE$AF)
ASE2$TC <- split_chr10_anno(ASE$TC)
saveRDS(ASE2, file=sprintf("%s/ASE.bygene.rds", dirpath))

#split chr10 into two chr entries
TPM2 <- split_chr10_anno(TPM)

#calculate TPM ratios and do global adjustment
TPM_normed <- norm_data_table(TPM2)
TPM_adjusted <- adjust_by_inv_var(TPM_normed)
saveRDS(TPM_adjusted, file=sprintf("%s/TPM.bygene.normed.rds", dirpath))



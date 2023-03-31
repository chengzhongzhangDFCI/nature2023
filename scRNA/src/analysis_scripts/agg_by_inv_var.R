##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

#################################
### Read in ASE and TPM files ###
#################################

TPM_normed <- readRDS(file=sprintf("%s/TPM.bygene.normed.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/ASE.bygene.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/controlSampleIDs.rds", dirpath))

##########################################################################
### Function to aggregate with inverse variance by genomic coordinates ###
##########################################################################

agg_by_tpm_inv_var <- function(genomic_region = 10000000, mat, anno_cols = 6) {
  
  #read in geneRanges
  geneRanges <- readRDS(sprintf("%s/geneRanges_Nikos_2.rds", dirpath))
  
  #Ensure that chr columns in character and sort mat
  mat <- data.table(mat)
  mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
  setkey(mat, chr, start, end)
  
  #Get the max length of each chromosome
  chr_lengths <- c()
  for (i in unique(geneRanges$chr)) {
    chr_lengths <- append(chr_lengths, max(geneRanges$end[which(geneRanges$chr == i)]))
  }
  chr_seq_all <- data.table(); ind <- 0;
  for (i in 1:length(unique(geneRanges$chr))){
    num_bins <- floor(chr_lengths[i]/genomic_region)
    chr_seq <- seq(0, (num_bins*genomic_region), by = genomic_region)
    chr_seq[length(chr_seq)] <- chr_lengths[i]
    ind <- ind + length(chr_seq) - 1
    chr_seq_mat <- cbind(chr = as.character(unique(geneRanges$chr)[i]), cbind(start = chr_seq[1:(length(chr_seq)-1)], end = chr_seq[2:length(chr_seq)]))
    chr_seq_all <- rbind(chr_seq_all, chr_seq_mat)
  }
  chr_seq_all$start <- as.numeric(chr_seq_all$start); chr_seq_all$end <- as.numeric(chr_seq_all$end); 
  chr_seq_all <- cbind(bin = c(1:ind), chr_seq_all)
  setkey(chr_seq_all, chr, start, end)
  
  #Prep data table for foverlaps
  changeCols <- c("start", "end") #stores (non-anno) columns that need to be converted to numeric values
  mat[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
  setkey(mat, chr, start, end)
  
  #foverlaps mat with genomic coordinate bins
  mat_bin <- foverlaps(chr_seq_all, mat, nomatch=NA) 
  
  #Relabel columns and make bin numeric
  setnames(mat_bin, old=c("start", "end", "i.start", "i.end"), new=c("gene_start", "gene_end", "bin_start", "bin_end"))
  mat_bin <- cbind(mat_bin[,c("chr", "bin", "bin_start", "bin_end", "gene_id", "gene_start", "gene_end", "width", "strand")], mat_bin[,-c("chr", "bin", "bin_start", "bin_end", "gene_id", "gene_start", "gene_end", "width", "strand")])
  
  #get the inv variance of each row #TODO: only get var over control cells?
  gene_inv_variance <- 1/rowVars(as.matrix(mat_bin[,..controlSampleIDs]), na.rm = T)
  
  #merge weights with mat for aggregation
  mat_w <- data.table(cbind(weight = gene_inv_variance, mat_bin))
  mat_w$weight[is.na(mat_w$weight)] <- 0
  
  #remove chrX
  mat_w_noX <- mat_w %>% filter(chr != "chrX")
  
  #address boundaries of variance calculations
  gene_inv_variance_noinf <- mat_w_noX$weight[-which(mat_w_noX$weight == Inf)]
  w_th <- quantile(probs = c(.95), gene_inv_variance_noinf)
  mat_w$weight[mat_w$weight > w_th] <- w_th #cap weights at 95% quantile excluding chrX
  
  #aggregate with mean weighted by inv variance
  cols <- c(1:(anno_cols+4))
  sum_cols <- colnames(mat_w[,-..cols])
  mat_w_binned <- data.table(mat_w %>% 
    group_by(bin) %>% 
    summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
  
  #merge mat with annos by bin
  setkey(chr_seq_all, bin)
  setkey(mat_w_binned, bin)
  mat_binned2 <- merge(chr_seq_all, mat_w_binned)
  mid <- (mat_binned2$start + mat_binned2$end) / 2
  mat_binned3 <- cbind(mat_binned2[,c(1:4)], cbind(mid, mat_binned2[,-c(1:4)]))
  
  #make sure all value columns are numeric
  changeCols <- colnames(mat_binned2)[-c(1:5)]
  mat_binned3[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
  
  return(mat_binned3)
}

############################################
### aggregate ASE by genomic coordinates ###
############################################

#use function to bin each 
ASE_bybin <- list()
ASE_bybin$AF <- agg_by_tpm_inv_var(mat = ASE$AF, genomic_region = 10000000)

saveRDS(ASE_bybin, file=sprintf("%s/ASE.inv_var.bybin.rds", dirpath))

############################################
### aggregate TPM by genomic coordinates ###
############################################

#use function to bin each 
TPM_bybin <- agg_by_tpm_inv_var(mat = TPM_normed, genomic_region = 10000000)

saveRDS(TPM_bybin, file=sprintf("%s/TPM.inv_var.bybin.rds", dirpath))

############################
#### Calculate AS-TPM v2 ###
############################

AS_TPM_bybin <- list()
AS_TPM_bybin$A <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * (2*ASE_bybin$AF[,-c(1:5)]) ))
AS_TPM_bybin$B <- cbind(ASE_bybin$AF[,c(1:5)], (TPM_bybin[,-c(1:5)] * (2*(1-ASE_bybin$AF[,-c(1:5)])) ))

saveRDS(AS_TPM_bybin, file=sprintf("%s/AS-TPM.inv_var.bybin.rds", dirpath))

###################################################################
### Function to aggregate with inverse variance over chromosome ###
###################################################################

agg_inv_var_byarm_or_chr <- function(mat, size = "chr", anno_cols = 6, simpleX = F) {
  
  #Ensure that chr columns in character and sort mat
  mat <- data.table(mat)
  mat$chr <- as.character(mat$chr) #turn factor chrs to character chrs
  setkey(mat, chr, start, end)
  
  #If size is arm read in armRanges and merge anno
  if(size == "arm") {
    
    message("aggregating matrix over chr arms based on variance")
    
    #read in arm annotations
    armRanges <- data.table(readRDS(file = sprintf("%s/armRanges.rds", dirpath)))
    setnames(armRanges, old = c("seqnames"), new = c("chr"))
    setkey(armRanges, chr, start, end)
    
    #foverlaps chr and arm annotations
    mat_arm <- foverlaps(mat, armRanges)
    cols <- c("chr", "start", "end", "width", "strand")
    mat_arm <- mat_arm[,-..cols]
    setnames(mat_arm, old = c("i.start", "i.end", "i.width", "i.strand"), new = c("start", "end", "width", "strand"))
    
    #change the order of the columns to match orginal form
    mat <- data.table(cbind(mat_arm[,c("gene_id", "arm", "start", "end", "width", "strand")], mat_arm[,-c("gene_id", "arm", "start", "end", "width", "strand")]))
    
    #get the inv variance of each row #TODO: only get var over control cells?
    gene_inv_variance <- 1/rowVars(as.matrix(mat[,..controlSampleIDs]), na.rm = T)
    
    #merge weights with mat for aggregation
    mat_w <- data.table(cbind(weight = gene_inv_variance, mat))
    mat_w$weight[is.na(mat_w$weight)] <- 0
    
    #remove chrX
    #mat_w_noX <- mat_w %>% filter(chr != "chrX")
    mat_w_noX <- mat_w %>% filter(!arm %in% c("Xp", "Xq")) ### CHANGE ###
    
    
    #address boundaries of variance calculations
    gene_inv_variance_noinf <- mat_w_noX$weight[-which(mat_w_noX$weight == Inf)]
    w_th <- quantile(probs = c(.95), gene_inv_variance_noinf)
    mat_w$weight[mat_w$weight > w_th] <- w_th #cap weights at 95% quantile excluding chrX
    
    #aggregate with mean weighted by inv variance
    cols <- c(1:(anno_cols+1))
    sum_cols <- colnames(mat_w[,-..cols])
    mat_byarm <- data.table(mat_w %>% 
                              group_by(arm) %>% 
                              summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))
    if(sum(is.na(mat_byarm$arm)) > 0) {mat_byarm_test <- mat_byarm[-is.na(mat_byarm$arm),]}
    
    if(simpleX) {
      #aggregate with simple mean
      mat_byarm_simpleX <- data.table(mat_w %>% 
                                        group_by(arm) %>% 
                                        summarise_at(sum_cols, funs(mean(., na.rm = T))))
      if(sum(is.na(mat_byarm_simpleX$arm)) > 0) {mat_byarm_simpleX <- mat_byarm_simpleX[-is.na(mat_byarm_simpleX$arm),]}

      #Use weighted mean arms but simple mean arms from chrX.
      mat_byarm <- rbind(mat_byarm[-which(mat_byarm$arm %in% c("Xp", "Xq"))], mat_byarm_simpleX[which(mat_byarm_simpleX$arm %in% c("Xp", "Xq"))])
    }
    
    #make sure all value columns are numeric
    changeCols <- colnames(mat_byarm)[-c(1)]
    mat_byarm[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
    
    #return agg by chr matrix
    return(mat_byarm)
  }
  
  if(size == "chr") {
    
    message("aggregating matrix over chrs based on variance")
    
    #get the inv variance of each row 
    cols <- c(1:anno_cols)
    gene_inv_variance <- 1/rowVars(as.matrix(mat[,..controlSampleIDs]), na.rm = T)

    #merge weights with mat for aggregation
    mat_w <- data.table(cbind(weight = gene_inv_variance, mat))
    mat_w$weight[is.na(mat_w$weight)] <- 0
    
    #remove chrX
    mat_w_noX <- mat_w %>% filter(chr != "chrX")
    
    #address boundaries of variance calculations
    gene_inv_variance_noinf <- mat_w_noX$weight[-which(mat_w_noX$weight == Inf)]
    w_th <- quantile(probs = c(.95), gene_inv_variance_noinf)
    mat_w$weight[mat_w$weight > w_th] <- w_th #cap weights at 95% quantile excluding chrX

    #aggregate with mean weighted by inv variance
    cols <- c(1:(anno_cols+1))
    sum_cols <- colnames(mat_w[,-..cols])
    mat_bychr <- data.table(mat_w %>% 
                              group_by(chr) %>% 
                              summarise_at(sum_cols, funs(weighted.mean(., weight, na.rm = T))))

    if(simpleX){
      #aggregate with simple mean
      mat_bychr_simpleX <- data.table(mat_w %>% 
                                        group_by(chr) %>% 
                                        summarise_at(sum_cols, funs(mean(., na.rm = T))))
      #Use weighted mean arms but simple mean arms from chrX.
      mat_bychr <- rbind(mat_bychr[-which(mat_bychr$chr %in% c("chrX"))], mat_bychr_simpleX[which(mat_bychr_simpleX$chr %in% c("chrX"))])
    }
    
    #make sure all value columns are numeric
    changeCols <- colnames(mat_bychr)[-c(1)]
    mat_bychr[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric for allele A
    
    #return agg by chr matrix
    return(mat_bychr)
  }

}

####################################
### aggregate ASE and TPM by chr ###
####################################

#use function to bin by chr 
ASE_bychr <- list()
ASE_bychr$AF <- agg_inv_var_byarm_or_chr(mat = ASE$AF, size = "chr", simpleX = T)
saveRDS(ASE_bychr, file=sprintf("%s/ASE.inv_var.bychr.rds", dirpath))

#use function to bin each 
TPM_bychr <- agg_inv_var_byarm_or_chr(mat = TPM_normed, size = "chr")
saveRDS(TPM_bychr, file=sprintf("%s/TPM.inv_var.bychr.rds", dirpath))

####################################
### aggregate ASE and TPM by arm ###
####################################

#recombine 10a and 10b for arm level analysis
ASE_tmp <- copy(ASE)
ASE_tmp$AF$chr[ASE_tmp$AF$chr %in% c("chr10a", "chr10b")] <- "chr10"
TPM_tmp <- copy(TPM_normed)
TPM_tmp$chr[TPM_tmp$chr %in% c("chr10a", "chr10b")] <- "chr10"

#use function to bin by arm 
ASE_byarm <- list()
ASE_byarm$AF <- agg_inv_var_byarm_or_chr(mat = ASE_tmp$AF, size = "arm", simpleX = T)

#use function to bin each 
TPM_byarm <- agg_inv_var_byarm_or_chr(mat = TPM_tmp, size = "arm")

#only take arms covered by both AF and TPM values
covered_arms <- intersect(ASE_byarm$AF$arm, TPM_byarm$arm)
covered_arms_noNA <- covered_arms[!is.na(covered_arms)]
ASE_byarm$AF <- ASE_byarm$AF[which(ASE_byarm$AF$arm %in% covered_arms_noNA),]
TPM_byarm <- TPM_byarm[which(TPM_byarm$arm %in% covered_arms_noNA),]

#save the objects
saveRDS(ASE_byarm, file=sprintf("%s/ASE.inv_var.byarm.rds", dirpath))
saveRDS(TPM_byarm, file=sprintf("%s/TPM.inv_var.byarm.rds", dirpath))

############################
#### Calculate AS-TPM v2 ###
############################

#calculate AS_TPM by chr
AS_TPM_bychr <- list()
AS_TPM_bychr$A <- data.table(cbind(ASE_bychr$AF[,c(1)], (TPM_bychr[,-c(1)] * (2*ASE_bychr$AF[,-c(1)]) )))
AS_TPM_bychr$B <- data.table(cbind(ASE_bychr$AF[,c(1)], (TPM_bychr[,-c(1)] * (2*(1-ASE_bychr$AF[,-c(1)]))) ))
saveRDS(AS_TPM_bychr, file=sprintf("%s/AS-TPM.inv_var.bychr.rds", dirpath))

#calculate AS_TPM by arm
AS_TPM_byarm <- list()
AS_TPM_byarm$A <- cbind(ASE_byarm$AF[,c(1)], (TPM_byarm[,-c(1)] * (2*ASE_byarm$AF[,-c(1)]) ))
AS_TPM_byarm$B <- cbind(ASE_byarm$AF[,c(1)], (TPM_byarm[,-c(1)] * (2*(1-ASE_byarm$AF[,-c(1)]))) )
saveRDS(AS_TPM_byarm, file=sprintf("%s/AS-TPM.inv_var.byarm.rds", dirpath))


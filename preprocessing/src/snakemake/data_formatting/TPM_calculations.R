##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

################
#### Imports ###
################

#Annotion list with information on every cell and IDs of control samples
anno <- data.table(read.csv( sprintf("%s/analysis_data/annotation_list.csv", dirpath)))
anno$Fastq_files <- as.character(anno$Fastq_files)
anno$WTA.plate <- as.character(anno$WTA.plate)
controlSampleIDs <- readRDS( sprintf("%s/analysis_data/controlSampleIDs.rds", dirpath))

#import gene annotations
geneRanges <- readRDS(sprintf("%s/analysis_data/geneRanges_Nikos.rds", dirpath))
setkey(geneRanges, "gene_id")

#import QC info for use in curating control cells
all_QC <- readRDS(file = sprintf("%s/all_QC.rds", dirpath))

###################################
### Read in RSEM output tables ###
##################################

#Generate the directories based on experiments
dir_paths <- list()
for (i in experiments) {dir_paths <- cbind(dir_paths, sprintf("%s/%s/RSEM/output", dirpath, i))}

#list all files in each directory folder
files <- c()
for (i in 1:length(dir_paths)){files <- append(files, paste(dir_paths[i], list.files(path = as.character(dir_paths[i]), pattern = ".genes.results"), sep = "/"))}

#List sample names
files_samples <- files
for (i in 1:length(anno$Fastq_files)){
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T, value=T, x = files)
  if (length(res) < 1) {res <- grep(pattern =  anno$Fastq_files[i], ignore.case = T, value=T, x = files)}
  files_samples[which(files == res[1])] <- anno$WTA.plate[i]
  files_samples[which(files == res[2])] <- sprintf("%s_L2", anno$WTA.plate[i]) #if there is a second file it is because there were two lanes
}
files_samples_v2 <- append(files_samples[which(files_samples %in% anno$WTA.plate)], files_samples[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])
files_v2 <- append(files[which(files_samples %in% anno$WTA.plate)], files[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])

####################################
### Parallel process RSEM tables ###
####################################

get_tpm_matrix <- function(file, cell_name) {
  
  #save start time
  start_time <- Sys.time()
  
  #Read in input file
  TPM <- data.table(read.table(file, header = T)) #Read-in input file
  if (nrow(TPM) < 1) {return(sprintf("Empty: %s", file))}
  TPM <- TPM[,c("gene_id", "TPM")] #only take relevant columns
  
  #annotate and clean data.table
  counts_wcell <- cbind(cell = cell_name, TPM) #add cell name
  
  #make directory for allelic data if it doesn't exist
  system(sprintf("mkdir -p %s/TPM_data", dirpath, cell_name))
  
  #save data frame with allelic variant coverage
  saveRDS(counts_wcell, sprintf("%s/TPM_data/%s.TPM.rds", dirpath, cell_name))
  
  #save end time
  end_time <- Sys.time()
  
  #how long function took to run
  time <- end_time - start_time
  message(sprintf("Time: %s", time))
  
  #return run-time and file name
  return( sprintf("Completed: %s", file) )
  
}

##############################################
### Run function and aggregate RSEM tables ###
##############################################

require(doParallel)
registerDoParallel(30)

#reform ASE information in parallel
foreach(i = c(1:length(files_v2))) %dopar% { 
  get_tpm_matrix(file = files_v2[i], cell_name = files_samples_v2[i])
}

#concatenate ASE data 
TPM_all_files <- data.frame()
for (i in 1:length(files_samples_v2)) { 
  if (file.exists(sprintf("%s/TPM_data/%s.TPM.rds", dirpath, files_samples_v2[i]))) {
    TPM_cov <- readRDS(sprintf("%s/TPM_data/%s.TPM.rds", dirpath, files_samples_v2[i]))
    TPM_all_files <- rbind(TPM_all_files, TPM_cov)
  }
}

#Save concatentated data
saveRDS(TPM_all_files, sprintf("%s/analysis_data/ASE.concatTPM.rds", dirpath))
TPM_all_files <- readRDS(sprintf("%s/analysis_data/ASE.concatTPM.rds", dirpath))

#cast the long data frame by cell
TPM_all_cells <- dcast(TPM_all_files, gene_id ~ cell, value.var = "TPM")
TPM_all_cells[is.na(TPM_all_cells)] <- 0
TPM_all_cells <- data.table(TPM_all_cells)
setkey(TPM_all_cells, "gene_id")

#merge with geneRanges
rsemtpm <- merge(geneRanges, TPM_all_cells)
setkey(rsemtpm, chr, start, end)
rsemtpm <- rsemtpm[which(!rsemtpm$chr %in% c("chrY", "chrM")),]

######################################
### Merge counts from 2 lane cells ###
######################################

#use grep to get cells with L2 in name
L2_cells <- grep(pattern = "L2", ignore.case = T, value = T, x = colnames(rsemtpm))
cells <- str_remove(L2_cells, pattern = "_L2") #get corresponding L1 cells
L1_L2_cells <- append(cells, L2_cells) #combine cell column names
TPM_2 <- cbind( rsemtpm[,-..L1_L2_cells], ((rsemtpm[,..cells] + rsemtpm[,..L2_cells])/2) ) #avg counts between lanes

#sort columns
col_order <- append(colnames(TPM_2)[c(1:6)], sort(colnames(TPM_2)[-c(1:6)]))
TPM_3 <- TPM_2[,..col_order]

#####################################
### Take only cells with ASE data ###
#####################################

ASE <- readRDS(file=sprintf("%s/analysis_data/ASE.bygene.rds", dirpath))

cols <- intersect(colnames(ASE$AF), colnames(TPM_3))

TPM <- TPM_3[,..cols] 

#save TPM intermediate file
saveRDS(TPM, file=sprintf("%s/analysis_data/raw_TPM.rds", dirpath))
write.csv(TPM, file=sprintf("%s/analysis_data/raw_TPM.csv", dirpath), row.names = F)
#TPM <- readRDS(file=sprintf("%s/aggregated_results/raw_TPM.rds", dirpath))

#################################
### Remove mono-allelic genes ###
#################################

TMP_AF_DF <- ASE$AF[!which(ASE$AF$chr %in% c("chr10", "chrX")),]
avg_AF_bygene <- rowMeans(TMP_AF_DF[,..controlSampleIDs], na.rm = T)
monoallelic_genes <- TMP_AF_DF$gene_id[-which(between(avg_AF_bygene, .3, .7))]

#special case for 10q where there is a third copy translocated to chrX
TMP_AF_10a_df <- ASE$AF[which(ASE$AF$chr %in% c("chr10")),]
TMP_AF_10a_df <- TMP_AF_10a_df[which(TMP_AF_10a_df$start < 61000000),]
avg_AF_10a_bygene <- rowMeans(TMP_AF_10a_df[,..controlSampleIDs], na.rm = T)
monoallelic_genes_10a <- TMP_AF_10a_df$gene_id[-which(between(avg_AF_10a_bygene, .3, .7))]

TMP_AF_10b_df <- ASE$AF[which(ASE$AF$chr %in% c("chr10")),]
TMP_AF_10b_df <- TMP_AF_10b_df[which(TMP_AF_10b_df$start >= 61000000),]
avg_AF_10b_bygene <- rowMeans(TMP_AF_10b_df[,..controlSampleIDs], na.rm = T)
monoallelic_genes_10b <- TMP_AF_10b_df$gene_id[-which(between(avg_AF_10b_bygene, .2, .4))]

#combine all monoallelic genes
monoallelic_genes_10 <- append(monoallelic_genes_10a, monoallelic_genes_10b)
monoallelic_genes_all <- append(monoallelic_genes, monoallelic_genes_10)

######################################
### Reduce DFs to reduced gene set ###
######################################

TPM <- TPM[!which(TPM$gene_id %in% monoallelic_genes_all),]

ASE$A <- ASE$A[!which(ASE$A$gene_id %in% monoallelic_genes_all),..cols] 
ASE$B <- ASE$B[!which(ASE$B$gene_id %in% monoallelic_genes_all),..cols] 
ASE$AF <- ASE$AF[!which(ASE$AF$gene_id %in% monoallelic_genes_all),..cols] 
ASE$TC <- ASE$TC[!which(ASE$TC$gene_id %in% monoallelic_genes_all),..cols] 

#Save TPM object by gene
TPM$chr <- as.character(unlist(TPM$chr))
saveRDS(TPM, file=sprintf("%s/analysis_data/TPM.nolim.rds", dirpath))
TPM <- readRDS(file=sprintf("%s/analysis_data/TPM.nolim.rds", dirpath))

################################################
### Eliminate genes with very low expression ###
################################################

rowmeans <- cbind(TPM[,c(1:6)], ctrl_means = rowMeans(TPM[,..controlSampleIDs]))
min_exp_genes <- unlist(rowmeans[which(rowmeans$ctrl_means >= 25),]$gene_id) #35 seems to be good medium point
#genes <- intersect(min_exp_genes, intersect(unique(ASE$AF$gene_id), unique(TPM$gene_id)))

TPM <- TPM[which(TPM$gene_id %in% min_exp_genes),]

########################
### Save Data tables ###
########################

#resave reduced ASE data structure
saveRDS(ASE, file=sprintf("%s/analysis_data/ASE.bygene.rds", dirpath))

#Save TPM object by gene
saveRDS(TPM, file=sprintf("%s/analysis_data/TPM.bygene.rds", dirpath))


############################################################################################################
#End of Script



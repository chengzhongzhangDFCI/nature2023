##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

#-------------------------------------------------------
#Read in and aggregate results from ASEReadCounter

########################################
### Read in Hets file and Phase file ###
########################################

#read in anno list
anno <- data.table(read.csv( sprintf("%s/analysis_data/annotation_list.csv", dirpath)))
anno$Fastq_files <- as.character(anno$Fastq_files)
anno$WTA.plate <- as.character(anno$WTA.plate)

#Read in centromere data and get ranges for each chr
centromeres <- readRDS( sprintf("%s/analysis_data/centromeres.rds", dirpath) )

#read in gene annotation gtf
gtf <- rtracklayer::import(sprintf("%s/../../data/refgenomes/Gencode/v25/gencode.v25.primary_assembly.annotation.gtf", scriptsdir)) #Adjust to appropriate location
gtf_df=as.data.table(gtf)
setnames(gtf_df, old = c("seqnames"), new = c("chr"))
setkey(gtf_df, "chr", "start", "end", "gene_id")

#Get genes ref df and exon/UTR df
geneRanges <- gtf_df[which(gtf_df$type == "gene"),c("gene_id", "chr", "start", "end", "width", "strand")]
setkey(geneRanges, "chr", "start", "end")
chrs <- paste0("chr", append(c(1:22), "X"))
geneRanges <- geneRanges[which(geneRanges$chr %in% chrs),]
coding_regions <- gtf_df[which(gtf_df$type %in% c("UTR", "exon")),c("chr", "start", "end", "type", "gene_type")]
setkey(coding_regions, "chr", "start", "end")
geneAnnos <- foverlaps(coding_regions, geneRanges)
setnames(geneAnnos, old = c("start", "end", "i.start", "i.end"), new = c("gene_start", "gene_end", "start", "end"))

#read in geneRanges
setkey(geneRanges, "chr", "start", "end")
setkey(geneAnnos, "chr", "start", "end")
saveRDS(geneRanges, sprintf("%s/analysis_data/geneRanges_Nikos.rds", dirpath))

#read in phase from Greg Brunette
phase_file <- sprintf("%s/../../data/RPE1_refs/RPE1_Haplotype_update.dat", scriptsdir)
hets_phase <- data.table(read_table(phase_file))

#Filter hets file by Greg's recommended parameters sample_count > 10 and MAF > 0.348 
hets_phase <- hets_phase[which(hets_phase$MAF > .348),]
hets_phase <- hets_phase[which(hets_phase$sample_count > 10),]

#Remove indeterminant sites from the beginning and extra columns
hets_phase <- hets_phase[which(!hets_phase$Haplotype == 0),]
hets_phase <- hets_phase[,c("chr", "pos", "Haplotype")]
setkey(hets_phase, "chr", "pos", "Haplotype")


#Generate the directories based on experiments
dir_paths <- list()
for (i in experiments) {dir_paths <- cbind(dir_paths, sprintf("%s/%s/ASE", dirpath, i))}

#list all files in each directory folder
files <- c()
for (i in 1:length(dir_paths)){files <- append(files, paste(dir_paths[i], list.files(path = as.character(dir_paths[i]), pattern = ".table"), sep = "/"))}

#List sample names
files_samples <- files
message(head(files_samples))
for (i in 1:length(anno$Fastq_files)){ 
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T, value=T, x = files)
  if (length(res) < 1) {res <- grep(pattern =  anno$Fastq_files[i], ignore.case = T, value=T, x = files)}
  files_samples[which(files == res[1])] <- anno$WTA.plate[i]
  files_samples[which(files == res[2])] <- sprintf("%s_L2", anno$WTA.plate[i])
}
files_samples_v2 <- append(files_samples[which(files_samples %in% anno$WTA.plate)], files_samples[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])
files_v2 <- append(files[which(files_samples %in% anno$WTA.plate)], files[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])

### TMP ### Remove 180201-6A_S7_L001 because it only has 4 SNPs
files_samples_v2_to_remove <- grep(pattern = "180201_6A", ignore.case = T, value=T, x = files_samples_v2)
files_v2_to_remove <- grep(pattern = "180201-6A_S7_L001", ignore.case = T, value=T, x = files_v2)
files_samples_v2 <- files_samples_v2[-which(files_samples_v2 %in% files_samples_v2_to_remove)]
files_v2 <- files_v2[-which(files_v2 %in% files_v2_to_remove)]

##################################################
### Function for data engineering each sample ###
##################################################

get_allelic_var_matrix <- function(file, cell_name, hets_phase_info) {
  
  #save start time
  start_time <- Sys.time()
  
  #Read in input file
  ASE <- data.table(read.table(file, header = T)) #Read-in input file
  if (nrow(ASE) < 1) {return(sprintf("Empty: %s", file))}
  setnames(ASE, old = c("contig", "position"), c("chr", "pos")) #Change column names
  ASE <- ASE[,c("chr", "pos", "refCount", "altCount")] #only take relevant columns
  
  #annotate and clean data.table
  setkey(ASE, "chr", "pos") #set columns to merge
  counts <- merge(hets_phase_info, ASE) #merge with hets phase information
  counts_hets <- counts[!which(counts$Haplotype == 0),] #remove indeterminate sites
  counts_wcell <- cbind(cell = cell_name, counts_hets) #add cell name
  
  #use haplotype and ref/alt coverage to get determine allelic coverage
  Allelic_var_cov <- c() #start new data frame
  rows <- nrow(counts_wcell) #save number of rows for progress calculations
  for (i in 1:nrow(counts_wcell)) { #iterate through cell's SNPs
    tmp <- counts_wcell[i] #store individual line to save time on calculations
    if (tmp$Haplotype == 1) { #+1 haplotype annotation means ref count is allele A
      Allelic_var_cov <- rbind(Allelic_var_cov, cbind(tmp[,c("cell", "chr", "pos")], cbind(allele = "A", counts = tmp$refCount))) #add allele annotation with ref as A
      Allelic_var_cov <- rbind(Allelic_var_cov, cbind(tmp[,c("cell", "chr", "pos")], cbind(allele = "B", counts = tmp$altCount))) #add allele annotation with alt as B
    }
    if (tmp$Haplotype == -1) { #-1 haplotype annotation means ref count is allele B
      Allelic_var_cov <- rbind(Allelic_var_cov, cbind(tmp[,c("cell", "chr", "pos")], cbind(allele = "B", counts = tmp$refCount))) #add allele annotation with ref as B
      Allelic_var_cov <- rbind(Allelic_var_cov, cbind(tmp[,c("cell", "chr", "pos")], cbind(allele = "A", counts = tmp$altCount))) #add allele annotation with alt as A
    }
    if((i %% 5000) == 0) {message(sprintf("%s %s done", round(i/rows, digits = 2), "%"))} #report progress
  }
  
  #make directory for allelic data if it doesnt exist
  system(sprintf("mkdir -p %s/ASE_data", dirpath, cell_name))
  
  #save data frame with allelic variant coverage
  saveRDS(Allelic_var_cov, sprintf("%s/ASE_data/%s.ASE.rds", dirpath, cell_name))
  
  #save end time
  end_time <- Sys.time()
  
  #how long function took to run
  time <- end_time - start_time
  message(sprintf("Time: %s", time))
  
  #return run-time and file name
  return( sprintf("Completed: %s", file) )
  
}

#####################################################
### Run ASE function for each ASEReadCounter file ###
#####################################################

require(doParallel)
registerDoParallel(20)

#reform ASE information in parallel
foreach(i = c(1:length(files_v2))) %dopar% { 
  get_allelic_var_matrix(file = files_v2[i], cell_name = files_samples_v2[i], hets_phase_info = hets_phase)
}

#concatenate ASE data 
ASE_all_files <- data.frame()
for (i in 1:length(files_samples_v2)) { 
  if (file.exists(sprintf("%s/ASE_data/%s.ASE.rds", dirpath, files_samples_v2[i]))) {
    Allelic_var_cov <- readRDS(sprintf("%s/ASE_data/%s.ASE.rds", dirpath, files_samples_v2[i]))
    ASE_all_files <- rbind(ASE_all_files, Allelic_var_cov)
  }
}

#cast the long data frame by cell
ASE_all_cells <- dcast(ASE_all_files, chr + pos + allele ~ cell, value.var = "counts")
saveRDS(ASE_all_cells, file=sprintf("%s/analysis_data/ASE_all_cells.rds", dirpath))
ASE_all_cells <- readRDS(file=sprintf("%s/analysis_data/ASE_all_cells.rds", dirpath))
ASE_all_cells[is.na(ASE_all_cells)] <- 0
ASE_all_cells <- data.table(ASE_all_cells)
setkey(ASE_all_cells, "chr", "pos")

#TEMP: ONLY KEEP HETS FILTER HETS #REMOVES 8 SITES, WHY?
hets_phase <- hets_phase[,c("chr", "pos")]
setkey(hets_phase, "chr", "pos")
ASE_all_cells <- merge(hets_phase, ASE_all_cells)

#make all count columns numeric
changeCols <- colnames(ASE_all_cells)[-c(1:3)] #stores (non-anno) columns that need to be converted to numeric values
ASE_all_cells[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric 

######################################
### Merge counts from 2 lane cells ###
######################################

#use grep to get cells with L2 in name
L2_cells <- grep(pattern = "L2", ignore.case = T, value = T, x = colnames(ASE_all_cells))
cells <- str_remove(L2_cells, pattern = "_L2") #get corresponding L1 cells
L1_L2_cells <- append(cells, L2_cells) #combine cell column names
L1_L2 <- cbind(cells, L2_cells) #combine cell column names
ASE_all_cells <- cbind( ASE_all_cells[,-..L1_L2_cells], ((ASE_all_cells[,..cells] + ASE_all_cells[,..L2_cells])/2) ) #avg counts between lanes

############################################################################################################
#Clean and structure ASE information by gene

#########################################
### Aggregate SNPs by gene and get AF ###
#########################################

#merge geneRanges with SNPs
setnames(ASE_all_cells, old = "pos", new = "start")
ASE_all_cells$start <- as.numeric(ASE_all_cells$start)
ASE_all_cells$end <- copy(ASE_all_cells$start)
setkey(ASE_all_cells, "chr", "start", "end")

#foverlap geneRanges and SNPs
ASE_bygene <- foverlaps(ASE_all_cells, geneAnnos)

#Relabel columns
setnames(ASE_bygene, old=c("i.start", "i.end"), new=c("SNP_pos", "SNP_pos_end"))
ASE_bygene <- ASE_bygene[,-c("SNP_pos_end")]

#Elimintate SNPs not inside of a gene
ASE_bygene2 <- ASE_bygene[which(!ASE_bygene$gene_id == "NA"),]

#############################################
### Aggregate SNPs by genomic coordinates ###
#############################################

#Change cols in numeric
changeCols <- colnames(ASE_bygene2)[-c(1:12)] #stores (non-anno) columns that need to be converted to numeric values
ASE_bygene2[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols] #converts all columns in "changeCols" to numeric 

#sort columns
col_order <- append(colnames(ASE_bygene2)[c(1:12)], sort(colnames(ASE_bygene2)[-c(1:12)]))
ASE_bygene2 <- ASE_bygene2[,..col_order]
setkey(ASE_bygene2, chr, gene_start, gene_end)

#Split matrix by allele
ASE_A <- data.table(ASE_bygene2[which(ASE_bygene2$allele == "A"),])
ASE_B <- data.table(ASE_bygene2[which(ASE_bygene2$allele == "B"),])

#aggregate SNPs to gene level with simple sum
ASE_A2 <- data.table(ASE_A[,-c(1,3:12)] %>%
  group_by(gene_id) %>%
  summarise_all(sum, na.rm = TRUE))
setkey(ASE_A2, gene_id)
setkey(geneRanges, gene_id)
ASE_A3 <- merge(geneRanges, ASE_A2)
setkey(ASE_A3, chr, start, end)

#aggregate SNPs to gene level with simple sum
ASE_B2 <- data.table(ASE_B[,-c(1,3:12)] %>%
  group_by(gene_id) %>%
  summarise_all(sum, na.rm = TRUE))
setkey(ASE_B2, gene_id)
setkey(geneRanges, gene_id)
ASE_B3 <- merge(geneRanges, ASE_B2)
setkey(ASE_B3, chr, start, end)

ASE <- list()
ASE$A <- ASE_A3 #Allele A SNPs counts across all genes
ASE$B <- ASE_B3 #Allele B SNPs counts across all genes

###################################################
### Calculate Allele Frequency and Total Counts ###
###################################################

#Calc Allele Frequency 
ASE$AF <- cbind(ASE$A[,c(1:6)], (ASE$A[,-c(1:6)] / (ASE$A[,-c(1:6)] + ASE$B[,-c(1:6)])))
ASE$AF[ASE$AF == "NaN"] <- NA #for the 0/0 points

#Calc Total Counts
ASE$TC <- cbind(ASE$A[,c(1:6)],  (ASE$A[,-c(1:6)] + ASE$B[,-c(1:6)]) )

saveRDS(ASE, file=sprintf("%s/analysis_data/ASE.bygene.rds", dirpath))
ASE <- readRDS(file=sprintf("%s/analysis_data/ASE.bygene.rds", dirpath))

############################################################################################################
#End of Script






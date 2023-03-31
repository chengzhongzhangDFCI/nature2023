#automated analysis script

##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

####################
### Read-in data ###
####################

#Read in reference data
ref_AS_TPM <- readRDS(file = sprintf("%s/ref_AS_TPM_unmelt.rds", dirpath))

#read in hap spec CN ref distributions
ref_hap_spec_CN <- readRDS(file = sprintf("%s/ref_hap_spec_CN.rds", dirpath))

#read in normed TPM by chr
TPM_bychr <- readRDS(file=sprintf("%s/TPM.inv_var.bychr.rds", dirpath))

#new ASE
ASE_bychr <- readRDS(file=sprintf("%s/ASE.inv_var.bychr.rds", dirpath))

#new AS-TPM
AS_TPM_bychr <- readRDS(file=sprintf("%s/AS-TPM.inv_var.bychr.rds", dirpath))
VAR_A <- copy(AS_TPM_bychr$A)
VAR_B <- copy(AS_TPM_bychr$B)

#read in allele and chr specific pvalues
AS_TPM_bychr_ctrl_pvals_chr_spec <- readRDS(file = sprintf("%s/pval_AS_TPM_bychr_ctrl_pvals_chr_spec.rds", dirpath))
AS_TPM_bychr_mono_pvals_chr_spec <- readRDS(file = sprintf("%s/pval_AS_TPM_bychr_mono_pvals_chr_spec.rds", dirpath))

#############################
### Adjust pvalues format ###
#############################

#change pvalues into datatable and sort by chr
AS_TPM_bychr_ctrl_pvals_chr_spec$A <- data.table(AS_TPM_bychr_ctrl_pvals_chr_spec$A); AS_TPM_bychr_ctrl_pvals_chr_spec$B <- data.table(AS_TPM_bychr_ctrl_pvals_chr_spec$B);
AS_TPM_bychr_mono_pvals_chr_spec$A <- data.table(AS_TPM_bychr_mono_pvals_chr_spec$A); AS_TPM_bychr_mono_pvals_chr_spec$B <- data.table(AS_TPM_bychr_mono_pvals_chr_spec$B);
setkey(AS_TPM_bychr_mono_pvals_chr_spec$A, chr); setkey(AS_TPM_bychr_mono_pvals_chr_spec$B, chr); 

#Add empty rows for chr10b and chrX so that it fits length of other df
AS_TPM_bychr_ctrl_pvals_chr_spec$A <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$A, t(c("chr10b", rep(NA,(ncol(AS_TPM_bychr_ctrl_pvals_chr_spec$A)-1)))), use.names=FALSE)
AS_TPM_bychr_ctrl_pvals_chr_spec$A <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$A, t(c("chrX", rep(NA,(ncol(AS_TPM_bychr_ctrl_pvals_chr_spec$A)-1)))), use.names=FALSE)
AS_TPM_bychr_ctrl_pvals_chr_spec$B <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$B, t(c("chr10b", rep(NA,(ncol(AS_TPM_bychr_ctrl_pvals_chr_spec$B)-1)))), use.names=FALSE)
AS_TPM_bychr_ctrl_pvals_chr_spec$B <- rbind(AS_TPM_bychr_ctrl_pvals_chr_spec$B, t(c("chrX", rep(NA,(ncol(AS_TPM_bychr_ctrl_pvals_chr_spec$B)-1)))), use.names=FALSE)
setkey(AS_TPM_bychr_ctrl_pvals_chr_spec$A, chr); setkey(AS_TPM_bychr_ctrl_pvals_chr_spec$B, chr); 

#Read in annotation list with sister and cousin information
anno <- data.table(read.csv( sprintf("%s/annotation_list.csv", dirpath)))
changeCols <- colnames(anno)[-c(13:16)]
anno[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols] 

#################################################
### Multiply TPM ratios by factor for base CN ###
#################################################

#calculate x adjustment factor
x_factor <- mean(as.numeric(unlist(VAR_B[which(VAR_B$chr == "chrX"),-c("chr")])), na.rm = T)

#split up matrix into 3 groups
VAR_A_chr10b <- cbind(VAR_A[which(VAR_A$chr == "chr10b"),c("chr")], VAR_A[which(VAR_A$chr == "chr10b"),-c("chr")]*1.5) #cause chr10b in varB should be centered at 2
VAR_A_chrx <- cbind(VAR_A[which(VAR_A$chr == "chrX"),c("chr")], VAR_A[which(VAR_A$chr == "chrX"),-c("chr")]/x_factor) #cause chrX in varB should be centered at 1
VAR_A_rest <- cbind(VAR_A[!which(VAR_A$chr %in% c("chr10b","chrX")),c("chr")], VAR_A[!which(VAR_A$chr %in% c("chr10b","chrX")),-c("chr")]) #no change
VAR_A <- rbind(VAR_A_chr10b, VAR_A_chrx, VAR_A_rest) %>% setkey(chr)

#split up matrix into 3 groups
VAR_B_chr10b <- cbind(VAR_B[which(VAR_B$chr == "chr10b"),c("chr")], VAR_B[which(VAR_B$chr == "chr10b"),-c("chr")]*1.5) #cause chr10b in varB should be centered at 1.5
VAR_B_chrx <- cbind(VAR_B[which(VAR_B$chr == "chrX"),c("chr")], VAR_B[which(VAR_B$chr == "chrX"),-c("chr")]/x_factor) #cause chrX in varB should be centered at 1
VAR_B_rest <- cbind(VAR_B[!which(VAR_B$chr %in% c("chr10b","chrX")),c("chr")], VAR_B[!which(VAR_B$chr %in% c("chr10b","chrX")),-c("chr")]) #no change
VAR_B <- rbind(VAR_B_chr10b, VAR_B_chrx, VAR_B_rest) %>% setkey(chr)

#calculate x adjustment factor
x_factor_total_R <- mean(as.numeric(unlist(VAR_B[which(VAR_B$chr == "chrX"),-c("chr")])), na.rm = T) + mean(as.numeric(unlist(VAR_A[which(VAR_A$chr == "chrX"),-c("chr")])), na.rm = T)

#split up matrix into 3 groups
TPM_bychr_chr10b <- cbind(TPM_bychr[which(TPM_bychr$chr == "chr10b"),c("chr")], TPM_bychr[which(TPM_bychr$chr == "chr10b"),-c("chr")]*3) #cause total expression of chr10b should be centered at 3
TPM_bychr_chrx <- cbind(TPM_bychr[which(TPM_bychr$chr == "chrX"),c("chr")], TPM_bychr[which(TPM_bychr$chr == "chrX"),-c("chr")]*x_factor_total_R) #cause total expression of chrX should be centered at 1 + expresion of silenced allele
TPM_bychr_rest <- cbind(TPM_bychr[!which(TPM_bychr$chr %in% c("chr10b","chrX")),c("chr")], TPM_bychr[!which(TPM_bychr$chr %in% c("chr10b","chrX")),-c("chr")]*2) #no change
TPM_bychr <- rbind(TPM_bychr_chr10b, TPM_bychr_chrx, TPM_bychr_rest) %>% setkey(chr)

#######################
### Strip Plot Data ###
#######################
#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#extract family information
setkey(anno_gen_1_2, WTA.plate)
families <- sort(table(anno_gen_1_2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()

names <- anno_gen_1_2$WTA.plate
ids <- anno_gen_1_2$Cell_IDs
relationship <- anno_gen_1_2$Relationship
families <- anno_gen_1_2$Pairs
family_ids <- anno_gen_1_2$Family_IDs
gen <- anno_gen_1_2$key_pairs
rupt_time <- anno_gen_1_2$MN_rupt_time_simple
event <- anno_gen_1_2$event2
MN_info <- anno_gen_1_2$MN.Daughter
fam_vals <- unlist(rle(families)$values)
fam_reps <- as.numeric(rle(families)$lengths)*(length(unique(TPM_bychr$chr))-length(exclude_chrs))
fam_ids_vals <- unlist(rle(family_ids)$values)
fam_ids_reps <- as.numeric(rle(family_ids)$lengths)*(length(unique(TPM_bychr$chr))-length(exclude_chrs))

chr_num <- str_remove(unique(TPM_bychr$chr), pattern = "chr")
chr_num[chr_num == "X"] <- 23
chr_num[chr_num == "10a"] <- 10.1; chr_num[chr_num == "10b"] <- 10.2;
chr_num <- as.numeric(chr_num)
chr_num <- unique(TPM_bychr$chr)

#Prepare Visual Objects
visual.data <- data.table(TPM = as.numeric(unlist(flatten(TPM_bychr[,..names]))))
visual.data$VAR_A <- as.numeric(unlist(flatten(VAR_A[,..names])))
visual.data$VAR_B <- as.numeric(unlist(flatten(VAR_B[,..names])))
visual.data$chr <- rep(unique(chr_num), length = length(visual.data$TPM))
visual.data$cell <- rep(names, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$ids <- rep(ids, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$family <- rep(fam_vals, fam_reps)
visual.data$family_ids <- rep(fam_ids_vals, fam_ids_reps)
visual.data$relationship <- rep(relationship, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$generation <- rep(gen, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$rupt_time <- rep(rupt_time, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$event <- rep(event, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$MN_info <- rep(MN_info, each=length(unique(chr_num))-length(exclude_chrs))

#######################################################
### Determine threshold for full transcription loss ###
#######################################################

reference_CI <- ref_AS_TPM %>% 
  group_by(CN) %>% 
  summarise(qs_a = quantile(aneuploid, c(.05,.95)), qs_n = quantile(normal, c(.05,.95)), prob = c(.05,.95))

full_loss_th <- unname(reference_CI$qs_a[intersect(which(reference_CI$prob == .95), which(reference_CI$CN == 0))])

###################################################
### Label transcription states based on pvalues ###
###################################################

#flatten allele A pvalues
visual.data$hapA_pval_CN1 <- as.numeric(unlist(flatten(AS_TPM_bychr_mono_pvals_chr_spec$A[,..names])))
visual.data$hapA_pval_CN2 <- as.numeric(unlist(flatten(AS_TPM_bychr_ctrl_pvals_chr_spec$A[,..names]))) #does not match length because no control values for 10b and x

#flatten allele B pvalues
visual.data$hapB_pval_CN1 <- as.numeric(unlist(flatten(AS_TPM_bychr_mono_pvals_chr_spec$B[,..names])))
visual.data$hapB_pval_CN2 <- as.numeric(unlist(flatten(AS_TPM_bychr_ctrl_pvals_chr_spec$B[,..names]))) #does not match length because no control values for 10b and x

#merge cell lineage data with possible events
setkey(visual.data, "chr")
possible_data2 <- visual.data[,c("chr", "relationship", "cell", "ids", "family", "family_ids", "event", "generation", "rupt_time", "MN_info", "TPM", "VAR_A", "VAR_B", "hapA_pval_CN1", "hapB_pval_CN1", "hapA_pval_CN2", "hapB_pval_CN2")]

#set alpha values
alpha <- .05 #normal alpha for trisomy and disomy as those distributions as samples size is 1
alpha_bonf <- .05/(length(unique(visual.data$chr))*2) #Bonferroni corrected pval (due to 45 homologs)

#get state classifications for allele A
gain_table_A <- possible_data2 %>% filter(chr != "chrX") %>% filter(hapA_pval_CN1 < alpha_bonf) %>% filter(VAR_A > 1) %>% bind_cols(., state_A = rep("gain", nrow(.)))
gain_table_A_X <- possible_data2 %>% filter(chr == "chrX") %>% filter(hapA_pval_CN1 < alpha_bonf) %>% filter(VAR_A > .23) %>% bind_cols(., state_A = rep("gain", nrow(.)))
loss_table_A <- possible_data2 %>% filter(chr != "chrX") %>% filter(hapA_pval_CN1 < alpha_bonf) %>% filter(VAR_A < 1) %>% bind_cols(., state_A = rep("loss", nrow(.)))
loss_table_A_X <- possible_data2 %>% filter(chr == "chrX") %>% filter(hapA_pval_CN1 < alpha_bonf) %>% filter(VAR_A < .23) %>% bind_cols(., state_A = rep("loss", nrow(.)))
ctrl_table_A <- possible_data2 %>% filter(hapA_pval_CN1 >= alpha_bonf) %>% bind_cols(., state_A = rep("normal", nrow(.)))

#get state classifications for allele B
gain_table_B <- possible_data2 %>% filter(chr != "chr10b") %>% filter(hapB_pval_CN1 < alpha_bonf) %>% filter(VAR_B > 1) %>% bind_cols(., state_B = rep("gain", nrow(.)))
gain_table_B_10b <- possible_data2 %>% filter(chr == "chr10b") %>% filter(hapB_pval_CN1 < alpha_bonf) %>% filter(VAR_B > 2) %>% bind_cols(., state_B = rep("gain", nrow(.)))
loss_table_B <- possible_data2 %>% filter(chr != "chr10b") %>% filter(hapB_pval_CN1 < alpha_bonf) %>% filter(VAR_B < 1) %>% bind_cols(., state_B = rep("loss", nrow(.)))
loss_table_B_10b <- possible_data2 %>% filter(chr == "chr10b") %>% filter(hapB_pval_CN1 < alpha_bonf) %>% filter(VAR_B < 2) %>% bind_cols(., state_B = rep("loss", nrow(.)))
ctrl_table_B <- possible_data2 %>% filter(hapB_pval_CN1 >= alpha_bonf) %>% bind_cols(., state_B = rep("normal", nrow(.)))

#bind all CN1 state classifcations together
ploidy_table_A <- rbind(gain_table_A, gain_table_A_X, loss_table_A, loss_table_A_X, ctrl_table_A)
ploidy_table_B <- rbind(gain_table_B, gain_table_B_10b, loss_table_B, loss_table_B_10b, ctrl_table_B)
ploidy_table <- merge(ploidy_table_A, ploidy_table_B, by = colnames(ploidy_table_A)[1:length(colnames(ploidy_table_A))-1])

#get CN2 state classifications for allele A and B
CN2_table_A <- possible_data2 %>% filter(hapA_pval_CN2 >= alpha_bonf) %>% bind_cols(., CN2_state_A = rep("full gain", nrow(.)))
CN2_table_A_intermediate <- possible_data2 %>% filter(hapA_pval_CN2 < alpha_bonf) %>% bind_cols(., CN2_state_A = rep("intermediate gain", nrow(.)))
CN2_table_B <- possible_data2 %>% filter(hapB_pval_CN2 >= alpha_bonf) %>% bind_cols(., CN2_state_B = rep("full gain", nrow(.)))
CN2_table_B_intermediate <- possible_data2 %>% filter(hapB_pval_CN2 < alpha_bonf) %>% bind_cols(., CN2_state_B = rep("intermediate gain", nrow(.)))

#bind all CN2 state classifications together
CN2_state_table_A <- rbind(CN2_table_A, CN2_table_A_intermediate)
CN2_state_table_B <- rbind(CN2_table_B, CN2_table_B_intermediate)
CN2_state_table <- merge(CN2_state_table_A, CN2_state_table_B, by = colnames(CN2_state_table_A)[1:length(colnames(CN2_state_table_A))-1])

#combine classifications
ploidy_table_2 <- merge(ploidy_table, CN2_state_table, by = intersect(colnames(ploidy_table), colnames(CN2_state_table)), all.x = T)
ploidy_table_2[is.na(ploidy_table_2)] <- "."

#combine state annotations
state_A2 <- rep(".", nrow(ploidy_table_2))
state_B2 <- rep(".", nrow(ploidy_table_2))
for(i in 1:nrow(ploidy_table_2)){
  #adjust states based on CIs and diploid distribution for allele A
  if(ploidy_table_2$state_A[i] == "loss" & ploidy_table_2$VAR_A[i] < full_loss_th){state_A2[i] <- "0"}
  if(ploidy_table_2$state_A[i] == "loss" & ploidy_table_2$VAR_A[i] > full_loss_th){state_A2[i] <- "1-"}
  if(ploidy_table_2$state_A[i] == "gain" & ploidy_table_2$CN2_state_A[i] == "intermediate gain"){state_A2[i] <- "1+"}
  if(ploidy_table_2$state_A[i] == "gain" & ploidy_table_2$CN2_state_A[i] == "full gain"){state_A2[i] <- "2"}
  if(ploidy_table_2$state_A[i] == "gain" & ploidy_table_2$CN2_state_A[i] == "."){state_A2[i] <- "1+"}
  if(ploidy_table_2$state_A[i] == "normal"){state_A2[i] <- "1"}
  #adjust states based on CIs and diploid distribution for allele B
  if(ploidy_table_2$state_B[i] == "loss" & ploidy_table_2$VAR_B[i] < full_loss_th){state_B2[i] <- "0"}
  if(ploidy_table_2$state_B[i] == "loss" & ploidy_table_2$VAR_B[i] > full_loss_th){state_B2[i] <- "1-"}
  if(ploidy_table_2$state_B[i] == "gain" & ploidy_table_2$CN2_state_B[i] == "intermediate gain"){state_B2[i] <- "1+"}
  if(ploidy_table_2$state_B[i] == "gain" & ploidy_table_2$CN2_state_B[i] == "full gain"){state_B2[i] <- "2"}
  if(ploidy_table_2$state_B[i] == "gain" & ploidy_table_2$CN2_state_B[i] == "."){state_B2[i] <- "1+"}
  if(ploidy_table_2$state_B[i] == "normal"){state_B2[i] <- "1"}
}

#reorder data
visual.data2 <- cbind(ploidy_table_2, state_A2, state_B2)
keycol <-c("generation", "family", "chr","relationship")
setorderv(visual.data2, keycol)

#Save visual.data2 storing all classification data
saveRDS(visual.data2, file = sprintf("%s/classification_data.rds", dirpath))


### cast the data into table ###
#Cast cell name values
seg_table_ids <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "ids")
setnames(seg_table_ids, old = c("A", "B", "c1", "c2"), new = c("A_id", "B_id", "c1_id", "c2_id"))
#Cast TPM values
seg_table_TPM <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "TPM")[,c("A", "B", "c1", "c2")]
setnames(seg_table_TPM, old = c("A", "B", "c1", "c2"), new = c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM"))
#cast allelic values
seg_table_VAR_A <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_A")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_A) <- c("A_AlleleA", "B_AlleleA", "c1_AlleleA", "c2_AlleleA") 
seg_table_VAR_B <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_B")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_B) <- c("A_AlleleB", "B_AlleleB", "c1_AlleleB", "c2_AlleleB") 
#cast CN state
seg_table_CN_A <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "state_A2")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_CN_A) <- c("A_AlleleA_state", "B_AlleleA_state", "c1_AlleleA_state", "c2_AlleleA_state") 
seg_table_CN_B <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "state_B2")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_CN_B) <- c("A_AlleleB_state", "B_AlleleB_state", "c1_AlleleB_state", "c2_AlleleB_state") 
#cast CN1 pval
seg_table_CN1_A <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "hapA_pval_CN1")[,c("A", "B", "c1", "c2")] 
colnames(seg_table_CN1_A) <- c("A_CN1_AlleleA", "B_CN1_AlleleA", "c1_CN1_AlleleA", "c2_CN1_AlleleA") 
seg_table_CN1_B <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "hapB_pval_CN1")[,c("A", "B", "c1", "c2")] 
colnames(seg_table_CN1_B) <- c("A_CN1_AlleleB", "B_CN1_AlleleB", "c1_CN1_AlleleB", "c2_CN1_AlleleB") 
#cast CN2 pval
seg_table_CN2_A <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "hapA_pval_CN2")[,c("A", "B", "c1", "c2")] 
colnames(seg_table_CN2_A) <- c("A_CN2_AlleleA", "B_CN2_AlleleA", "c1_CN2_AlleleA", "c2_CN2_AlleleA") 
seg_table_CN2_B <- dcast(visual.data2, generation + family + family_ids + MN_info + event + rupt_time + chr ~ relationship, value.var = "hapB_pval_CN2")[,c("A", "B", "c1", "c2")] 
colnames(seg_table_CN2_B) <- c("A_CN2_AlleleB", "B_CN2_AlleleB", "c1_CN2_AlleleB", "c2_CN2_AlleleB") 
#Bring all casted dfs together
seg_table <- cbind(seg_table_ids, seg_table_TPM, seg_table_VAR_A, seg_table_VAR_B, seg_table_CN_A, seg_table_CN_B, seg_table_CN1_A, seg_table_CN1_B, seg_table_CN2_A, seg_table_CN2_B)
#save
write.csv(seg_table, file = sprintf("%s/all_chr_profiles.csv", dirpath), row.names = F)

################################
### Fix notation for sorting ###
################################

all_family_events_2 <- data.table(seg_table)

#reorder columns for easier analysis
all_family_events_2$chr[all_family_events_2$chr == "chrX"] <- "chr23";
all_family_events_2$chr[all_family_events_2$chr == "chr10a"] <- "chr10.1"; all_family_events_2$chr[all_family_events_2$chr == "chr10b"] <- "chr10.2";
all_family_events_2$chr <- as.numeric(str_replace(string = all_family_events_2$chr, pattern = "chr", replacement = ""))
setkey(all_family_events_2, "generation", "family", "chr")

######################################
### Estimate CN pattern annotation ###
######################################

CN_values <- round(all_family_events_2[,c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM")])
CN_pattern_column <- paste(CN_values$A_TPM, CN_values$B_TPM, CN_values$c1_TPM, CN_values$c2_TPM)
CN_pattern_column_2 <- sub(sub(sub(x = CN_pattern_column, pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = "")

all_family_events_2$CN_pattern <- CN_pattern_column_2

###############################################
### Aggregate state classifications for fam ###
###############################################

agg_states_A_anno_col <- c()
agg_states_B_anno_col <- c()
for(i in 1:nrow(all_family_events_2)){
  agg_states_A <- unlist(unname(all_family_events_2[i,c("A_AlleleA_state", "B_AlleleA_state", "c1_AlleleA_state", "c2_AlleleA_state")]))
  agg_states_B <- unlist(unname(all_family_events_2[i,c("A_AlleleB_state", "B_AlleleB_state", "c1_AlleleB_state", "c2_AlleleB_state")]))
  agg_states_A_anno_col <- append(agg_states_A_anno_col, paste(agg_states_A[!is.na(agg_states_A)],  collapse = " "))
  agg_states_B_anno_col <- append(agg_states_B_anno_col, paste(agg_states_B[!is.na(agg_states_B)],  collapse = " "))
}

all_family_events_2$transcription_pattern_A <- as.character(agg_states_A_anno_col)
all_family_events_2$transcription_pattern_B <- as.character(agg_states_B_anno_col)

#############################
### Annotate MN haplotype ###
#############################

#if else representation of MN hap reasoning
hap_column <- c()
for(i in 1:nrow(all_family_events_2)) {
  tmp_hap <- c()
  tmp_event <- all_family_events_2[i,]
  if(length(setdiff(unlist(strsplit(tmp_event$transcription_pattern_A, split =  " ")), c("1"))) > 0){tmp_hap <- append(tmp_hap, "A")}
  if(length(setdiff(unlist(strsplit(tmp_event$transcription_pattern_B, split =  " ")), c("1"))) > 0){tmp_hap <- append(tmp_hap, "B")}
  if(is_empty(tmp_hap)){tmp_hap <- "normal"}
  if(length(tmp_hap) > 1){tmp_hap <- "NA"; print("Could not determine MN haplotpye")}
  hap_column <- append(hap_column, tmp_hap)
}

#add hap annotations to df
all_family_events_2$MN_hap <- hap_column

###############################
### Hap specific CN pattern ###
###############################

hap_spec_CN_pattern_column <- data.table()
for(i in 1:nrow(all_family_events_2)){
  if(all_family_events_2$MN_hap[i] %in% c("A", "B")){
    haps <- c(sprintf("A_Allele%s", all_family_events_2$MN_hap[i]), sprintf("B_Allele%s", all_family_events_2$MN_hap[i]), sprintf("c1_Allele%s", all_family_events_2$MN_hap[i]), sprintf("c2_Allele%s", all_family_events_2$MN_hap[i]))
    hap_spec_CN <- round(all_family_events_2[i,..haps])
    colnames(hap_spec_CN) <- c("A", "B", "c1", "c2")
    hap_spec_CN_pattern_column <- rbind(hap_spec_CN_pattern_column, hap_spec_CN)
  } else {hap_spec_CN_pattern_column <-  rbind(hap_spec_CN_pattern_column, data.table(NA, NA, NA, NA), use.names=FALSE)}
}
hap_spec_CN_pattern <- paste(hap_spec_CN_pattern_column$A, hap_spec_CN_pattern_column$B, hap_spec_CN_pattern_column$c1, hap_spec_CN_pattern_column$c2)
hap_spec_CN_pattern_2 <- sub(sub(sub(sub(x = hap_spec_CN_pattern, pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = "")

all_family_events_2$MNhap_CN_pattern <- hap_spec_CN_pattern_2

##########################
### AS_TPM for MN cell ###
##########################

hap_spec_CN_column <- c()
hap_state_spec_CN_column <- c()
hap_CN1_spec_CN_column <- c()
hap_CN2_spec_CN_column <- c()
for(i in 1:nrow(all_family_events_2)){
  if(all_family_events_2$MN_hap[i] %in% c("A", "B")){
    hap <- c(sprintf("A_Allele%s", all_family_events_2$MN_hap[i]))
    hap_state <- c(sprintf("A_Allele%s_state", all_family_events_2$MN_hap[i]))
    hap_pval_CN1 <- c(sprintf("A_CN1_Allele%s", all_family_events_2$MN_hap[i]))
    hap_pval_CN2 <- c(sprintf("A_CN2_Allele%s", all_family_events_2$MN_hap[i]))
    
    hap_spec_CN_column <- append(hap_spec_CN_column, all_family_events_2[i,..hap])
    hap_state_spec_CN_column <- append(hap_state_spec_CN_column, all_family_events_2[i,..hap_state])
    hap_CN1_spec_CN_column <- append(hap_CN1_spec_CN_column, all_family_events_2[i,..hap_pval_CN1])
    hap_CN2_spec_CN_column <- append(hap_CN2_spec_CN_column, all_family_events_2[i,..hap_pval_CN2])
    
  } else {
    hap_spec_CN_column <- append(hap_spec_CN_column, "NA")
    hap_state_spec_CN_column <- append(hap_state_spec_CN_column, "NA")
    hap_CN1_spec_CN_column <- append(hap_CN1_spec_CN_column, "NA")
    hap_CN2_spec_CN_column <- append(hap_CN2_spec_CN_column, "NA")
    
  }
}

all_family_events_2$MNhap_TPM <- as.character(hap_spec_CN_column)
all_family_events_2$MNhap_state <- as.character(hap_state_spec_CN_column)
all_family_events_2$MNhap_CN1_pval <- as.character(hap_CN1_spec_CN_column)
all_family_events_2$MNhap_CN2_pval <- as.character(hap_CN2_spec_CN_column)

#########################################
### Remove diploid/normal chromosomes ###
#########################################

#remove cells that are completely diploid
all_family_events_3 <- all_family_events_2 %>% filter(MN_hap != "normal") %>% data.table()
seg_table_diploid_events <- all_family_events_2 %>% filter(MN_hap == "normal")
write.csv(all_family_events_2, file = sprintf("%s/all_events.csv", dirpath), row.names = F)
write.csv(seg_table_diploid_events, file = sprintf("%s/diploid_events.csv", dirpath), row.names = F)

#Check if there are any completely diploid famliies
diploid_families <- all_family_events_2 %>% 
  group_by(family) %>%
  mutate(dummy = as.integer(n_distinct(MN_hap) == 1)) %>%
  filter(dummy == 1)

##########################################
### Reorder, round and save all events ###
##########################################

#Reorder cols
ordered_cols <- c("generation","MN_info","event","rupt_time","chr","family","family_ids","A_id","A_TPM","A_AlleleA","A_AlleleA_state","A_AlleleB","A_AlleleB_state","B_id","B_TPM","B_AlleleA","B_AlleleA_state","B_AlleleB","B_AlleleB_state","c1_id","c1_TPM","c1_AlleleA","c1_AlleleA_state","c1_AlleleB","c1_AlleleB_state","c2_id","c2_TPM","c2_AlleleA","c2_AlleleA_state","c2_AlleleB","c2_AlleleB_state","transcription_pattern_A","transcription_pattern_B", "MN_hap", "MNhap_TPM", "MNhap_state", "MNhap_CN1_pval", "MNhap_CN2_pval")
all_family_events_3 <- all_family_events_3[,..ordered_cols]
setkey(all_family_events_3, "generation", "family", "chr")

#save event table
write.csv(all_family_events_3, file = sprintf("%s/all_family_events.csv", dirpath), row.names = F)

#arm values for manual entry
AS_TPM_byarm <- readRDS(file=sprintf("%s/AS-TPM.inv_var.byarm.rds", dirpath))
AS_TPM_byarm_ctrl_pvals_arm_spec <- readRDS(file = sprintf("%s/AS_TPM_pval_matrix_control_byarm_armspec.rds", dirpath))



#summary visuals generator

##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#Here we make figures that summarize our interpretations from manual-review of the data
#This script was written by Nikos Mynhier 

####################
### loading data ###
####################

message("Loading data...")

#Read in results of manual review
all_family_events_results <- data.table(read.csv( sprintf("%s/manual_review_results.csv", dirpath)))
all_family_events_results$chr[all_family_events_results$chr == 'X'] <- 23
all_family_events_results$chr <- as.numeric(all_family_events_results$chr)

#load cell annotations
anno <- data.table(read.csv( sprintf("%s/annotation_list.csv", dirpath)))
changeCols <- colnames(anno)[-c(13:16)]
anno[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols] 
anno <- anno[which(anno$exclusion == ""),]

#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#Get MN cells and chrs
MN_Cell_n_chr_1 <- all_family_events_results[which(all_family_events_results$Interpretation %in% c("defect", "partial defect", "no defect")),c("family", "chr", "Interpretation", "generation")]
MN_Cell_n_chr_1$chr[-which(MN_Cell_n_chr_1$chr %in% c("10.1", "10.2"))] <- round(MN_Cell_n_chr_1$chr[-which(MN_Cell_n_chr_1$chr %in% c("10.1", "10.2"))])
MN_Cell_n_chr_1$chr <- as.character(MN_Cell_n_chr_1$chr)
MN_Cell_n_chr_2 <- data.table(anno_gen_1_2[which(anno_gen_1_2$Sister1 == 1),c("WTA.plate", "chr_of_interest", "Family_IDs")])
MN_Cell_n_chr_2$chr_of_interest <- str_remove(MN_Cell_n_chr_2$chr_of_interest, "chr")
MN_Cell_n_chr_2$chr_of_interest[MN_Cell_n_chr_2$chr_of_interest == "X"] <- 23
MN_Cell_n_chr_2$chr_of_interest <- as.character(MN_Cell_n_chr_2$chr_of_interest)
setnames(MN_Cell_n_chr_2, old = c("chr_of_interest", "Family_IDs"), new = c("chr", "family"))
setkey(MN_Cell_n_chr_1, "chr", "family")
setkey(MN_Cell_n_chr_2, "chr", "family")
MN_Cell_n_chr_2 <- MN_Cell_n_chr_2[!which(MN_Cell_n_chr_2$family %in% c("F73", "F79")),] #manually remove 2 cells that are too highly aneuploid
MN_Cell_n_chr_3 <- merge(MN_Cell_n_chr_1, MN_Cell_n_chr_2, all.x = T)
for(i in unique(MN_Cell_n_chr_3$family)){ #removes cases with no chromosome of interest
  tmp <- MN_Cell_n_chr_3[which(MN_Cell_n_chr_3$family == i),]
  tmp_cell <- tmp$WTA.plate[!is.na(tmp$WTA.plate)]
  if(!is_empty(tmp_cell)){MN_Cell_n_chr_3[which(MN_Cell_n_chr_3$family == i),]$WTA.plate <- tmp_cell} else {MN_Cell_n_chr_3 <- MN_Cell_n_chr_3[-which(MN_Cell_n_chr_3$family == i),]}
}
MN_Cell_n_chr <- copy(MN_Cell_n_chr_3)
MN_Cell_n_chr$WTA.plate[intersect(which(MN_Cell_n_chr$chr == 19), which(MN_Cell_n_chr$family == "F25"))] <- "170209_A8"
MN_Cells <- unique(MN_Cell_n_chr$WTA.plate)

#load classification data frame
visual.data.results.all <- readRDS(file = sprintf("%s/classification_data.rds", dirpath))

#Read in reference data
ref_TPM <- readRDS(file = sprintf("%s/ref_TPM.rds", dirpath))
ref_AS_TPM <- readRDS(file = sprintf("%s/ref_AS_TPM.rds", dirpath))
ref_hap_spec_CN <- readRDS(file = sprintf("%s/ref_hap_spec_CN.rds", dirpath))
ref_AS_TPM_byarm <- readRDS(file = sprintf("%s/ref_AS_TPM_byarm.rds", dirpath))

#new AS-TPM
AS_TPM_bychr <- readRDS(file=sprintf("%s/AS-TPM.inv_var.bychr.rds", dirpath))

#pvals for AS-TPM
AS_TPM_bychr_mono_pvals_chr_spec <- readRDS(file = sprintf("%s/pval_AS_TPM_bychr_mono_pvals_chr_spec.rds", dirpath))

##################################
### Hap Spec Curated ref distr ###
##################################

#main data to be plotted 
boxplots.vals <- cbind(ref_AS_TPM[,c("value", "CN", "variable")], color = ref_AS_TPM$CN)
setnames(boxplots.vals, old = c("value", "CN", "variable"), new = c("TPM", "Group", "allele"))
boxplots.vals$Group[boxplots.vals$Group == 0] <- "loss"; boxplots.vals$Group[boxplots.vals$Group == 1] <- "control"; boxplots.vals$Group[boxplots.vals$Group == 2] <- "gain";
boxplots.vals$allele[which(boxplots.vals$Group == "control")] <- "normal"
boxplots.vals$TPM <- as.numeric(boxplots.vals$TPM)
boxplots.vals$Group <- as.factor(boxplots.vals$Group)

#Plot the boxplots and strip plots
boxplot.visual.2 <- ggplot() + 
  #change theme
  theme_bw() + 
  theme(aspect.ratio = 1/3, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #MN experimental samples 
  #geom_boxplot(data = visual.data.results.all, mapping = aes(x = state, y = TPM, alpha = .5), outlier.alpha = 0)  +
  #selected samples
  geom_boxplot(data = boxplots.vals, mapping = aes(x = Group, y = TPM, alpha = .5, fill = allele), outlier.alpha = 0)  +
  geom_hline(yintercept = c(0,1,2), linetype="dotted") + 
  #all samples
  ylim(c(0,2)) + scale_x_discrete(limits = c("loss", "control", "gain"), labels=c( "ref monosomy", "control", "ref trisomy")) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

ggsave(plot = boxplot.visual.2, filename = "ref_distr_AS_TPM.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 15, height = 4, dpi = 300, units = "in")


####################################
### Plot TPM distributions bychr ###
####################################

#subset the AS ref distr
ctrl_AS_TPM_AB_melted <- ref_hap_spec_CN[which(ref_hap_spec_CN$CN == 1)]

#melt ctrl AS_TPM values for plotting
#ctrl_AS_TPM_AB_melted <- melt(ctrl_AS_TPM_AB, measure.vars = c("A", "B"))
ctrl_AS_TPM_AB_melted$vals[which(ctrl_AS_TPM_AB_melted$chr == "10b")] <- ctrl_AS_TPM_AB_melted$vals[which(ctrl_AS_TPM_AB_melted$chr == "10b")]*1.5
x_scale_f <- (mean(ctrl_AS_TPM_AB_melted[intersect(which(ctrl_AS_TPM_AB_melted$chr == "23"), which(ctrl_AS_TPM_AB_melted$hap == "B")),]$vals))
ctrl_AS_TPM_AB_melted$vals[which(ctrl_AS_TPM_AB_melted$chr == "23")] <- ctrl_AS_TPM_AB_melted$vals[which(ctrl_AS_TPM_AB_melted$chr == "23")]/(x_scale_f)

#subset the combined TPM ref distr
ctrl_TPM <- ref_TPM[which(ref_TPM$CN == 2)]

#curate ctrl TPM values for plotting
ctrl_TPM_melted <- cbind(ctrl_TPM, hap = c("A&B"))
ctrl_TPM_melted <- ctrl_TPM_melted[,c("ID", "CN", "chr", "hap", "vals")]
ctrl_TPM_melted$vals <- ctrl_TPM_melted$vals*2

#combine ref distributions
ref_all_chrs <- rbind(ctrl_AS_TPM_AB_melted, ctrl_TPM_melted)
ref_all_chrs$chr[ref_all_chrs$chr == "10a"] <- 10.1; ref_all_chrs$chr[ref_all_chrs$chr == "10b"] <- 10.2; 

#Get TPM values for each MN cell based on anno list
MN_Cells_TPM_A <- cbind(AS_TPM_bychr$A[,1], hap = "A", AS_TPM_bychr$A[,..MN_Cells])
MN_Cells_TPM_B <- cbind(AS_TPM_bychr$A[,1], hap = "B", AS_TPM_bychr$B[,..MN_Cells])
MN_Cells_TPM_A[which(MN_Cells_TPM_A$chr == "chr10b"),] <- c("chr10b", "A", MN_Cells_TPM_A[which(MN_Cells_TPM_A$chr == "chr10b"),-c(1:2)]*1.5)
MN_Cells_TPM_B[which(MN_Cells_TPM_B$chr == "chr10b"),] <- c("chr10b", "B", MN_Cells_TPM_B[which(MN_Cells_TPM_B$chr == "chr10b"),-c(1:2)]*1.5)
MN_Cells_TPM_A[which(MN_Cells_TPM_A$chr == "chrX"),] <- c("chrX", "A", MN_Cells_TPM_A[which(MN_Cells_TPM_A$chr == "chrX"),-c(1:2)]/x_scale_f)
MN_Cells_TPM_B[which(MN_Cells_TPM_B$chr == "chrX"),] <- c("chrX", "B", MN_Cells_TPM_B[which(MN_Cells_TPM_B$chr == "chrX"),-c(1:2)]/x_scale_f)
MN_Cells_TPM <- rbind(MN_Cells_TPM_A, MN_Cells_TPM_B)

#melt MN cells TPM so that it can be added to plot
melted_MN_TPM <- data.table(melt(MN_Cells_TPM, measure.vars = c(colnames(MN_Cells_TPM)[-c(1:2)])))
melted_MN_TPM$chr <- str_remove(melted_MN_TPM$chr, "chr")
melted_MN_TPM$chr[melted_MN_TPM$chr == "X"] <- "23"
melted_MN_TPM$chr[melted_MN_TPM$chr == "10a"] <- 10.1; melted_MN_TPM$chr[melted_MN_TPM$chr == "10b"] <- 10.2; 

#make seperate df with only chrs of interest in MN cell 
setnames(MN_Cell_n_chr, new=c("variable"), old=c("WTA.plate"))
setkey(MN_Cell_n_chr, "variable", "chr"); setkey(melted_MN_TPM, "variable", "chr"); 
melted_MN_TPM_chr_spec <- merge(MN_Cell_n_chr, melted_MN_TPM) #NOTE: 6 cases removed here for not having a chr of interest. 61->55

#merge with anno to get family info
setkey(melted_MN_TPM_chr_spec, "variable", "chr", "family")
anno_m <- copy(anno_gen_1_2)[,c("WTA.plate", "chr_of_interest", "Pairs")]
colnames(anno_m) <- c("variable", "chr", "family"); anno_m$chr[anno_m$chr == "chrX"] <- "chr23"; anno_m$chr <- str_remove(anno_m$chr, "chr");
setkey(anno_m, "variable", "chr", "family")
melted_MN_TPM_chr_spec_2 <- merge(melted_MN_TPM_chr_spec, anno_m, all.x = T)
melted_MN_TPM_chr_spec_2$chr <- as.character(as.numeric(melted_MN_TPM_chr_spec_2$chr))

#get MN hap from all family events
haplotype_info <- all_family_events_results[,c("family", "chr", "MN_hap", "generation", "MNhap_DNA_CN_pattern", "Interpretation", "rupt_time")]
haplotype_info$chr[is.na(haplotype_info$chr)] <- 1 #If we can determine is the chr of interest we record none
haplotype_info$chr[-which(haplotype_info$chr %in% c(10.1, 10.2))] <- round(haplotype_info$chr[-which(haplotype_info$chr %in% c(10.1, 10.2))]) #we round the chrs' # (not 10) to 1 digit
haplotype_info$chr <- as.character(haplotype_info$chr)
setkey(haplotype_info, "family", "chr", "generation", "Interpretation")
setkey(melted_MN_TPM_chr_spec_2, "family", "chr", "generation", "Interpretation")
melted_MN_TPM_chr_spec_3 <- merge(melted_MN_TPM_chr_spec_2, haplotype_info)
melted_MN_TPM_chr_spec_4 <- melted_MN_TPM_chr_spec_3[which(melted_MN_TPM_chr_spec_3$hap == melted_MN_TPM_chr_spec_3$MN_hap),]

#Add base CN state annotation
DNA_CN_patterns_labels <- data.table(t(data.frame(c("1 0", "Gen2", 1),c("1 0 1", "Gen2", 1),c("1 0 1 1", "Gen2", 1),c("1 1", "Gen1", 1),c("1 1", "Gen2", 2),c("1 1 1 1 | 1 0 1 1", "Gen2", "arm"),c("2 0", "Gen1", 2),c("2 0 | 1 1", "Gen2", "arm"),c("2 1 0", "Gen2", 2),c("2 1 0 0", "Gen2", 2),c("2 1 0 0 | 1 1 0 1", "Gen2", "arm"),c("2 2", "Gen1", 2),c("2 2 2 2 | 1 0 1 1", "Gen2", "arm"),c("2 1", "Gen2", 2),c("2 1 0 0 | 1 1 1 0", "Gen2", "arm"),c("2 1 | 1 1", "Gen2", "arm"),c("2 1 2 2", "Gen2", 2), c("1 1 | 2 1", "Gen2", "arm"), c("1 2", "Gen2", 2), c("3 1", "Gen1", 3), c("1 1 | 1 1", "Gen1", "arm"))))
colnames(DNA_CN_patterns_labels) <- c("MNhap_DNA_CN_pattern", "generation", "DNA_CN_base");
setkey(DNA_CN_patterns_labels, "MNhap_DNA_CN_pattern", "generation")
setkey(melted_MN_TPM_chr_spec_4, "MNhap_DNA_CN_pattern", "generation")
melted_MN_TPM_chr_spec_5 <- merge(melted_MN_TPM_chr_spec_4, DNA_CN_patterns_labels, all.x = T)

#Add a column to melted_MN_TPM_w_fam_2 which shows distance from mean TPM
melted_MN_TPM_grouped <- data.table( melted_MN_TPM %>% group_by(chr, hap) %>% mutate(avg_val = mean(value)) %>% mutate(sd_val = sd(value)) %>% mutate(zscore=abs(value-avg_val)/sd_val) ) 

#confidence intervals for ref data
control_CI <- ref_all_chrs[which(ref_all_chrs$hap == "A&B"),] %>% group_by(chr) %>% mutate(lower = quantile(x = vals, .05)) %>% mutate(upper = quantile(x = vals, .95)) %>% mutate(median = median(vals))
control_CI_2 <- data.table(control_CI[!duplicated(control_CI$upper),])[,c("chr", "hap", "upper", "lower", "median")]
hap_spec_CI <- ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),] %>% group_by(chr, hap) %>% mutate(lower = quantile(x = vals, .05)) %>% mutate(upper = quantile(x = vals, .95)) %>% mutate(median = median(vals))
hap_spec_CI_2 <- data.table(hap_spec_CI[!duplicated(hap_spec_CI$upper),])[,c("chr", "hap", "upper", "lower", "median")]

#merge upper and lower bounds with melted_MN_TPM_grouped
setkey(hap_spec_CI_2, "chr", "hap")
setkey(melted_MN_TPM_grouped, "chr", "hap")
melted_MN_TPM_grouped <- merge(melted_MN_TPM_grouped, hap_spec_CI_2)

#add an annotation for greater than 99% CI
melted_MN_TPM_grouped$significance <- NA
melted_MN_TPM_grouped$pval <- 2*pnorm(-abs(melted_MN_TPM_grouped$zscore))
sig_rows <- which(melted_MN_TPM_grouped$pval < (.05/48))
melted_MN_TPM_grouped$significance[sig_rows] <- "significant"

#Subset version of melted_MN_TPM_grouped only with significant rows
melted_MN_TPM_grouped_sig <- melted_MN_TPM_grouped[which(melted_MN_TPM_grouped$significance == "significant"),] #2.576 corresponds to 99% CI

#Make the reference distributions
make_path <- sprintf("mkdir -p %s/../visual_results/summary_visuals", dirpath)
system(make_path)

# Version 1 # 
# ignore all MN chrs only grey and box plot. Annotate the outliers

all_MN_chrs <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(data = hap_spec_CI_2, mapping = aes(x = chr, ymin = lower, ymax = upper, group = interaction(chr, hap)), width=.4, position=position_dodge(.9))  +
  #geom_errorbar(data = control_CI_2, mapping = aes(x = chr, ymin = lower, ymax = upper), width=.2)  +
  geom_point(data = melted_MN_TPM_grouped, mapping = aes(x = chr, y = value, group = interaction(chr, hap), fill = hap, alpha = zscore, color = significance), position=position_jitterdodge(jitter.width=.1)) + #position=position_jitterdodge()
  scale_x_discrete(limits = as.character(c(1:9, 10.1, 10.2, 11:23)), labels=as.character(c(1:9, "10a", "10b", 11:23))) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

ggsave(plot = all_MN_chrs, filename = "all_MN_chrs_distr.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 15, height = 10, dpi = 300, units = "in")

#outlier_chrs <- ggplot() + 
#  theme_bw() + 
#  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#  geom_errorbar(data = hap_spec_CI_2, mapping = aes(x = chr, ymin = lower, ymax = upper, group = interaction(chr, hap)), width=.4, position=position_dodge(.9))  +
#  geom_errorbar(data = control_CI_2, mapping = aes(x = chr, ymin = lower, ymax = upper), width=.2)  +
#  geom_point(data = melted_MN_TPM_grouped_sig, mapping = aes(x = chr, y = value, group = interaction(chr, hap), fill = hap, color = significance), position=position_jitterdodge()) + #position=position_jitterdodge()
#  scale_x_discrete(limits = as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_grouped_sig$chr)), labels=as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_grouped_sig$chr) )) +
#  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

# Version 2 # 
#Just show MN chrs on boxplots

#chr_intr_distr <- ggplot() + 
#  theme_bw() + 
#  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#  geom_boxplot(data = ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals, group = interaction(chr, hap)), outlier.alpha = 0)  +
#  geom_boxplot(data = ref_all_chrs[which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals), outlier.alpha = 0)  +
#  geom_jitter(data = melted_MN_TPM_chr_spec_5, mapping = aes(x = as.character(chr), y = value, group = interaction(chr, hap), color = generation, fill=factor(ifelse((melted_MN_TPM_chr_spec_5$Interpretation == "no defect"), NA, generation)), shape = DNA_CN_base), position=position_jitterdodge(), size = 2.5) + 
#  scale_x_discrete(limits = as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_chr_spec_5$chr)), labels= as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_chr_spec_5$chr))) +
#  scale_color_manual(values = c("Gen1" = "blue", "Gen2" = "red")) + 
#  scale_shape_manual(values=c(21,24,22)) +
#  scale_fill_manual(values = c("Gen1" = "blue", "Gen2" = "red"), na.value=NA, guide="none") +
#  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))


####################################################################################################
#Plots for David#

#First plot: Just control distrs
chr_spec_distrs <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1/6, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  geom_boxplot(data = ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals, group = interaction(chr, hap)), outlier.alpha = 0)  +
  geom_boxplot(data = ref_all_chrs[which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals), outlier.alpha = 0)  +
  scale_x_discrete(limits = as.character(c(1:9, 10.1, 10.2, 11:23)), labels=as.character(c(1:9, "10a", "10b", 11:23))) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("Control Distributions"))

ggsave(plot = chr_spec_distrs, filename = "chr_spec_ref_distr.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 15, height = 10, dpi = 300, units = "in")

#Second plot: Gen1 results over ctrl distrs
melted_MN_TPM_chr_spec_5_gen1 <- melted_MN_TPM_chr_spec_5 %>% filter(generation == "Gen1", DNA_CN_base %in% c("1", "2"))

plot_g1_summary <- function(rt) {
  ggplot() + 
    theme_bw() + 
    theme(aspect.ratio = 1/length(unique(melted_MN_TPM_chr_spec_5_gen1[which(melted_MN_TPM_chr_spec_5_gen1$rupt_time == rt)]$chr)), 
          text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
    geom_crossbar(data = hap_spec_CI_2, mapping = aes(x = chr, y = median, ymin = lower, ymax = upper, group = interaction(chr, hap)), width=.35, position=position_dodge(.6), fill = "grey90", color="grey60")  +
    geom_crossbar(data = control_CI_2, mapping = aes(x = chr, y = median, ymin = lower, ymax = upper), width=.5, fill = "grey90", color="grey60")  +
    geom_jitter(data = melted_MN_TPM_chr_spec_5_gen1[which(melted_MN_TPM_chr_spec_5_gen1$rupt_time == rt)], mapping = aes(x = as.character(chr), y = value, group = interaction(chr, hap), fill=Interpretation, shape = DNA_CN_base), position=position_jitterdodge(jitter.width=.2), size = 2.5) + 
    scale_x_discrete(limits = as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_chr_spec_5_gen1$chr[which(melted_MN_TPM_chr_spec_5_gen1$rupt_time == rt)]))) +  
    scale_y_continuous(limits = c(0,2.5)) +
    scale_shape_manual(values=c(21,24,22)) +
    scale_fill_manual(values = c("defect" = "blue", "no defect" = "red", "partial defect" = "green"), na.value=NA, guide="none") +
    labs(x = "Group", y = "TPM Ratio") 
}

plot_g1_1D_scatter <- function(rt) {
  ggplot() + 
    theme_bw() + 
    theme(aspect.ratio = 5,text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
    geom_point(data = melted_MN_TPM_chr_spec_5_gen1[which(melted_MN_TPM_chr_spec_5_gen1$rupt_time == rt)], mapping = aes(x = rt, y = value, fill=Interpretation, shape = DNA_CN_base), size = 2.5) + 
    scale_y_continuous(limits = c(0,2.5)) +
    scale_shape_manual(values=c(21,24,22)) +
    scale_fill_manual(values = c("defect" = "blue", "no defect" = "red", "partial defect" = "green"), na.value=NA, guide="none") +
    labs(x = "Group", y = "TPM Ratio", title = sprintf("Generation 1")) 
}

manual_order_rts <- c("intact", "ruptured")
gen1_summary <- lapply(manual_order_rts, plot_g1_summary)
gen1_summary_2 <- lapply(manual_order_rts, plot_g1_1D_scatter)

#Third plot: Gen2 results over ctrl distrs
melted_MN_TPM_chr_spec_5_gen2 <- melted_MN_TPM_chr_spec_5 %>% filter(generation == "Gen2", DNA_CN_base != "arm")

plot_g2_summary <- function(rt) {
  ggplot() + 
    theme_bw() + 
    theme(aspect.ratio = 1/length(unique(melted_MN_TPM_chr_spec_5_gen2[which(melted_MN_TPM_chr_spec_5_gen2$rupt_time == rt)]$chr)), 
          text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
    geom_crossbar(data = hap_spec_CI_2, mapping = aes(x = chr, y = median, ymin = lower, ymax = upper, group = interaction(chr, hap)), width=.35, position=position_dodge(.6), fill = "grey90", color="grey60")  +
    geom_crossbar(data = control_CI_2, mapping = aes(x = chr, y = median, ymin = lower, ymax = upper), width=.5, fill = "grey90", color="grey60")  +
    geom_jitter(data = melted_MN_TPM_chr_spec_5_gen2[which(melted_MN_TPM_chr_spec_5_gen2$rupt_time == rt)], mapping = aes(x = as.character(chr), y = value, group = interaction(chr, hap), fill=Interpretation, shape = DNA_CN_base), position=position_jitterdodge(jitter.width=.2), size = 2.5) + 
    scale_x_discrete(limits = as.character(intersect(c(1:9, 10.1, 10.2, 11:23), melted_MN_TPM_chr_spec_5_gen2$chr[which(melted_MN_TPM_chr_spec_5_gen2$rupt_time == rt)]))) +  
    scale_y_continuous(limits = c(0,2.5)) +
    scale_shape_manual(values=c(21,24,22)) +
    scale_fill_manual(values = c("defect" = "blue", "no defect" = "red"), na.value=NA, guide="none") +
    labs(x = "Group", y = "TPM Ratio") #, title = rt

}

plot_g2_1D_scatter <- function(rt) {
  ggplot() + 
    theme_bw() + 
    theme(aspect.ratio = 5,text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
    geom_point(data =  melted_MN_TPM_chr_spec_5_gen2[which(melted_MN_TPM_chr_spec_5_gen2$rupt_time == rt)], mapping = aes(x = rt, y = value, fill=Interpretation, shape = DNA_CN_base), size = 2.5) + 
    scale_y_continuous(limits = c(0,2.5)) +
    scale_shape_manual(values=c(21,24,22)) +
    scale_fill_manual(values = c("defect" = "blue", "no defect" = "red", "partial defect" = "green"), na.value=NA, guide="none") +
    labs(x = "Group", y = "TPM Ratio", title = "Generation 2") 
}

manual_order_rts <- c("intact", "ruptured")
gen2_summary <- lapply(manual_order_rts, plot_g2_summary)
gen2_summary_2 <- lapply(manual_order_rts, plot_g2_1D_scatter)

do.call("grid.arrange", c(list(chr_spec_distrs), gen1_summary, gen2_summary, nrow=5))
do.call("grid.arrange", c(list(chr_spec_distrs), gen1_summary_2, gen2_summary_2, ncol=5))

pdf(file = sprintf("%s/../visual_results/summary_visuals/analysis_results.pdf", dirpath), width = 11, height = 8.5)
do.call("grid.arrange", c(gen1_summary, gen2_summary, nrow=4))
dev.off()

pdf(file = sprintf("%s/../visual_results/summary_visuals/analysis_results_1D.pdf", dirpath), width = 11, height = 8.5)
do.call("grid.arrange", c(gen1_summary_2, gen2_summary_2, ncol=4))
dev.off()

g2_defect_DNA_CN <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1.5,text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") +
  geom_point(data =  melted_MN_TPM_chr_spec_5_gen2, mapping = aes(x = interaction(Interpretation, DNA_CN_base), y = value, fill=Interpretation, shape = DNA_CN_base), size = 2.5) + 
  #scale_x_discrete(limits =  as.character(unique(interaction(melted_MN_TPM_chr_spec_5_gen2$Interpretation, melted_MN_TPM_chr_spec_5_gen2$DNA_CN_base))), labels = c("No Defect 2:2", "Defect 2:2", "No Defect 3:1", "Defect 3:1")) +  
  scale_x_discrete(limits =  c( "defect.1", "no defect.1", "defect.2", "no defect.2"), labels = c("Defect 2:2", "No Defect 2:2", "Defect 3:1", "No Defect 3:1")) +  
  scale_y_continuous(limits = c(0,2.5)) +
  scale_shape_manual(values=c(21,24,22)) +
  scale_fill_manual(values = c("defect" = "blue", "no defect" = "red", "partial defect" = "green"), na.value=NA, guide="none") +
  labs(x = "Group", y = "TPM Ratio", title = "Generation 2") 

ggsave(plot = g2_defect_DNA_CN, filename = "g2_summary_alt.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 15, height = 10, dpi = 300, units = "in")


############################################
### Summarize Defect vs No Defect by chr ###
############################################

#Get generation 2 events including arm events
gen2_events <- melted_MN_TPM_chr_spec_5[which(melted_MN_TPM_chr_spec_5$generation == "Gen2"),]
gen2_events_2 <- gen2_events %>% filter(rupt_time == "intact", Interpretation == "no defect")

#Plot barchart of gen2 defect
gen2_barchart <- ggplot(gen2_events, aes(x = as.character(gen2_events$chr), fill = Interpretation)) + 
  scale_x_discrete(limits=as.character(c(1,2,4:9,10.2,11:13,19,23))) +
  scale_fill_grey() + theme_classic() +
  labs(x = "Chromosome", y = "Number of MN events") +
  geom_bar()

ggsave(plot = gen2_barchart, filename = "gen2_defect_frequency.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 8, height = 6)

all_cells_barchart <- ggplot(melted_MN_TPM_chr_spec_5, aes(x = as.character(melted_MN_TPM_chr_spec_5$chr))) + 
  scale_x_discrete(limits=as.character(c(1:9,10.1,10.2,11:13,15,19,22:23))) +
  scale_fill_grey() + theme_classic() +
  labs(x = "Chromosome", y = "Number of MN events") +
  geom_bar()

ggsave(plot = all_cells_barchart, filename = "MN_event_frequency.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf",  width = 8, height = 6)

###################################################
### Summarize Defect vs No Defect in each group ###
###################################################

events_summary <- all_family_events_results %>% 
  select(generation, rupt_time, Interpretation, family) %>% 
  filter(Interpretation %in% c("defect", "partial defect", "no defect"))

#Split by generation
events_summary_g1 <- events_summary %>% filter(generation == "Gen1")
events_summary_g2 <- events_summary %>% filter(generation == "Gen2")

#plot barchart for each generation
gen1_rupt_time <- ggplot(data = events_summary_g1, mapping = aes(x = rupt_time, fill = Interpretation)) + 
  geom_bar(position = position_stack(reverse = TRUE)) + 
  scale_x_discrete(limits=as.character(c("intact", "ruptured"))) +
  scale_fill_grey() + theme_classic() +
  labs(x = "Rupture Timing", y = "Number of MN chromosomes")

ggsave(plot = gen1_rupt_time, filename = "gen1_defect_frequency_bytime.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf")

gen2_rupt_time <- ggplot(data = events_summary_g2, mapping = aes(x = rupt_time, fill = Interpretation)) + 
  geom_bar() + 
  scale_x_discrete(limits=as.character(c("intact", "ruptured"))) +
  scale_fill_grey() + theme_classic() +
  labs(x = "Rupture Timing", y = "Number of MN chromosomes")

ggsave(plot = gen2_rupt_time, filename = "gen2_defect_frequency_bytime.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf")

#####################################################

#Pull all exp families
MN_families <- data.table(anno_gen_1_2[,c("WTA.plate", "chr_of_interest", "Family_IDs", "Pairs")])

#Get TPM values for each MN cell based on anno list
experimental_cells <- MN_families$WTA.plate
Exp_Cells_TPM_A <- cbind(AS_TPM_bychr$A[,1], hap = "A", AS_TPM_bychr$A[,..experimental_cells])
Exp_Cells_TPM_B <- cbind(AS_TPM_bychr$A[,1], hap = "B", AS_TPM_bychr$B[,..experimental_cells])
Exp_Cells_TPM_A[which(Exp_Cells_TPM_A$chr == "chr10b"),] <- c("chr10b", "A", Exp_Cells_TPM_A[which(Exp_Cells_TPM_A$chr == "chr10b"),-c(1:2)]*1.5)
Exp_Cells_TPM_B[which(Exp_Cells_TPM_B$chr == "chr10b"),] <- c("chr10b", "B", Exp_Cells_TPM_B[which(Exp_Cells_TPM_B$chr == "chr10b"),-c(1:2)]*1.5)
Exp_Cells_TPM_A[which(Exp_Cells_TPM_A$chr == "chrX"),] <- c("chrX", "A", Exp_Cells_TPM_A[which(Exp_Cells_TPM_A$chr == "chrX"),-c(1:2)]/x_scale_f)
Exp_Cells_TPM_B[which(Exp_Cells_TPM_B$chr == "chrX"),] <- c("chrX", "B", Exp_Cells_TPM_B[which(Exp_Cells_TPM_B$chr == "chrX"),-c(1:2)]/x_scale_f)
Exp_Cells_TPM <- rbind(Exp_Cells_TPM_A, Exp_Cells_TPM_B)

#melt Exp cells TPM so that it can be added to plot
melted_exp_cell_TPM <- data.table(melt(Exp_Cells_TPM, measure.vars = c(colnames(Exp_Cells_TPM)[-c(1:2)])))
melted_exp_cell_TPM$chr <- str_remove(melted_exp_cell_TPM$chr, "chr")
melted_exp_cell_TPM$chr[melted_exp_cell_TPM$chr == "X"] <- "23"
melted_exp_cell_TPM$chr[melted_exp_cell_TPM$chr == "10a"] <- 10.1; melted_exp_cell_TPM$chr[melted_exp_cell_TPM$chr == "10b"] <- 10.2
melted_exp_cell_TPM <- melted_exp_cell_TPM %>% setnames(old = c("variable", "value"), new = c("cell", "TPM"))

#Get TPM values for each MN cell based on anno list
experimental_cells <- MN_families$WTA.plate
Exp_Cells_pval_A <- cbind(AS_TPM_bychr_mono_pvals_chr_spec$A[,1], hap = "A", AS_TPM_bychr_mono_pvals_chr_spec$A[,..experimental_cells])
Exp_Cells_pval_B <- cbind(AS_TPM_bychr_mono_pvals_chr_spec$A[,1], hap = "B", AS_TPM_bychr_mono_pvals_chr_spec$B[,..experimental_cells])
Exp_Cells_pval <- rbind(Exp_Cells_pval_A, Exp_Cells_pval_B)

#melt Exp cells pvals so that it can be added to plot
melted_exp_cell_pval <- data.table(melt(Exp_Cells_pval, measure.vars = c(colnames(Exp_Cells_pval)[-c(1:2)])))
melted_exp_cell_pval$chr <- str_remove(melted_exp_cell_pval$chr, "chr")
melted_exp_cell_pval$chr[melted_exp_cell_pval$chr == "X"] <- "23"
melted_exp_cell_pval$chr[melted_exp_cell_pval$chr == "10a"] <- 10.1; melted_exp_cell_pval$chr[melted_exp_cell_pval$chr == "10b"] <- 10.2
melted_exp_cell_pval <- melted_exp_cell_pval %>% setnames(old = c("variable", "value"), new = c("cell", "pval"))

#calculate confidence intervals
hap_spec_CI <- ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),] %>% group_by(chr, hap) %>% mutate(lower = quantile(x = vals, .05)) %>% mutate(upper = quantile(x = vals, .95)) %>% mutate(median = median(vals))
hap_spec_CI_2 <- data.table(hap_spec_CI[!duplicated(hap_spec_CI$upper),])[,c("chr", "hap", "upper", "lower", "median")]

#merge upper and lower bounds with melted_MN_TPM_grouped
setkey(melted_exp_cell_TPM, "chr", "hap", "cell")
setkey(melted_exp_cell_pval, "chr", "hap", "cell")
melted_exp_cell_grouped <- merge(melted_exp_cell_TPM, melted_exp_cell_pval)

#add an annotation for greater than 99% CI
melted_exp_cell_grouped$significance <- NA
sig_rows <- which(melted_exp_cell_grouped$pval < (.05/48))
melted_exp_cell_grouped$significance[sig_rows] <- "significant"

all_exp_chrs <- ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data = melted_exp_cell_grouped, mapping = aes(x = chr, y = TPM, group = interaction(chr, hap), fill = hap, alpha = .5, color = significance), position=position_jitterdodge(jitter.width=.1)) + #position=position_jitterdodge()
  geom_errorbar(data = hap_spec_CI_2, mapping = aes(x = chr, ymin = lower, ymax = upper, group = interaction(chr, hap)), width=.4, position=position_dodge(.9))  +
  scale_x_discrete(limits = as.character(c(1:9, 10.1, 10.2, 11:23)), labels=as.character(c(1:9, "10a", "10b", 11:23))) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

ggsave(plot = all_exp_chrs, filename = "all_exp_chrs_distr.pdf", path = sprintf("%s/../visual_results/summary_visuals", dirpath), device = "pdf", width = 15, height = 10, dpi = 300, units = "in")




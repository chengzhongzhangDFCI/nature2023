plot_AS_TPM_ratio_on_reference <- function(ref_hap_spec_CN, ref_TPM, AS_TPM_bychr, 
                                           destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results",
                                           myfamily_file = "F258", myid_files = c("F258.1", "F258.2", "F258.3", "F258.4"), myids = c("210503_9A", "210503_9B", "210503_9C", "210503_9D"), 
                                           MN_related = c("210503_9C", "210503_9D"), MN_unrelated = c("210503_9A", "210503_9B"), save.vis = F, return.vis = F) {
  
  #Conditional line for gen1 families with 2 cells
  gen1_families <- anno$Family_IDs[which(anno$key_pairs == "Gen1")]
  if(myfamily_file %in% gen1_families){MN_unrelated <- MN_related[2]; MN_related <- MN_related[1]}
  
  #######################################
  ### Visualize control distributions ###
  #######################################
  
  #subset the AS ref distr
  ctrl_AS_TPM_AB_melted <- ref_hap_spec_CN[which(ref_hap_spec_CN$CN == 1)]
  
  #melt ctrl AS_TPM values for plotting
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
  
  #plot the control distributions 
  chr_spec_distrs <- ggplot(data = ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals, group = interaction(chr, hap))) + 
    geom_boxplot(data = ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals, group = interaction(chr, hap)), outlier.alpha = 0)  +
    geom_boxplot(data = ref_all_chrs[which(ref_all_chrs$hap == "A&B"),], mapping = aes(x = chr, y = vals), outlier.alpha = 0)  +
    scale_colour_manual(values = c("gray40", "gray70")) +
    scale_x_discrete(limits = as.character(c(1:9, 10.1, 10.2, 11:23)), labels=as.character(c(1:9, "10a", "10b", 11:23))) +
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
    coord_cartesian(clip="off") +
    geom_rangeframe(x=0, y=seq(0, 2.5, along.with = ref_all_chrs[-which(ref_all_chrs$hap == "A&B"),]$vals)) + 
    theme_tufte() +
    theme(aspect.ratio = 1/3, text=element_text(size=18), axis.ticks.x = element_blank(), axis.ticks.length=unit(.25, "cm"), 
          axis.text = element_text(color="black"), axis.title = element_blank()) 

  ################################################################
  ### plot a single family on top of the control distributions ###
  ################################################################
  
  #Get TPM values for each MN cell based on anno list
  family_TPM_A <- cbind(AS_TPM_bychr$A[,1], hap = "A", AS_TPM_bychr$A[,..myids])
  family_TPM_B <- cbind(AS_TPM_bychr$A[,1], hap = "B", AS_TPM_bychr$B[,..myids])
  family_TPM_A[which(family_TPM_A$chr == "chr10b"),] <- c("chr10b", "A", family_TPM_A[which(family_TPM_A$chr == "chr10b"),-c(1:2)]*1.5)
  family_TPM_B[which(family_TPM_B$chr == "chr10b"),] <- c("chr10b", "B", family_TPM_B[which(family_TPM_B$chr == "chr10b"),-c(1:2)]*1.5)
  family_TPM_A[which(family_TPM_A$chr == "chrX"),] <- c("chrX", "A", family_TPM_A[which(family_TPM_A$chr == "chrX"),-c(1:2)]/x_scale_f)
  family_TPM_B[which(family_TPM_B$chr == "chrX"),] <- c("chrX", "B", family_TPM_B[which(family_TPM_B$chr == "chrX"),-c(1:2)]/x_scale_f)
  family_TPM <- rbind(family_TPM_A, family_TPM_B)
  
  #melt family cells TPM so that it can be added to plot
  melted_family_TPM <- data.table(melt(family_TPM, measure.vars = c(colnames(family_TPM)[-c(1:2)])))
  melted_family_TPM$chr <- str_remove(melted_family_TPM$chr, "chr")
  melted_family_TPM$chr[melted_family_TPM$chr == "X"] <- "23"
  melted_family_TPM$chr[melted_family_TPM$chr == "10a"] <- 10.1; melted_family_TPM$chr[melted_family_TPM$chr == "10b"] <- 10.2; 
  
  #add MN related information to melted TPM
  melted_family_TPM$MN_relation <- rep(NA, nrow(melted_family_TPM))
  melted_family_TPM$MN_relation[which(melted_family_TPM$variable %in% MN_related)] <- "MN Related"
  melted_family_TPM$MN_relation[which(melted_family_TPM$variable %in% MN_unrelated)] <- "MN Unrelated"
  
  #Add individual family
  #chr_spec_distrs_2 <- chr_spec_distrs +  geom_jitter(data = melted_family_TPM, mapping = aes(x = as.character(chr), y = value, group = interaction(chr, hap), color = MN_relation), position=position_jitterdodge(),  size = 2.5)
  chr_spec_distrs_1 <- chr_spec_distrs +  geom_point(data = melted_family_TPM, mapping = aes(x = as.character(chr), y = value, group = interaction(chr, hap), color = MN_relation),  position=position_dodge(width=0.75), size = 2.5)
  
  #save as pdf for print use
  if(save.vis){ggsave(plot = chr_spec_distrs_1, filename = sprintf("%s.family_summary.pdf", myfamily_file), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")}
  
  #print and return visuals if true
  if(return.vis){return(chr_spec_distrs_1)}
  
}
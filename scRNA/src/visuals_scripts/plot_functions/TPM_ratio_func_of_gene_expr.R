CumulativeTPM_and_ratioTPM_plots <- function(adt, geneRanges, controlSampleIDs,
                                             destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results",
                                             myfamily_file = "F258", myid_file = "F258.3", myid = "210503_9A", chr = "chr5", 
                                             binsize = 10, save.rds = F) {
  
  cols <- c("gene_id","chr","start","end","width","strand",myid)
  
  adt_all <- adt[,..cols]
  
  adt_control <- data.table(cbind(adt_all, control_TPM = rowMeans(adt[,..controlSampleIDs])))
  
  setkey(adt_control, "control_TPM")
  
  adt_control <- adt_control[which(adt_control$control_TPM > 0),]
  
  adt_control$ratio <- (adt_control[,..myid] / adt_control[,c("control_TPM")])
  
  rows <- which(adt_control$chr == chr)
  adt_control_chr <- adt_control[rows,]
  
  #########################################################
  ### Plot TPM ratio as function of avg gene expression ###
  #########################################################
  
  #bin genes
  num_to_rep <- ceiling(nrow(adt_control_chr)/binsize)
  bin <- rep(c(1:num_to_rep), each=binsize)[1:nrow(adt_control_chr)]
  adt_control_bin <- cbind(adt_control_chr, bin)
  
  adt_bin <- adt_control_bin[,-c(1:6)] %>%
    group_by(bin) %>%
    summarise_all(mean, na.rm = TRUE)
  
  adt_bin$control_TPM_logged <- log2(adt_bin$control_TPM)
  ticks <- as.character(c(1,seq(100,4096,100)))
  ticks[!ticks %in% c("1","100","1000","4000")] <- ""
  
  p1 <- ggplot(adt_bin, aes(x = control_TPM, y = ratio)) + 
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
    geom_point(color = "black", aes(x = control_TPM, y = ratio)) + 
    coord_trans(x="log2") +
    scale_x_continuous(breaks = c(1,seq(100,4096,100)), limits = c(1,4096), labels = ticks) +
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), limits = c(0,2.5), labels = c("0","","1","","2","")) +
    ggtitle(sprintf("ID: %s | Chr: %s", myid_file, chr)) +
    xlab("Avg Expression Level") +
    ylab("TPM Ratio to Normal") +
    geom_rangeframe(x=seq(1, 4096, along.with = adt_bin$control_TPM) , y=seq(0, 2.5, along.with = adt_bin$ratio)) + 
    theme_tufte() +
    #coord_cartesian(clip="off") +
    theme(aspect.ratio = 1/5, text=element_text(size=18), axis.text = element_text(color="black"), 
          axis.title = element_text(size=10), axis.ticks.length=unit(.25, "cm")) 
  
  ###############################################
  ### Plot cumulative TPM for sample vs ctrl ###
  ###############################################
  
  adt_control_chr$samp_cumsum <- cumsum(adt_control_chr[,..myid])
  adt_control_chr$ctrl_cumsum <- cumsum(adt_control_chr$control_TPM)
  
  hypot <- sqrt(max(adt_control_chr$ctrl_cumsum)^2 + max(adt_control_chr$samp_cumsum)^2)/2
  
  p2 <- ggplot(adt_control_chr, aes(x = ctrl_cumsum, y = samp_cumsum)) + 
    ylim(c(0, max(adt_control_chr$ctrl_cumsum))) + 
    xlim(c(0, max(adt_control_chr$ctrl_cumsum))) + 
    geom_point(color = "black", aes(x = ctrl_cumsum, y = samp_cumsum)) + 
    geom_segment(aes(x = 0, y = 0, xend = max(adt_control_chr$ctrl_cumsum), yend = max(adt_control_chr$ctrl_cumsum)),linetype = 2, size = 1) + #slope 1 line
    geom_segment(aes(x = 0, y = 0, xend = max(adt_control_chr$ctrl_cumsum), yend = max(adt_control_chr$ctrl_cumsum)/2),linetype=2, size = .5) + #slope .5 line
    geom_segment(aes(x = 0, y = 0, xend = max(adt_control_chr$ctrl_cumsum)/1.5, yend = max(adt_control_chr$ctrl_cumsum)),linetype=2, size = .5) + #slope 1.5 line
    theme_tufte() +
    ggtitle(sprintf("ID: %s | Chr: %s", myid_file, chr)) +
    xlab("Control Cumulative Sum") +
    ylab("Sample Cumulative Sum") +
    #coord_cartesian(clip="off") +
    geom_rangeframe(x = seq(0, max(adt_control_chr$ctrl_cumsum), along.with = adt_control_chr$ctrl_cumsum), y = seq(0, max(adt_control_chr$ctrl_cumsum), along.with = adt_control_chr$ctrl_cumsum) ) + 
    theme(aspect.ratio = 1/3, text=element_text(size=18), axis.text = element_text(color="black"), axis.title = element_text(size=10)) 

  #save raw visuals as pdf for print use
  ggsave(plot = p1, filename = sprintf("%s.%s.%s.TPMratiobyExpr.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
  ggsave(plot = p2, filename = sprintf("%s.%s.%s.CumulativeTPM.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 3, height = 3, dpi = 300, units = "in")
  
  if(save.rds == T) {
    #save raw visuals in r format for grid arrangement later
    saveRDS(p1, file = sprintf("%s/%s.%s.%s.TPMratiobyExpr.rds", destDir, myfamily_file, myid_file, chr))
    saveRDS(p2, file = sprintf("%s/%s.%s.%s.CumulativeTPM.rds", destDir, myfamily_file, myid_file, chr))
    return( data.table(expr_file = sprintf("%s/%s.%s.%s.TPMratiobyExpr.rds", destDir, myfamily_file, myid_file, chr), cum_file = sprintf("%s/%s.%s.%s.CumulativeTPM.rds", destDir, myfamily_file, myid_file, chr)) )
  }

}
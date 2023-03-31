plot_barplots_of_AllelicAndExpBiasPerSamples <- function(AS_TPM_bychr, single_sample = F, plotAxisText = T, return_plot = F, save_plot = F,
                                                         ids = c("210503_9A", "210503_9B", "210503_9C", "210503_9D"), chr = "chr12", 
                                                         name_ids = c("F258.1","F258.2","F258.3","F258.4"), myfamily_file = "F258",
                                                         family_relationship = c(),
                                                         destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results") {

  #get subset of the AS_TPM_bychr for specified cells
  cols <- c("chr", ids) #set columns to get from AS_TPM data.table
  cell_spec_AS_TPM_A <- melt(AS_TPM_bychr$A[,..cols], measure.vars = ids) #melt (reshape) the matrix into a table
  cell_spec_AS_TPM_B <- melt(AS_TPM_bychr$B[,..cols], measure.vars = ids) #melt (reshape) the matrix into a table
  cell_spec_AS_TPM <- rbind( cbind(cell_spec_AS_TPM_A, "hap" = "A"), cbind(cell_spec_AS_TPM_B, "hap" = "B")) #combine A and B alleles into one table
  setnames(cell_spec_AS_TPM, old = c("variable", "value"), new = c("cell", "expression"))
  cell_spec_AS_TPM_2 <- cell_spec_AS_TPM %>% group_by(chr, cell) %>% mutate(total = sum(expression))
  
  if(length(family_relationship) != 0){
    #Add relationship information
    relationlist <- data.table( cbind(ids, name_ids, family_relationship))
    cell_spec_AS_TPM_2 <- data.table(cell_spec_AS_TPM_2)
    setnames(relationlist, old=c("ids"), new=c("cell"))
    setkey(relationlist, cell)
    setkey(cell_spec_AS_TPM_2, cell)
    cell_spec_AS_TPM_2 <- merge(cell_spec_AS_TPM_2, relationlist)
    setkey(cell_spec_AS_TPM_2, family_relationship, cell, chr)
    setkey(relationlist, family_relationship, cell)
  }
  
  #use numeric version of chr
  str_replacement_1 <- str_replace(string = cell_spec_AS_TPM_2$chr, pattern = "chr", replacement = ""); str_replacement_2 <- str_replace(string = str_replacement_1, pattern = "a", replacement = ".1"); str_replacement_3 <- str_replace(string = str_replacement_2, pattern = "b", replacement = ".2"); str_replacement_4 <- str_replace(string = str_replacement_3, pattern = "X", replacement = "23");
  cell_spec_AS_TPM_2$chr <- as.numeric(str_replacement_4)
  
  #manually adjust chrX and chr10b 
  x_norm_factor <- rowMeans(AS_TPM_bychr$B[which(AS_TPM_bychr$B$chr == "chrX"),..controlSampleIDs])
  cell_spec_AS_TPM_2$expression[grep("23", cell_spec_AS_TPM_2$chr)] <- cell_spec_AS_TPM_2$expression[grep("23", cell_spec_AS_TPM_2$chr)]/x_norm_factor #manually divide X by 2
  cell_spec_AS_TPM_2$expression[grep("10.2", cell_spec_AS_TPM_2$chr)] <- cell_spec_AS_TPM_2$expression[grep("10.2", cell_spec_AS_TPM_2$chr)]*1.5 #manually divide X by 2
  
  if(length(family_relationship) != 0){
    #plot stacked barchart
    p <- ggplot(cell_spec_AS_TPM_2, aes(x = factor(family_relationship, levels = c("A", "B", "c1", "c2"), ordered = T), y = expression, fill = hap)) + 
      geom_bar(stat = "identity", position  = 'stack', width = 0.8, size = 0.7, colour = "black") +
      facet_wrap(~chr, nrow = 1, strip.position = "bottom") +
      ggtitle(sprintf("IDs: %s", paste(relationlist$name_ids, collapse = " | "))) +
      scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3), labels = c("0","","1","","2","","3"), limits = c(0, 4)) + #if boxs are getting clipped off remove or expand these limits
      scale_fill_manual(values = c("dodgerblue3", "firebrick2")) + 
      theme_tufte() +
      coord_cartesian(clip="off") +
      theme(aspect.ratio = 4, text=element_text(size=17), axis.ticks.x = element_blank(), # 
            axis.text.x = element_blank(), axis.title = element_blank(), panel.spacing = unit(0, "lines"))
  } else {
    #plot stacked barchart
    p <- ggplot(cell_spec_AS_TPM_2, aes(x = cell, y = expression, fill = hap)) + 
      geom_bar(stat = "identity", position  = 'stack', width = 0.8, size = 0.7, colour = "black") +
      facet_wrap(~chr, nrow = 1, strip.position = "bottom") +
      ggtitle(sprintf("IDs: %s", paste(name_ids, collapse = " | "))) +
      scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5,3), labels = c("0","","1","","2","","3"), limits = c(0, 4)) + #if boxs are getting clipped off remove or expand these limits
      scale_fill_manual(values = c("dodgerblue3", "firebrick2")) + 
      theme_tufte() +
      coord_cartesian(clip="off") +
      theme(aspect.ratio = 4, text=element_text(size=17), axis.ticks.x = element_blank(), # 
            axis.text.x = element_blank(), axis.title = element_blank(), panel.spacing = unit(0, "lines"))
  }
  
  
  #save the plots
  if(save_plot == T) {
    if(single_sample == T){
      ggsave(plot = p, filename = sprintf("%s.ExprRatio.pdf", name_ids), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
    } else {
      ggsave(plot = p, filename = sprintf("%s.ExprRatio.pdf", myfamily_file), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
    }
  }
  if(return_plot == T) {
    plot(p)
    return(p)
  }
  
}

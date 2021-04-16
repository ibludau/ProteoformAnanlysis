plotProteoformVolcanoLines <- function(traces, name="proteoform_pseudo_volcano_lines", score_cutoff = 0.1, adj_pval_cutoff =  0.1, PDF=TRUE){
  if (PDF){
    pdf(paste0(name,".pdf"), width=4, height=4)
  }
  scores <- unique(traces$trace_annotation[,.(proteoform_score,proteoform_score_pval,proteoform_score_pval_adj), by=.(protein_id)])
  scores <- scores[!is.na(scores$proteoform_score)]
  scores[, sig := ifelse((proteoform_score_pval_adj <= adj_pval_cutoff) & (proteoform_score >= score_cutoff), "sig", "non_sig"), by='protein_id']
  p <- ggplot(scores, aes(x=proteoform_score, y=-log10(proteoform_score_pval_adj), colour=sig)) +
    geom_point(alpha=0.5, show.legend = FALSE) +
    geom_vline(xintercept = score_cutoff, lty="dashed", color="darkgrey") +
    geom_hline(yintercept = -log10(adj_pval_cutoff), lty="dashed", color="darkgrey") +
    scale_colour_manual(values = c("grey", "green4")) +
    xlab("proteoform score") +
    ylab("-log10 (adj. p-value) ") +
    theme_classic()
  print(p)
  if (PDF){
    dev.off()
  }
}
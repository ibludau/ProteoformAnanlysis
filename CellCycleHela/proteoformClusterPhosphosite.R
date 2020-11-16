plotPeptideClusterSigPhosphosites <- function(traces,protein, PDF=FALSE, closeGaps=FALSE){
  traces <- subset(traces, protein, trace_subset_type="protein_id")
  dt <- subset(traces$trace_annotation,protein_id==protein)
  setkeyv(dt, c("protein_id","PeptidePositionStart"))
  dt[,PeptidePositionStartRank := seq_len(.N), by="protein_id"]
  cbPalette <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#756bb1","#1c9099")
  dt[,cluster:=ifelse(cluster==0, NA, cluster)]
  dt$cluster <- as.factor(dt$cluster)
  dt$phophoPeptide <- as.factor(dt$phophoPeptide)
  if (PDF){
    pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=3)
  }
  if (closeGaps) {
    q <- ggplot(dt,aes(x=PeptidePositionStartRank,
                       y=1,
                       fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) 
    p <- ggplot(dt,aes(x=PeptidePositionStartRank,
                       y=1,
                       fill=phophoPeptide)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) 
    f <- ggarrange(p, q, 
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
    print(annotate_figure(f, fig.lab = paste0(protein)))
    #print(f + ggtitle(paste0(protein," : ",unique(dt$Gene_names))))
  } else {
    q <- ggplot(dt,aes(x=PeptidePositionStart,
                       y=1,
                       fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette)
    p <- ggplot(dt,aes(x=PeptidePositionStart,
                       y=1,
                       fill=phophoPeptide)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette)
    f <- ggarrange(p, q, 
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
    print(annotate_figure(f, fig.lab = paste0(protein)))
    #print(f + ggtitle(paste0(protein," : ",unique(dt$Gene_names)))) 
  }
  if (PDF){
    dev.off()
  }
}

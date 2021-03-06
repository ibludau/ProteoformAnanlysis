#' ---
#' title: "Mouse tissue proteoform annotation "
#' author: "Isabell Bludau"
#' output:
#'   html_document:
#'     toc: true
#'   html_notebook:
#'     toc: true
#'   pdf_document:
#'     toc: true
#' ---

#rmarkdown::render("GetMouseTissueProteoforms.R", "pdf_document")

#' ## Load CCprofiler and setup the work environment
library(devtools)
options(warn=-1)

install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
library('CCprofiler')
library('data.table')
library('ggplot2')

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Mouse/")
knitr::opts_knit$set(
  root.dir = '/Users/isabell/Desktop/projects/ProteoformProject/Mouse/')

options(knitr.table.format = function() {
  if (knitr::is_latex_output()) 
    "latex" else "pandoc"
})

#' ## Read input data 
input_data <- fread("input_data.txt") 
knitr::kable(head(input_data, n=3))

input_data_annotation <- fread("input_data_annotation.txt") 
knitr::kable(head(input_data_annotation, n=3))

#' ## Craete traces object:
traces <- importPCPdata(input_data = input_data, 
                        fraction_annotation = input_data_annotation)

#' ## Annotate traces with uniprot information
trace_annotation <- fread(
  "../InputNonFrac/uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_+AND+revie--.tab")

traces <- annotateTraces(traces,
                         trace_annotation,
                         traces_id_column = "protein_id",
                         trace_annotation_id_column = "Entry")

#' ## Inspect traces list:
summary(traces)

#' ## Annotate peptide sequence positions
traces_iso_pos <- annotatePeptideSequences(
  traces, 
  fasta_file="../InputNonFrac/uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_+AND+revie--.fasta")

#' ## Select representative peptides
traces_start_end <- summarizeAlternativePeptideSequences(
  traces_iso_pos, topN=1,start="PeptidePositionStart",end="PeptidePositionEnd", verbose=F)

## Remove traces with 0 standard deviation (can't be clustered)
zerovar <- apply(getIntensityMatrix(traces_start_end), 1, var) > 0
traces_zerovar <- subset(traces_start_end,
                        trace_subset_ids = names(zerovar[zerovar]))

#' ## Remove single peptide genes
traces_multiPep <- filterSinglePeptideHits(traces_zerovar)

#' ## Remove outlier protein 'A2ASS6' with too many peptides
valid_proteins <- unique(traces_multiPep$trace_annotation$protein_id)
valid_proteins <- valid_proteins[valid_proteins != 'A2ASS6']
traces_filtered <- subset(traces_multiPep, 
                          trace_subset_ids = valid_proteins, 
                          trace_subset_type = "protein_id")

#' ## Estimate proteoform scores by correlation clustering
traces_corr <- calculateGeneCorrMatrices(traces_filtered)

traces_clustered <- clusterPeptides(
  traces_corr,
  method = "average", plot = T, PDF=T,
  name="ProteoformClusters_mouse")

traces_clusteredInN <- cutClustersInNreal(traces_clustered,
                                           clusterN = 2,
                                           min_peptides_per_cluster = 2)

traces_scored <- calculateProteoformScore(traces_clusteredInN)

traces_scored$trace_annotation[, proteoform_score := ifelse(is.na(proteoform_score), 0, proteoform_score)]
traces_scored$trace_annotation[, proteoform_score_pval_adj := ifelse(is.na(proteoform_score_pval_adj), 1, proteoform_score_pval_adj)]

plotProteoformPvalHist(traces_scored, name="traces_scored_pval_hist", PDF=T)

plotProteoformVolcano(traces_scored, name="traces_scored_volcano", PDF=T)

source("../CCprofilerAnalysis/thesis/PaperAnalysis/CellCycleHela/plotProteoformVolcanoLines.R")
plotProteoformVolcanoLines(traces_scored, name="traces_scored_volcano", score_cutoff = 0.1, adj_pval_cutoff = 0.1, PDF=T)

#' ## Proteoform annotation
#traces_proteoforms <- annotateTracesWithProteoforms(
#  traces_scored, score_cutoff = 0)
traces_proteoforms <- annotateTracesWithProteoforms(
  traces_scored, score_cutoff = 0, adj_pval_cutoff =  1)

#' ## Evaluate sequence location of the determined proteoforms
traces_location <- evaluateProteoformLocation(
  traces_proteoforms, name=paste0("NormalizedSD_",2))

#' ## Save traces 
saveRDS(traces_location,"traces_location.rds")


#' ## Plot proteoform profile for an example protein
plotPeptideCluster(
  traces_location, "Q9JKS4", closeGaps=T, PDF=F)

plotSub <- subset(
  traces_location, trace_subset_ids = "Q9JKS4", 
  trace_subset_type = "protein_id")

plotSub <- plotSub[c("traces","trace_type","trace_annotation","fraction_annotation")]
class(plotSub) <- 'traces'

plot(plotSub, legend = T, colour_by = "proteoform_id")

#' ## Plot all proteoform profiles
sigProteins <- subset(traces_location$trace_annotation, (proteoform_score>0.1) & (proteoform_score_pval_adj<=0.1))
sigProteins <- unique(sigProteins[order(proteoform_score, decreasing = T)]$protein_id)

pdf("allProteoformClusters_score01_pval01.pdf",height=4, width=6)
for (p in sigProteins){
  plotPeptideCluster(
    traces_location,p, closeGaps=T, PDF=F)
  plotSub <- subset(
    traces_location, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plotSub <- plotSub[c("traces","trace_type","trace_annotation","fraction_annotation")]
  class(plotSub) <- 'traces'
  plot(plotSub, legend = T,
       name=paste0(p, "\n score=", round(plotSub$trace_annotation$proteoform_score[1], digits = 3),
                   "\n adj.pval=", round(plotSub$trace_annotation$proteoform_score_pval_adj[1], digits = 3)),
       colour_by = "proteoform_id")
}
dev.off()

## Plot all proteoform profiles with high sequence proximity
intCases <- subset(traces_location$trace_annotation, (proteoform_score>0.1) & (proteoform_score_pval_adj<=0.1) & genomLocation_pval_lim_min<=0.05)
intCases <- unique(intCases[order(proteoform_score, decreasing = T)]$protein_id)

pdf("allProteoformClusters_score01_pval01_pval_lim_min_005.pdf",height=4, width=6)
for (p in intCases){
  plotPeptideCluster(
    traces_location, p, closeGaps=T, PDF=F)
  plotSub <- subset(
    traces_location, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plotSub <- plotSub[c("traces","trace_type","trace_annotation","fraction_annotation")]
  class(plotSub) <- 'traces'
  plot(plotSub, legend = T, 
       name=paste0(p, "\n score=", round(plotSub$trace_annotation$proteoform_score[1], digits = 3),
                   "\n adj.pval=", round(plotSub$trace_annotation$proteoform_score_pval_adj[1], digits = 3)),
       colour_by = "proteoform_id")
}
dev.off()


# @ToDo write tables with proteins for enrichment analysis
getProteoformStats <- function(traces, adj_pval_cutoff, localization_threshold=0.05){
  tab <- traces$trace_annotation
  
  proteins <- unique(tab$protein_id)
  n_proteins <- length(proteins)
  write.table(proteins,"ProteinTables/proteins_total.txt", sep="\t", quote = F, col.names = F, row.names = F)
  
  tab_sub <- subset(tab, (proteoform_score_pval_adj<=adj_pval_cutoff) & (proteoform_score >= 0.1))
  
  proteins_proteoforms <- unique(tab_sub$protein_id)
  n_proteins_proteoforms <- length(proteins_proteoforms)
  write.table(proteins_proteoforms,paste0("ProteinTables/proteoforms_adj_pval_cutoff_",adj_pval_cutoff,".txt"), sep="\t", quote = F, col.names = F, row.names = F)
  
  proteins_notExtreme <- unique(tab_sub[(genomLocation_pval_lim_min > localization_threshold) & (genomLocation_pval_min > localization_threshold)]$protein_id)
  proteins_extreme <- unique(tab_sub[(genomLocation_pval_lim_min <= localization_threshold) | (genomLocation_pval_min <= localization_threshold)]$protein_id)
  proteins_moreExtreme <- unique(tab_sub[(genomLocation_pval_lim_min <= localization_threshold) & (genomLocation_pval_min > localization_threshold)]$protein_id)
  proteins_equallyExtreme <- unique(tab_sub[(genomLocation_pval_min <= localization_threshold)]$protein_id)
  
  n_proteins_notExtreme <- length(proteins_notExtreme)
  n_proteins_extreme <- length(proteins_extreme)
  n_proteins_moreExtreme <- length(proteins_moreExtreme)
  n_proteins_equallyExtreme <- length(proteins_equallyExtreme)
  
  dt <- data.table(adj_pval_cutoff = adj_pval_cutoff,
                   localization_threshold = localization_threshold,
                   n_proteins=n_proteins,
                   n_proteins_proteoforms=n_proteins_proteoforms,
                   p_proteins_proteoforms=n_proteins_proteoforms/n_proteins,
                   n_proteins_notExtreme=n_proteins_notExtreme,
                   p_proteins_notExtreme=n_proteins_notExtreme/n_proteins_proteoforms,
                   n_proteins_extreme=n_proteins_extreme,
                   p_proteins_extreme=n_proteins_extreme/n_proteins_proteoforms,
                   n_proteins_moreExtreme=n_proteins_moreExtreme,
                   p_proteins_moreExtreme=n_proteins_moreExtreme/n_proteins_proteoforms,
                   n_proteins_equallyExtreme=n_proteins_equallyExtreme,
                   p_proteins_equallyExtreme=n_proteins_equallyExtreme/n_proteins_proteoforms)
  return(dt)
}


res <- data.table()
for (i in seq(0,1,0.05)) {
  for (j in seq(0, 0.5, 0.05)){
    res_i <- getProteoformStats(traces_location, adj_pval_cutoff=i, localization_threshold=j)
    res <- rbind(res, res_i)
  }
}

protein_stats = subset(res, select=c("adj_pval_cutoff", "localization_threshold", "n_proteins","n_proteins_proteoforms","n_proteins_extreme","n_proteins_equallyExtreme","n_proteins_moreExtreme"))

knitr::kable(protein_stats)
write.table(protein_stats, "protein_proteoform_stats.txt", 
            quote=F, sep="\t", row.names = F, col.names = T)


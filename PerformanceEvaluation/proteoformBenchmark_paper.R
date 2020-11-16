#' ---
#' title: "Proteoform annotation benchmark"
#' author: "Isabell Bludau"
#' output:
#'   html_document:
#'     toc: true
#'   html_notebook:
#'     toc: true
#'   pdf_document:
#'     toc: true
#' ---

#' # Overview
#' The goal of this workflow is to benchmark proteoform detection in CCprofiler. 
#' For this purpose, we perform an in silico generation of mixture proteins with 
#' multiple proteoforms by mixing peptides of different proteins. These mixture 
#' proteins are then used to evaluate the ability of our algorithm to 
#' (a) correctly determine the mixture proteins and 
#' (b) correctly group the peptides of mixture proteins by their parental proteins.
# rmarkdown::render("proteoformBenchmark.R", "pdf_document")

#' # Load CCprofiler package, set working directory
#source(
#  "/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")
#library(CCprofiler)
library(data.table)
library(ggplot2)
#devtools::load_all("~/Desktop/projects/ProteoformProject/CCprofiler/")
install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
library(CCprofiler)

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Benchmark")
knitr::opts_knit$set(
  root.dir = '/Users/isabell/Desktop/projects/ProteoformProject/Benchmark')
options(knitr.table.format = function() {
  if (knitr::is_latex_output()) 
    "latex" else "pandoc"
})

#' Source functions required for benchmarking
source("../CCprofilerAnalysis/thesis/generateMixtureTraces_final.R")

#' Load data
peptide_traces <- readRDS("HEK_raw_pepTraces.rds")

#' Basic preprocessing
# Remove decoys from the peptdide traces
targets <- peptide_traces$trace_annotation$protein_id[grep("DECOY",
           peptide_traces$trace_annotation$protein_id, invert = T)]
peptide_traces <- subset(peptide_traces, trace_subset_ids = targets,
                         trace_subset_type = "protein_id",)
# Remove single peptide hits
peptide_traces <- filterSinglePeptideHits(peptide_traces)

#' To benchmark the sensitivity of detecting genes with multiple proteoforms 
#' we generate 5000 mixture traces (peptide_traces_mixed) and append them to the normal 
#' peptide traces (negative_traces). The resulting combined traces set (benchmark_traces) 
#' is subsequently used for proteoform detection.
#' The recovery of true mixture proteins is evaluated across different score cutoffs.
#' 
#' We mix peptides from 2 parental proteins. For mixing, only 50% 
#' (or minimally 2) peptides are selected per parental protein.   

i = 2

peptide_traces_mixed <- generateMixtureTraces(peptide_traces, 
                                              trueInteractions = NULL,
                                              n_proteins = i, 
                                              n_mixtureTraces = 3000, 
                                              peptide_ratio = 0.5, 
                                              min_peptide_count = 2, 
                                              seed = 123)

negative_traces <- copy(peptide_traces)
negative_traces$trace_annotation[,n_mixed_proteins := 1]
negative_traces$trace_annotation[,is_mixed:=FALSE]

benchmark_traces <- copy(negative_traces)
benchmark_traces$traces <- rbind(
  negative_traces$traces, peptide_traces_mixed$traces)
benchmark_traces$trace_annotation <- rbind(
  negative_traces$trace_annotation, peptide_traces_mixed$trace_annotation)
 
mixed_prot = unique(benchmark_traces$trace_annotation[is_mixed==TRUE]$protein_id)[44]
plot_sub <- subset(benchmark_traces, trace_subset_ids = mixed_prot, trace_subset_type = "protein_id")
pdf("benchmark_trace_mixed.pdf", height = 3.5, width=4)
p <- plot(plot_sub, colour_by = 'Entry_name', monomer_MW = FALSE)
print(p)
dev.off()

single_prot = unique(benchmark_traces$trace_annotation[is_mixed==FALSE]$protein_id)[444]
plot_sub <- subset(benchmark_traces, trace_subset_ids = single_prot, trace_subset_type = "protein_id")
pdf("benchmark_trace_single.pdf", height = 3, width=4)
p <- plot(plot_sub, colour_by = 'Entry_name', monomer_MW = FALSE)
print(p)
dev.off()


benchmark_traces_corr <- calculateGeneCorrMatrices(benchmark_traces)

benchmark_traces_clustered <- clusterPeptides(
  benchmark_traces_corr,
  method = "average", plot = T, PDF=T,
  name="ProteoformClusters_benchmark")

benchmark_traces_cut <- cutClustersInNreal(benchmark_traces_clustered,
                                           clusterN = i,
                                           min_peptides_per_cluster = 2)

benchmark_traces_scored <- calculateProteoformScore(benchmark_traces_cut)

benchmark_traces_score_dist = subset(benchmark_traces_scored$trace_annotation, select=c(proteoform_score,is_mixed,protein_id))
benchmark_traces_score_dist <- unique(benchmark_traces_score_dist)

pdf('proteoform_score_dist.pdf', width=5, height=4.5)
p <- ggplot(benchmark_traces_score_dist, aes(x=proteoform_score, fill=is_mixed, color=is_mixed)) + 
  geom_histogram(stat = "bin", binwidth=0.03, position='identity', alpha=0.3) +
  #scale_x_continuous(limits = c(0, 1.25), expand = c(0.02, 0.02)) +
  #scale_y_continuous(limits = c(0, 15000), expand = c(0.04, 0.04)) +
  xlab('score') +
  scale_color_manual(values=c("#999999", "green4")) +
  scale_fill_manual(values=c("#999999", "green4")) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.8))
print(p)
dev.off()

max_score = round(max(benchmark_traces_scored$trace_annotation$proteoform_score, na.rm = TRUE), digits=2) + 0.05
summaryStats_mixedProteinClustering <- data.table()

for (score_cutoff_i in seq(0, max_score, 0.05)){
  cluster_stats <- evaluateProteoformClusteringOfMixtureTraces(benchmark_traces_scored, score_cutoff = score_cutoff_i)
  
  # remove proteins with only outlier peptides
  cluster_stats_noOutliers <- subset(cluster_stats, n_proteins==n_proteins_noOutliers)
  
  # count single proteins
  n_single_proteins <- length(unique(subset(cluster_stats_noOutliers, n_proteins==1)$protein))
  # count correct and incorrect single proteins
  n_single_proteins_correct <- length(unique(subset(cluster_stats_noOutliers, n_proteins==1 & n_proteoforms==1)$protein))
  n_single_proteins_incorrect <- length(unique(subset(cluster_stats_noOutliers, n_proteins==1 & n_proteoforms==2)$protein))
  
  # count mixed proteins
  n_mixed_proteins <- length(unique(subset(cluster_stats_noOutliers, n_proteins==2)$protein))
  # count correct and incorrect single proteins
  n_mixed_proteins_correct <- length(unique(subset(cluster_stats_noOutliers, n_proteins==2 & n_proteoforms==2)$protein))
  n_mixed_proteins_incorrect <- length(unique(subset(cluster_stats_noOutliers, n_proteins==2 & n_proteoforms==1)$protein))
  
  # Calculate stats
  TPR = n_mixed_proteins_correct/n_mixed_proteins
  FPR = n_single_proteins_incorrect/n_single_proteins
  
  # Count mistakes in correct proteoform assignment
  n_mixed_proteins_correct_noMistakes <- length(unique(subset(cluster_stats_noOutliers, n_proteins==2 & n_proteoforms==2 & n_mistakes==0)$protein))
  p_mixed_proteins_correct_noMistakes <- n_mixed_proteins_correct_noMistakes/n_mixed_proteins_correct

  # Report results
  summaryStats_mixedProteinClustering <- rbind(
    summaryStats_mixedProteinClustering, 
    data.table(
      score=score_cutoff_i,
      n_single_proteins=n_single_proteins,
      n_single_proteins_correct=n_single_proteins_correct,
      n_single_proteins_incorrect=n_single_proteins_incorrect,
      n_mixed_proteins=n_mixed_proteins,
      n_mixed_proteins_correct=n_mixed_proteins_correct,
      n_mixed_proteins_incorrect=n_mixed_proteins_incorrect,
      TPR=TPR,
      FPR=FPR,
      n_mixed_proteins_correct_noMistakes=n_mixed_proteins_correct_noMistakes,
      p_mixed_proteins_correct_noMistakes=p_mixed_proteins_correct_noMistakes
      ))
}

setnames(summaryStats_mixedProteinClustering, 'p_mixed_proteins_correct_noMistakes', 'p_correct' )

pdf('ROC_curve.pdf', width=5, height=4.5)
p <- ggplot(summaryStats_mixedProteinClustering, aes(x=FPR, y=TPR)) +
  geom_abline(intercept = 0, slope = 1, color='grey') +
  geom_line(color='grey35') +
  #scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0.02)) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0.02)) +
  geom_point(aes(color=score, size=p_correct)) +
  scale_color_gradient(low="palegreen1", high="darkgreen", limits=c(0,max_score)) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.42))
print(p)
dev.off()

write.table(summaryStats_mixedProteinClustering,'summaryStats_mixedProteinClustering.tsv',sep='\t',quote=F,row.names = F, col.names = T)

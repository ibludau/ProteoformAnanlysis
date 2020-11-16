#' ## Load CCprofiler and setup the work environment
library(devtools)
options(warn=-1)

#install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
devtools::load_all("~/Desktop/projects/ProteoformProject/CCprofiler/")
#library('CCprofiler')
library('data.table')
library('ggplot2')

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Simulation/")

simulate_peptide_profiles <- function(n_pep=10,
                                      n_pep_mixed=2,
                                      ratio=2,
                                      prot_mean=10,
                                      prot_sd=1,
                                      pep_sd=0.5,
                                      n_samples=100,
                                      f_samples_mixed=0.5,
                                      prot_i=1){
  
  prot_ints = rnorm(n_samples,prot_mean,prot_sd)
  if (n_pep_mixed>=n_pep){
    stop("n_pep must be > n_pep_mixed")
  }
  n_pep_single = n_pep-n_pep_mixed
  n_samples_single = round(n_samples*(1-f_samples_mixed))
  n_samples_mixed = round(n_samples*(f_samples_mixed))
  
  pep_ints = list()
  for(i in seq(1,n_pep,1)){pep_ints[[i]] = unlist(lapply(prot_ints, function(x){rnorm(1,x,pep_sd)}))}
  if (n_pep_mixed > 0) {
    for(i in seq(1,n_pep_mixed,1)){
      pep_ints[[i]][1:n_samples_mixed] = pep_ints[[i]][1:n_samples_mixed]*ratio
    }
  }
  
  pep_ints_dt <- data.table(t(as.data.table(pep_ints)))
  colnames(pep_ints_dt) = as.character(seq(1,n_samples,1))
  pep_ints_dt[,id:=paste0('prot_',prot_i,'_pep_',seq(1,n_pep,1))]
  
  return(pep_ints_dt)
}

simulate_traces <- function(n_proteins_single=100,
                            n_proteins_mixed=100,
                            n_pep=10,
                            n_pep_mixed=0.5,
                            ratio=2,
                            prot_mean=10,
                            prot_sd=1,
                            pep_sd=0.5,
                            n_samples=100,
                            f_samples_mixed=0.5){
  
  traces_mixed <- data.table()
  traces_ann_mixed <- data.table()
  for (i in seq(1,n_proteins_mixed,1)){
    traces_i = simulate_peptide_profiles(n_pep=n_pep,
                                         n_pep_mixed=n_pep_mixed,
                                         ratio=ratio,
                                         prot_mean=prot_mean,
                                         prot_sd=prot_sd,
                                         pep_sd=pep_sd,
                                         n_samples=n_samples,
                                         f_samples_mixed=f_samples_mixed,
                                         prot_i = paste0("prot_mixed_",i))
    traces_ann_i = data.table(id=traces_i$id,
                              protein_id=paste0("prot_mixed_",i), 
                              is_mixed=TRUE,
                              n_mixed_proteins=2,
                              Entry_name=c(rep(paste0("prot_mixed_",i,"_proteoform_1"),n_pep_mixed),
                                           rep(paste0("prot_mixed_",i,"_proteoform_2"),n_pep-n_pep_mixed)))
    traces_mixed <- rbind(traces_mixed, traces_i)
    traces_ann_mixed <- rbind(traces_ann_mixed, traces_ann_i)
  }
  
  traces_single <- data.table()
  traces_ann_single <- data.table()
  for (i in seq(1,n_proteins_single,1)){
    traces_i = simulate_peptide_profiles(n_pep=n_pep,
                                         n_pep_mixed=0,
                                         ratio=ratio,
                                         prot_mean=prot_mean,
                                         prot_sd=prot_sd,
                                         pep_sd=pep_sd,
                                         n_samples=n_samples,
                                         f_samples_mixed=f_samples_mixed,
                                         prot_i = paste0("prot_single_",i))
    traces_ann_i = data.table(id=traces_i$id,
                              protein_id=paste0("prot_single_",i), 
                              is_mixed=FALSE,
                              n_mixed_proteins=1,
                              Entry_name=paste0("prot_single_",i))
    traces_single <- rbind(traces_single, traces_i)
    traces_ann_single <- rbind(traces_ann_single, traces_ann_i)
  }
  
  traces <- rbind(traces_single, traces_mixed)
  traces_ann <- rbind(traces_ann_single, traces_ann_mixed)
  
  fraction_ann = data.table(id=seq(1,n_samples,1))
  
  traces_obj <- list(traces=traces,
                     trace_type="peptide",
                     trace_annotation=traces_ann,
                     fraction_annotation=fraction_ann)
  class(traces_obj) <- "traces"
  
  return(traces_obj)
}


source("../CCprofilerAnalysis/thesis/generateMixtureTraces_final.R")

evaluateProteoformClusteringOfMixtureTraces <- function(traces, score_cutoff){
  all_proteins <- unique(traces$trace_annotation$protein_id)
  mistakes <- lapply(all_proteins, function(x){
    sub <- subset(traces$trace_annotation, protein_id==x)
    
    sub_noOutliers <- subset(sub, cluster != 100)
    
    n_outlierPeptides <- nrow(subset(sub, cluster == 100))
    
    unique_proteins <- unique(sub$Entry_name)
    n_proteins <- length(unique_proteins)
    
    unique_proteins_noOutliers <- unique(sub$Entry_name)
    n_proteins_noOutliers <- length(unique_proteins_noOutliers)
    
    if (nrow(sub_noOutliers) ==0) {
      
      n_proteoforms <- 1
      
      if (n_proteins == 1) {
        n_mistakes = 0
      } else {
        n_mistakes = nrow(sub)
      }
      
    } else {
      
      sub_noOutliers[, proteoform_id := ifelse(proteoform_score >= score_cutoff, paste(protein_id, cluster, sep ='_'), protein_id)]
      
      unique_proteoforms <- unique(sub_noOutliers$proteoform_id)
      n_proteoforms <- length(unique_proteoforms)
      
      n_mistakes <- sum(unlist(lapply(unique_proteoforms,function(p){
        subSub <- subset(sub_noOutliers,proteoform_id==p)
        nrow(subSub)-max(table(subSub$Entry_name))
      })))
      
    }
    
    
    return(data.table(protein=x, 
                      n_outlierPeptides=n_outlierPeptides, 
                      n_proteins=n_proteins, 
                      n_proteins_noOutliers=n_proteins_noOutliers, 
                      n_proteoforms=n_proteoforms, 
                      n_mistakes=n_mistakes))
  })
  mistakes <- do.call(rbind,mistakes)
  return(mistakes)
}

n_pep_mixed_c <- c(2,5)
pep_sd_c <- c(1.5,2,2.5,3,3.5,4) 

parameter_grid <- as.data.table(expand.grid(n_pep_mixed_c,pep_sd_c))
colnames(parameter_grid) <- c("n_pep_mixed_i","pep_sd_i")


summaryStats_mixedProteinClustering_all <- data.table()

for (i in seq(1,nrow(parameter_grid),1)) {
  traces <- simulate_traces(n_proteins_single=1000,
                            n_proteins_mixed=1000,
                            n_pep=10,
                            n_pep_mixed=parameter_grid$n_pep_mixed_i[i],
                            ratio=2,
                            prot_mean=10,
                            prot_sd=1,
                            pep_sd=parameter_grid$pep_sd_i[i],
                            n_samples=100,
                            f_samples_mixed=0.5)
  
  # test plots
  prot_mixed_1 <- subset(traces, trace_subset_ids = "prot_mixed_1", trace_subset_type = "protein_id")
  pdf(paste0("prot_mixed_",parameter_grid$n_pep_mixed_i[i],"_",parameter_grid$pep_sd_i[i],".pdf"), height = 3, width=4)
    print(plot(prot_mixed_1, colour_by = 'Entry_name'))
  dev.off()
  
  prot_single_1 <- subset(traces, trace_subset_ids = "prot_single_1", trace_subset_type = "protein_id")
  pdf(paste0("prot_single_",parameter_grid$n_pep_mixed_i[i],"_",parameter_grid$pep_sd_i[i],".pdf"), height = 3, width=4)
    print(plot(prot_single_1, colour_by = 'Entry_name'))
  dev.off()
  
  # Perform benchmark
  benchmark_traces <- copy(traces)
  
  benchmark_traces_corr <- calculateGeneCorrMatrices(benchmark_traces)
  
  benchmark_traces_clustered <- clusterPeptides(
    benchmark_traces_corr,
    method = "average", plot = F, PDF=F,
    name="ProteoformClusters_benchmark")
  
  benchmark_traces_cut <- cutClustersInNreal(benchmark_traces_clustered,
                                             clusterN = 2,
                                             min_peptides_per_cluster = 2)
  
  benchmark_traces_scored <- calculateProteoformScore(benchmark_traces_cut)
  
  benchmark_traces_scored$trace_annotation[which(is.na(benchmark_traces_scored$trace_annotation[]$proteoform_score))]$proteoform_score=0
  
  benchmark_traces_score_dist = subset(benchmark_traces_scored$trace_annotation, select=c(proteoform_score,is_mixed,protein_id))
  benchmark_traces_score_dist <- unique(benchmark_traces_score_dist)
  
  pdf(paste0("proteoform_score_dist_",parameter_grid$n_pep_mixed_i[i],"_",parameter_grid$pep_sd_i[i],".pdf"), width=2.5, height=2.25)
  p <- ggplot(benchmark_traces_score_dist, aes(x=proteoform_score, fill=is_mixed, color=is_mixed)) + 
    geom_histogram(stat = "bin", binwidth=0.03, position='identity', alpha=0.3) +
    xlab('score') +
    xlim(-0.1,1) +
    scale_color_manual(values=c("#999999", "green4")) +
    scale_fill_manual(values=c("#999999", "green4")) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  print(p)
  dev.off()
  
  max_score = round(max(benchmark_traces_scored$trace_annotation$proteoform_score, na.rm = TRUE), digits=2) + 0.05
  summaryStats_mixedProteinClustering <- data.table()
  
  for (score_cutoff_i in seq(0, max_score, 0.01)){
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
  
  summaryStats_mixedProteinClustering[,n_pep_mixed := parameter_grid$n_pep_mixed_i[i]]
  summaryStats_mixedProteinClustering[,pep_sd := parameter_grid$pep_sd_i[i]]
  
  #write.table(summaryStats_mixedProteinClustering,paste0('summaryStats_mixedProteinClustering_',parameter_grid$n_pep_mixed_i[i],"_",parameter_grid$pep_sd_i[i],'.tsv'),sep='\t',quote=F,row.names = F, col.names = T)
  
  summaryStats_mixedProteinClustering_all <- rbind(summaryStats_mixedProteinClustering_all,summaryStats_mixedProteinClustering)
}

summaryStats_mixedProteinClustering_all
write.table(summaryStats_mixedProteinClustering_all,'summaryStats_mixedProteinClustering_all.tsv',sep='\t',quote=F,row.names = F, col.names = T)

summaryStats_mixedProteinClustering_all[,group:=paste0(n_pep_mixed,"_",pep_sd)]

summaryStats_mixedProteinClustering_all$n_pep_mixed = as.character(summaryStats_mixedProteinClustering_all$n_pep_mixed)
summaryStats_mixedProteinClustering_all$pep_sd = as.character(summaryStats_mixedProteinClustering_all$pep_sd)


summaryStats_sub <- summaryStats_mixedProteinClustering_all[n_pep_mixed %in% c(2,5)]

pdf("simulation_ROC_npep_sd.pdf", width=5, height=4.5)
ggplot(summaryStats_sub, aes(x=FPR, y=TPR, colour=n_pep_mixed, shape=pep_sd, group=group)) + 
  geom_abline(intercept = 0, slope = 1, color='grey') +
  geom_line(color='grey35') +
  geom_point() +
  geom_line() +
  theme_classic() +
  xlim(0,1) +
  ylim(0,1) +
  scale_color_manual(values=c("#008080","#32cd32")) +
  theme(legend.position = c(0.9, 0.42))
dev.off()


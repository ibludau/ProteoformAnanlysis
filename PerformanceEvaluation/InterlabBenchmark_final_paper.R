#' ## Load CCprofiler and setup the work environment
library(devtools)
options(warn=-1)

install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
#devtools::load_all("~/Desktop/projects/ProteoformProject/CCprofiler/")
library('CCprofiler')
library('data.table')
library('ggplot2')

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Revisions/")

#####################################
### Create dataset ##################
#####################################

dt <- fread("site02_global_q_0.01_applied_to_local_global.txt")

# rename protein
dt[,ProteinName:=gsub("1/","",ProteinName)]
dt <- subset(dt, ProteinName != "iRT_protein") # remove iRT peptides
dt <- subset(dt, ProteinName != "sp|AQUA30|AQUA30") # remove AQUA peptides

# rename runs
dt[,run:=gsub("Site2_AQUA_HEK_","",run_id)]
dt[,run:=gsub("_180714.mzXML.gz","",run)]

# determine day
dt[,day:=gsub("S.*_SW_","",run)]
dt[,day:=gsub("_.*","",day)]

# subset to important columns
dt_sub <- subset(dt, select=c("run","day","ProteinName","FullPeptideName","peptide_group_label", "Intensity"))

# aggregate charge states
dt_sub[,pep_int:=sum(Intensity), by=c("ProteinName","FullPeptideName","run")]
dt_sub <- unique(subset(dt_sub, select=c("run","day","ProteinName","FullPeptideName","pep_int")))

# subset to sufficient data per condition
dt_sub[,n_day1 := sum(day=="day1"), by=c("FullPeptideName")]
dt_sub[,n_day3 := sum(day=="day3"), by=c("FullPeptideName")]
dt_sub[,n_day5 := sum(day=="day5"), by=c("FullPeptideName")]
dt_sub[,min_n_per_day := min(n_day1,n_day3,n_day5),by=c("FullPeptideName")]
dt_sub <- subset(dt_sub, min_n_per_day==7)


# subset to proteins with > 4 peptides:
dt_sub[,n_pep:=length(unique(FullPeptideName)),by=c("ProteinName")]
dt_sub <- dt_sub[n_pep>=4]

# median normalization
dt_sub[,log2_int:=log2(pep_int)]
dt_sub[,median_perRun:=median(log2_int), by="run"]
dt_sub[,median_median:=mean(median_perRun)]
dt_sub[,diff_median:=median_median-median_perRun]
dt_sub[,norm_log2_int := log2_int+diff_median]
dt_sub[,norm_int := 2^norm_log2_int]

# introduce variation
set.seed(1)
dt_sub[,diff_fac_3 := runif(1, min = 1, max = 6), by="ProteinName"]
set.seed(2)
dt_sub[,diff_fac_5 := runif(1, min = 1, max = 6), by="ProteinName"]
#dt_sub[,diff_norm_int:=ifelse(day %in% c("day3","day5"), diff_fac*norm_int, norm_int)]
dt_sub[,diff_norm_int:=ifelse(day == "day3", diff_fac_3*norm_int, norm_int)]
dt_sub[,diff_norm_int:=ifelse(day == "day5", diff_fac_5*norm_int, diff_norm_int)]

# randomly select 1000 proteins to perturb
set.seed(22)
proteins_to_perturb = sample(unique(dt_sub$ProteinName),1000)
dt_sub[,perturbed_protein := ifelse(ProteinName %in% proteins_to_perturb, TRUE, FALSE), by="ProteinName"]

# determine reduction factor for each protein
# sample from uniform distribution
set.seed(44)
dt_sub[,red_fac:=ifelse(perturbed_protein, runif(1, min = 0.01, max = 0.90), 1), by="ProteinName"]

interlab_benchmark_data <- copy(dt_sub)

saveRDS(interlab_benchmark_data, "interlab_benchmark_data.rds")


##################################
### Introduce pep. perturbation ##
##################################

rm(list=ls())
gc()
interlab_benchmark_data <- readRDS("interlab_benchmark_data.rds")

generatePerturbedProfiles <- function(input_data, nf_peptides_to_perturb="random"){
  dt_input <- copy(input_data)
  
  set.seed(66)
  
  if (nf_peptides_to_perturb == "random") {
    dt_input[,max_perturbed_peptides := ceiling(n_pep*0.5), by="ProteinName"]
    dt_input[,n_perturbed_peptides := max(sample(seq(2,max(2,max_perturbed_peptides),1), 1),2), by="ProteinName"]
  } else if (nf_peptides_to_perturb >= 1) {
    dt_input[,n_perturbed_peptides := nf_peptides_to_perturb, by="ProteinName"]
  } else if (nf_peptides_to_perturb < 1) {
    dt_input[,n_perturbed_peptides := ceiling(n_pep*nf_peptides_to_perturb), by="ProteinName"]
    dt_input[,n_perturbed_peptides := ifelse(n_perturbed_peptides<2,2,n_perturbed_peptides), by="ProteinName"]
  }
  
  dt_input[,perturbed_peptides := paste(unique(FullPeptideName)[c(1:n_perturbed_peptides)],collapse=";"), by = c("ProteinName")]
  
  # reduce day 5 by reduction factor for peptides in perturbed peptides
  dt_input[,perturbed_peptide:=((FullPeptideName %in% unlist(strsplit(perturbed_peptides,";"))) & (perturbed_protein==TRUE)), by=c("FullPeptideName","ProteinName")]
  
  dt_input[,mod_pep_int:=ifelse(((day=="day5") & (perturbed_peptide) & (perturbed_protein)), diff_norm_int*red_fac, diff_norm_int)]
  
  setnames(dt_input, c("FullPeptideName","ProteinName","run","mod_pep_int"), c("peptide_id","protein_id","filename","intensity"))
  
  # format for CCprofiler
  input_data <- subset(dt_input, select=c("peptide_id","protein_id","filename","intensity"))
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  # Plot
  pdf(paste0("plot_example_perturbed_proteins_",nf_peptides_to_perturb,".pdf"), width = 5, height = 3)
  col <- c("TRUE"="magenta4", "FALSE"="#6699cc")
  perturbed_proteins <- unique(dt_input[perturbed_protein==TRUE]$protein_id)
  for (i in c(1:100)){
    x <- dt_input[protein_id==perturbed_proteins[i]]
    
    x$filename <- factor(x$filename, levels=input_data_annotation$filename)
    
    filename_label <- c(rep("day1", 7),rep("day3", 7),rep("day5", 7))
    rep_label <- rep(seq(1,7,1),3)
    filename_label <- paste(filename_label, rep_label, sep="_")
    
    p <- ggplot(x, aes(x=filename, y=intensity, color=perturbed_peptide, group=peptide_id)) + 
      geom_point(alpha=0.9) +
      geom_line(alpha=0.7) +
      theme_classic() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = col, name="perturbed peptide") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_x_discrete(labels=filename_label) +
      xlab("sample") +
      ggtitle(paste0(perturbed_proteins[i]," \nperturbation factor = ",round(unique(x$red_fac),digits = 3)))
    print(p)
  }
  dev.off()
  
  pdf(paste0("plot_example_non_perturbed_proteins_",nf_peptides_to_perturb,".pdf"), width = 5, height = 3)
  proteins_not_perturbed = unique(dt_input[perturbed_protein == FALSE]$protein_id)
  col <- c("TRUE"="magenta4", "FALSE"="#6699cc")
  for (i in c(1:50)){
    x <- dt_input[protein_id==proteins_not_perturbed[i]]
    
    x$filename <- factor(x$filename, levels=input_data_annotation$filename)
    
    filename_label <- c(rep("day1", 7),rep("day3", 7),rep("day5", 7))
    rep_label <- rep(seq(1,7,1),3)
    filename_label <- paste(filename_label, rep_label, sep="_")
    
    p <- ggplot(x, aes(x=filename, y=intensity, color=perturbed_peptide, group=peptide_id)) + 
      geom_point(alpha=0.9) +
      geom_line(alpha=0.7) +
      theme_classic() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = col, name="perturbed peptide") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_x_discrete(labels=filename_label) +
      xlab("sample") +
      ggtitle(paste0(perturbed_proteins[i]," \nperturbation factor = ",round(unique(x$red_fac),digits = 3)))
    print(p)
  }
  dev.off()
  
  
  return(dt_input)
}

random_perm_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb="random")
perm_1pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=1)
perm_2pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=2)
perm_025pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=0.25)
perm_050pep_data <- generatePerturbedProfiles(input_data=interlab_benchmark_data, nf_peptides_to_perturb=0.5)

#####################################
### CCprofiler analysis #############
#####################################

library(ggpubr)
library(gridExtra)
library(grid)

performCCprofilerAnalysis <- function(data, name, score_cutoff = 0.1, adj_pval_cutoff =  0.1){
  
  dt_input <- copy(data)
  
  input_data_annotation <- unique(subset(dt_input,select=c("filename","day")))
  setorderv(input_data_annotation, c("day","filename"))
  input_data_annotation[,fraction_number:=.I]
  input_data_annotation <- subset(input_data_annotation,select=c("filename","fraction_number"))
  
  traces <- importPCPdata(input_data = dt_input,
                          fraction_annotation = input_data_annotation)
  
  trace_annotation <- unique(subset(dt_input, select=c("peptide_id",
                                                       "n_pep","n_perturbed_peptides",
                                                       "perturbed_protein","perturbed_peptide","red_fac")))
  
  traces <- annotateTraces(traces,
                           trace_annotation,
                           traces_id_column = "id",
                           trace_annotation_id_column = "peptide_id")
  
  ## Remove traces with 0 standard deviation (can't be clustered)
  zerovar <- apply(getIntensityMatrix(traces), 1, var) > 0
  traces_zerovar <- subset(traces,
                           trace_subset_ids = names(zerovar[zerovar]))
  
  #' ## Remove single peptide genes
  traces_multiPep <- filterSinglePeptideHits(traces_zerovar)
  
  #' ## Estimate proteoform scores by correlation clustering
  traces_corr <- calculateGeneCorrMatrices(traces_multiPep)
  
  traces_clustered <- clusterPeptides(
    traces_corr,
    method = "average", plot = F, PDF=F,
    name=paste0("ProteoformClusters_interlab_",name))
  
  traces_clusteredInN <- cutClustersInNreal(traces_clustered, clusterN = 2,
                                            min_peptides_per_cluster = 2)
  
  traces_scored <- calculateProteoformScore(traces_clusteredInN)
  
  scores <- unique(traces_scored$trace_annotation[,.(proteoform_score,proteoform_score_z,proteoform_score_pval_adj,perturbed_protein), by=.(protein_id)])
  scores <- scores[!is.na(scores$proteoform_score)]
  
  
  
  library("ggExtra")
  scatter <- ggplot(scores, 
                    aes(x=proteoform_score,y=-log10(proteoform_score_pval_adj), 
                        color=perturbed_protein)) + 
    scale_color_manual(values=c("#ff7700", "green4"), name= "protein with proteoform") +
    geom_point(alpha=0.5) +
    geom_vline(xintercept = score_cutoff, colour='darkgrey', lty='dashed') +
    geom_hline(yintercept = -log10(adj_pval_cutoff), colour='darkgrey', lty='dashed') +
    theme_classic() +
    ylab("-log10 ( adj. pvalue )") +
    xlab("proteoform score") +
    theme(legend.position="bottom") 
  
  scatter_m <- ggMarginal(scatter, type = "histogram", groupFill = TRUE, groupColour = TRUE, 
             bins=75, position = "identity", size = 4)
  
  pdf(paste0("pseudo_volcano_hist_",name,".pdf"), height=5, width=5)
    print(scatter_m)
  dev.off()
  
  traces_proteoforms <- annotateTracesWithProteoforms(
    traces_scored, score_cutoff = score_cutoff, adj_pval_cutoff =  adj_pval_cutoff)
  
  return(traces_proteoforms)
}

performFDRbenchmark <- function(traces, name="traces", score_cutoff = 0.1, adj_pval_cutoff =  0.1){
  score_thresholds = seq(0,0.35,0.05)
  pval_thresholds = sort(unique(c(c(1 %o% 10^(-12:1)),seq(0,1,0.05))))
  
  n_row = nrow(expand.grid(score_thresholds,pval_thresholds))
  
  res <- data.table(score_thresholds=rep(0,n_row),
                    pval_thresholds=rep(0,n_row),
                    P = length(unique(traces$trace_annotation[perturbed_protein==TRUE]$protein_id)),
                    N = length(unique(traces$trace_annotation[perturbed_protein==FALSE]$protein_id)),
                    n_prot_proteoforms = 0,
                    TP = 0,
                    FP = 0,
                    N_prot_noMistakes = 0,
                    N_prot_withMistakes = 0,
                    TP_noMistakes = 0,
                    TP_withMistakes = 0
  )
  
  idx = 0
  for (i in seq(1,length(score_thresholds),1)) {
    for (j in seq(1,length(pval_thresholds),1)) {
      idx = idx+1
      res$score_thresholds[idx] <- score_thresholds[i]
      res$pval_thresholds[idx] <- pval_thresholds[j]
      
      traces_proteoforms <- annotateTracesWithProteoforms(
        traces, score_cutoff = score_thresholds[i],  adj_pval_cutoff = pval_thresholds[j])
      
      res$n_prot_proteoforms[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1)]$protein_id))
      res$TP[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein)]$protein_id))
      res$FP[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein==FALSE)]$protein_id))
      
      traces_proteoforms$trace_annotation[, n_proteoforms_per_perturbed_group := length(unique(proteoform_id)), by=c("protein_id","perturbed_peptide")]
      traces_proteoforms$trace_annotation[, n_perturbed_groups_per_proteoform := length(unique(perturbed_peptide)), by=c("protein_id","proteoform_id")]
      
      traces_proteoforms$trace_annotation[, correct_peptide := ((n_proteoforms_per_perturbed_group==1) & (n_perturbed_groups_per_proteoform==1))]
      
      traces_proteoforms$trace_annotation[, protein_without_mistakes := all(correct_peptide), by="protein_id"]
      
      res$N_prot_noMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[protein_without_mistakes==TRUE]$protein_id))
      res$N_prot_withMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[protein_without_mistakes==FALSE]$protein_id))
      
      res$TP_noMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein) & (protein_without_mistakes==TRUE)]$protein_id))
      res$TP_withMistakes[idx] <- length(unique(traces_proteoforms$trace_annotation[(n_proteoforms>1) & (perturbed_protein) & (protein_without_mistakes==FALSE)]$protein_id))
    }
  }
  
  res[,FDR := FP/(FP+TP)]
  res[,FN := P-TP]
  res[,F1 := TP/(TP+(0.5*(FP+FN)))]
  
  res[,TPR := TP/P]
  res[,FPR := FP/N]
  
  res[,percent_TP_perfect := TP_noMistakes/TP]
  res[,percent_perfect := N_prot_noMistakes/(P+N)]
  
  res$score_thresholds <- round(res$score_thresholds, digits = 2)
  
  pdf(paste0("FDR_benchmark_",name,".pdf"), width=6.5, height=5)
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkgrey")
  res$score_thresholds_fac = as.factor(res$score_thresholds)
  p <- ggplot(res, aes(x=pval_thresholds,y=FDR, group=score_thresholds, colour=score_thresholds_fac)) + 
    geom_point() + geom_line() + 
    xlim(0,0.5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = okabe) +
    xlab("adj. p-value threshold") + 
    labs(colour = "proteoform score \n threshold") +
    theme_classic()
  print(p)
  p <- ggplot(res, aes(x=pval_thresholds,y=FDR, group=score_thresholds, colour=score_thresholds_fac)) + 
    geom_point() + geom_line() + 
    xlim(0,0.17) +
    ylim(0,0.17) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = okabe) +
    xlab("adj. p-value threshold") + 
    labs(colour = "proteoform score \n threshold") +
    theme_classic()
  print(p)
  dev.off()
  
  res[,selected_pval := ifelse(pval_thresholds==adj_pval_cutoff,"TRUE","FALSE")]
  pval_col <- c("TRUE"="black","FALSE"="#00ff0000")
  
  pdf(paste0("ROC_percentPerfect_",name,"_pval_",adj_pval_cutoff,".pdf"), width=6.5, height=5)
  p <- ggplot(res[score_thresholds==score_cutoff], aes(x=FPR,y=TPR)) + 
    geom_line(colour="grey",alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, color='grey') +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1"),
                         breaks=seq(0,1,0.25)) +
    geom_point(aes(fill=percent_perfect, color=selected_pval,stroke=selected_pval), shape = 21, size = 3, stroke=1) +
    scale_color_manual(values = pval_col) +
    geom_abline(intercept = 0, slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    theme_classic()
  print(p)
  dev.off()
  
  return(res)
}

random_perm_scored <- performCCprofilerAnalysis(data=random_perm_data, name="random", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
plotProteoformPvalHist(random_perm_scored, name="random_perm_scored_pval_hist", PDF=T)
plotProteoformVolcano(random_perm_scored, name="random_perm_scored_volcano", PDF=T)

perm_1pep_scored <- performCCprofilerAnalysis(data=perm_1pep_data, name="perm_1pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_2pep_scored <- performCCprofilerAnalysis(data=perm_2pep_data, name="perm_2pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_025pep_scored <- performCCprofilerAnalysis(data=perm_025pep_data, name="perm_025pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
perm_050pep_scored <- performCCprofilerAnalysis(data=perm_050pep_data, name="perm_050pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)

fdr_bench_random_perm <- performFDRbenchmark(random_perm_scored, name="random_perm", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_1pep <- performFDRbenchmark(perm_1pep_scored, name="perm_1pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_2pep <- performFDRbenchmark(perm_2pep_scored, name="perm_2pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_025pep <- performFDRbenchmark(perm_025pep_scored, name="perm_025pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)
fdr_bench_perm_050pep <- performFDRbenchmark(perm_050pep_scored, name="perm_050pep", score_cutoff = 0.1, adj_pval_cutoff =  0.1)


res_copf <- rbind(fdr_bench_random_perm[,nPerm := "random"],
                  fdr_bench_perm_1pep[,nPerm := "1"],
                  fdr_bench_perm_2pep[,nPerm := "2"],
                  fdr_bench_perm_025pep[,nPerm := "25%"],
                  fdr_bench_perm_050pep[,nPerm := "50%"])

saveRDS(res_copf, "res_copf.rds")

#####################################
### PeCorA analysis #################
#####################################

install_github("jessegmeyerlab/PeCorA")
library(PeCorA)

performPecoraAnalysis <- function(data){
  input_data_pecora <- copy(data)
  
  setnames(input_data_pecora,c("peptide_id","protein_id","intensity","day"), c("Peptide","Protein","Normalized.Area","Condition"))
  input_data_pecora[,"Peptide.Modified.Sequence":=Peptide]
  
  ann <- unique(subset(input_data_pecora, select=c("filename","Condition")))
  setorderv(ann,c("Condition","filename"))
  ann[,BioReplicate:=c(1:7), by="Condition"]
  
  input_data_pecora <- merge(input_data_pecora,ann,by=c("filename","Condition"))
  
  input_data_pecora <- subset(input_data_pecora, select=c("Protein","Peptide","Peptide.Modified.Sequence","Condition","BioReplicate","Normalized.Area"))
  
  pecora_df <- setDF(input_data_pecora)
  
  scaled_peptides <- PeCorA_preprocessing(pecora_df,
                                          area_column_name=6,
                                          threshold_to_filter=min(pecora_df$Normalized.Area),
                                          control_name="day1")
  
  disagree_peptides <- PeCorA(scaled_peptides)
  
  pecora_res <- data.table(disagree_peptides)
  pecora_res[,"peptide_id":=gsub("_all","",peptide)]
  
  trace_annotation <- unique(subset(data, select=c("peptide_id","n_pep","n_perturbed_peptides",
                                                   "perturbed_protein","perturbed_peptide","red_fac")))
  
  pecora_res <- merge(pecora_res,trace_annotation,by="peptide_id")
  
  pecora_res[, adj_adj_pval := p.adjust(adj_pval, method = "BH")]
  
  return(pecora_res)
}


pecoraFDRbenchmark <- function(data, name="pecora", pval_col = "adj_pval", pval_cutoff =  0.05){
  
  pecora_res <- copy(data)
  
  all_pos_prot <- length(unique(pecora_res[(perturbed_protein==TRUE)]$protein))
  all_neg_prot <- length(unique(pecora_res[(perturbed_protein==FALSE)]$protein))
  
  all_pos_pep <- length(unique(pecora_res[perturbed_peptide==TRUE]$peptide))
  all_neg_pep <- length(unique(pecora_res[perturbed_peptide==FALSE]$peptide))
  
  score_thresholds = sort(unique(c(c(1 %o% 10^(-12:1)),seq(0,1,0.05))))
  
  pecora_res_df <- data.table(score_thresholds=score_thresholds,
                              P = all_pos_prot,
                              N = all_neg_prot,
                              n_prot_proteoforms = 0,
                              TP = 0,
                              FP = 0,
                              N_prot_noMistakes = 0,
                              N_prot_withMistakes = 0,
                              TP_noMistakes = 0,
                              TP_withMistakes = 0
  )
  
  for (i in seq(1,length(score_thresholds),1)) {
    pecora_res_sig <- subset(pecora_res, get(pval_col)<=score_thresholds[i])
    
    pecora_res_df$n_prot_proteoforms[i] <- length(unique(pecora_res_sig$peptide))
    
    TP_pep = length(unique(pecora_res_sig[perturbed_peptide==TRUE]$peptide))
    FP_pep = length(unique(pecora_res_sig[perturbed_peptide==FALSE]$peptide))
    TP_prot = length(unique(pecora_res_sig[perturbed_protein==TRUE]$protein))
    FP_prot = length(unique(pecora_res_sig[perturbed_protein==FALSE]$protein))
    
    pecora_res_df$TP[i] <- TP_prot
    pecora_res_df$FP[i] <- FP_prot
    
    pecora_res[, TP_perturbed := (perturbed_peptide==TRUE) & (get(pval_col)<=score_thresholds[i]), by= "peptide"]
    pecora_res[, FP_perturbed := (perturbed_peptide==FALSE) & (get(pval_col)<=score_thresholds[i]), by= "peptide"]
    
    pecora_res[, correct_assignment := (((perturbed_peptide==TRUE) & (get(pval_col)<=score_thresholds[i])) | ((perturbed_peptide==FALSE) & (get(pval_col)>score_thresholds[i]))) , by= "peptide"]
    
    pecora_res[, perfect_prot := all(correct_assignment), by="protein"] 
    pecora_res[, min_adj_pval := min(get(pval_col)), by="protein"]
    
    pecors_score <- unique(subset(pecora_res, select=c("protein","perturbed_protein","perfect_prot","min_adj_pval")))
    pecora_res_df$N_prot_noMistakes[i] <- sum(pecors_score$perfect_prot)
    pecora_res_df$N_prot_withMistakes[i] <- nrow(pecors_score[perfect_prot==FALSE])
    
    pecora_res_df$TP_noMistakes[i] <- sum(pecors_score[(perturbed_protein==TRUE) & (min_adj_pval<=score_thresholds[i])]$perfect_prot)
    pecora_res_df$TP_withMistakes[i] <- nrow(pecors_score[(perturbed_protein==TRUE) & (min_adj_pval<=score_thresholds[i]) & (perfect_prot==FALSE)])
    
  }
  
  pecora_res_df[,TPR:=TP/P]
  pecora_res_df[,FPR:=FP/N]
  pecora_res_df[,percent_TP_perfect := TP_noMistakes/TP]
  pecora_res_df[,percent_perfect := N_prot_noMistakes/(P+N)]
  
  pecora_res_df[,FDR := FP/(FP+TP)]
  pecora_res_df[,FN := P-TP]
  pecora_res_df[,F1 := TP/(TP+(0.5*(FP+FN)))]
  
  pecora_res_df[,pval_thresholds:=score_thresholds]
  pecora_res_df[,selected_pval := ifelse(pval_thresholds==pval_cutoff,TRUE,FALSE)]
  pval_col <- c("TRUE"="black","FALSE"="#00ff0000")
  
  pdf(paste0("pecora_ROC_percentPerfect_",name,"_pval_",pval_cutoff,".pdf"), width=6.5, height=5)
  p <- ggplot(pecora_res_df, aes(x=FPR,y=TPR)) + 
    geom_line(colour="grey",alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, color='grey') +
    scale_fill_gradientn(limits = c(0,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1"),
                         breaks=seq(0,1,0.25)) +
    geom_point(aes(fill=percent_perfect, color=selected_pval,stroke=selected_pval), shape = 21, size = 3, stroke=1) +
    scale_color_manual(values = pval_col) +
    geom_abline(intercept = 0, slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    theme_classic()
  print(p)
  dev.off()
  
  return(pecora_res_df)
}


random_perm_pecora <- performPecoraAnalysis(data=random_perm_data)
fdr_bench_random_perm_pecora <- pecoraFDRbenchmark(data=random_perm_pecora, name="random_perm", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_random_perm_pecora_adj <- pecoraFDRbenchmark(data=random_perm_pecora, name="random_perm_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_1pep_pecora <- performPecoraAnalysis(data=perm_1pep_data)
fdr_bench_perm_1pep_pecora <- pecoraFDRbenchmark(data=perm_1pep_pecora, name="perm_1pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_1pep_pecora_adj <- pecoraFDRbenchmark(data=perm_1pep_pecora, name="perm_1pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_2pep_pecora <- performPecoraAnalysis(data=perm_2pep_data)
fdr_bench_perm_2pep_pecora <- pecoraFDRbenchmark(data=perm_2pep_pecora, name="perm_2pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_2pep_pecora_adj <- pecoraFDRbenchmark(data=perm_2pep_pecora, name="perm_2pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_025pep_pecora <- performPecoraAnalysis(data=perm_025pep_data)
fdr_bench_perm_025pep_pecora <- pecoraFDRbenchmark(data=perm_025pep_pecora, name="perm_025pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_025pep_pecora_adj <- pecoraFDRbenchmark(data=perm_025pep_pecora, name="perm_025pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

perm_050pep_pecora <- performPecoraAnalysis(data=perm_050pep_data)
fdr_bench_perm_050pep_pecora <- pecoraFDRbenchmark(data=perm_050pep_pecora, name="perm_050pep", pval_col = "adj_pval", pval_cutoff =  0.1)
fdr_bench_perm_050pep_pecora_adj <- pecoraFDRbenchmark(data=perm_050pep_pecora, name="perm_050pep_adj", pval_col = "adj_adj_pval", pval_cutoff =  0.1)

res_pecora <- rbind(fdr_bench_random_perm_pecora[,c("nPerm","pval_col") := list("random","adj_pval")],
                    fdr_bench_perm_1pep_pecora[,c("nPerm","pval_col") := list("1","adj_pval")],
                    fdr_bench_perm_2pep_pecora[,c("nPerm","pval_col") := list("2","adj_pval")],
                    fdr_bench_perm_025pep_pecora[,c("nPerm","pval_col") := list("25%","adj_pval")],
                    fdr_bench_perm_050pep_pecora[,c("nPerm","pval_col") := list("50%","adj_pval")],
                    fdr_bench_random_perm_pecora_adj[,c("nPerm","pval_col") := list("random","adj_adj_pval")],
                    fdr_bench_perm_1pep_pecora_adj[,c("nPerm","pval_col") := list("1","adj_adj_pval")],
                    fdr_bench_perm_2pep_pecora_adj[,c("nPerm","pval_col") := list("2","adj_adj_pval")],
                    fdr_bench_perm_025pep_pecora_adj[,c("nPerm","pval_col") := list("25%","adj_adj_pval")],
                    fdr_bench_perm_050pep_pecora_adj[,c("nPerm","pval_col") := list("50%","adj_adj_pval")])

saveRDS(res_pecora, "res_pecora.rds")


#####################################
### Compare PeCorA & CCprofiler #####
#####################################

performBenchmarkComparison <- function(ccprofiler_res, pecora_res, name="copf_pecora_benchmark", ccprofiler_score_cutoff = 0.1){
  
  res <- copy(ccprofiler_res)
  pecora_res_df <- copy(pecora_res)
  
  res[,method:="COPF"]
  pecora_res_df[,method:="PeCorA"]
  
  res[,score_thresholds_fac:=NULL]
  pecora_res_df <- setcolorder(pecora_res_df, names(res))
  pecora_res_df[,method:=ifelse(pval_col=="adj_adj_pval", "PeCorA_adj", "PeCorA")]
  pecora_res_df[,pval_col:=NULL]
  
  res_filtered_01 <- res[score_thresholds == ccprofiler_score_cutoff]
  res_filtered_01[,method:="COPF_0.1"]
  
  res_filtered_00 <- res[score_thresholds == 0]
  res_filtered_00[,method:="COPF"]
  
  res_combi <- rbind(res_filtered_01,res_filtered_00,pecora_res_df)
  
  cols <- c("PeCorA" = "blue4", "COPF" = "green4", "PeCorA_adj" = "blue1", "COPF_0.1" = "green2")
  
  pdf(paste0("ROC_benchmark_ccprofiler_pecora_",name,"_score_cutoff_",ccprofiler_score_cutoff,".pdf"), height=4, width=9)
  p <- ggplot(res_combi, aes(x=FPR, y= TPR, shape=method, color=percent_TP_perfect)) + 
    geom_abline(intercept = 0, slope = 1, color='grey') +
    geom_line(alpha=0.5) +
    geom_point() +
    scale_color_gradientn(limits = c(0,1),
                          colours=c("navyblue", "darkmagenta", "darkorange1"),
                          breaks=seq(0,1,0.25)) +
    geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
    xlim(0,1) +
    ylim(0,1) +
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom")
  print(p)
  p <- ggplot(res_combi, aes(x=FPR, y= TPR, shape=method, color=percent_TP_perfect)) + 
    geom_line(alpha=0.5) +
    geom_point() +
    scale_color_gradientn(limits = c(0,1),
                          colours=c("navyblue", "darkmagenta", "darkorange1"),
                          breaks=seq(0,1,0.25)) +
    geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
    xlim(0,0.1) +
    ylim(0,1) +
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom")
  print(p)
  dev.off()
  
  pdf(paste0("FDR_benchmark_ccprofiler_pecora_",name,"_score_cutoff_",ccprofiler_score_cutoff,".pdf"), height=4, width=9)
  p <- ggplot(res_combi, aes(x=pval_thresholds, y= FDR, colour=method, shape=method, group=method)) + 
    geom_abline(intercept = 0, slope = 1, color='grey') +
    xlab("adj. p-value") +
    ylab("empirical FDR") +
    geom_line()+
    geom_point()+
    scale_color_manual(values = cols) +
    xlim(0,1)+
    ylim(0,1)+
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom")
  print(p)
  p <- ggplot(res_combi, aes(x=pval_thresholds, y= FDR, colour=method, shape=method, group=method)) + 
    geom_abline(intercept = 0, slope = 1, color='grey') +
    xlab("adj. p-value") +
    ylab("empirical FDR") +
    geom_line()+
    geom_point()+
    scale_color_manual(values = cols) +
    geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
    xlim(0,1)+
    ylim(0,1)+
    facet_wrap(. ~ nPerm) +
    theme_classic() +
    theme(legend.position="bottom")
  print(p)
  dev.off()
  
  return(res_combi)
}


res_combi <- performBenchmarkComparison(ccprofiler_res=res_copf, pecora_res=res_pecora, name="copf_pecora_benchmark", ccprofiler_score_cutoff = 0.1)

saveRDS(res_combi, "res_combi.rds")


paper_sel <- c("1","2","50%")
res_combi_paper <- performBenchmarkComparison(ccprofiler_res=res_copf[nPerm %in% paper_sel], 
                                        pecora_res=res_pecora[nPerm %in% paper_sel][pval_col=="adj_pval"], 
                                        name="PAPER", ccprofiler_score_cutoff = 0.1)

res_combi_paper <- performBenchmarkComparison(ccprofiler_res=res_copf[nPerm %in% paper_sel], 
                                              pecora_res=res_pecora[nPerm %in% paper_sel], 
                                              name="SUPPLEMENT", ccprofiler_score_cutoff = 0.1)


##############
##############
##############

pdf("ROC_benchmark_ccprofiler_pecora_nPermComp.pdf", width=6, height=4)
ggplot(res_combi[nPerm != 'random'], aes(x=FPR, y= TPR, color=nPerm, group=nPerm)) + 
  geom_abline(intercept = 0, slope = 1, color='grey') +
  geom_line(alpha=0.5) +
  geom_point() +
  scale_color_manual(values =c("#A9C9A4","#699864","#308014","#2F4F2F")) +
  #geom_point(data=res_combi[(selected_pval==TRUE)], color="red", shape=1, size=3) +
  xlim(0,1) +
  ylim(0,1) +
  facet_wrap(. ~ method) +
  theme_classic()
dev.off()

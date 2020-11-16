generateMixtureTraces <- function(traces, trueInteractions=NULL, n_proteins = 2, n_mixtureTraces = 500, peptide_ratio = 0.5, min_peptide_count = 2, seed = 123, n_iterations=10) {
  unique_proteins <- unique(traces$trace_annotation$protein_id)
  max_mixtureTraces <- choose(length(unique_proteins), n_proteins)
  if (max_mixtureTraces < n_mixtureTraces) {
    message(paste0("Can only generate maximally ", length(max_mixtureTraces), " mixture traces."))
  }
  n_mixtureTraces <- min(n_mixtureTraces, ncol(max_mixtureTraces))
  
  set.seed(seed)
  if (is.null(trueInteractions)) {
    protein_combinations_noInteraction_max <- replicate(n_mixtureTraces, sample(unique_proteins, n_proteins, replace = FALSE))
  } else {
    #trueInteractions[,idx:=.I]
    #trueInteractions[,id:= paste(sort(c(a,b)),collapse = ";"), by=idx]
    protein_combinations_noInteraction <- matrix()
    it <- 0
    while((ncol(protein_combinations_noInteraction) < n_mixtureTraces) & (it < n_iterations)) {
      it <- it+1
      protein_combinations <- replicate(n_mixtureTraces*2, sample(unique_proteins, n_proteins, replace = FALSE))
      protein_combinations_id <- lapply(c(1:ncol(protein_combinations)),
                                        function(j){
                                          binary_comb <- combn(protein_combinations[,j], m = 2)
                                          lapply(c(1:ncol(binary_comb)),function(x){paste(sort(binary_comb[,x]),collapse = ";")})
                                        })
      noInteraction_idxs_TF <- lapply(protein_combinations_id,function(x){any(unlist(x) %in% trueInteractions$id)})
      
      noInteraction_idxs <- which(noInteraction_idxs_TF=="FALSE")
      protein_combinations_noInteraction <- protein_combinations[,noInteraction_idxs]
    }
    if(ncol(protein_combinations_noInteraction) < n_mixtureTraces) {
      message(paste0("Maximum number of iterations (",it,") reached. Only ",
                     ncol(protein_combinations_noInteraction)," mixtur proteins could be generated."))
      protein_combinations_noInteraction_max <- protein_combinations_noInteraction
    } else {
      protein_combinations_noInteraction_max <- protein_combinations_noInteraction[,1:n_mixtureTraces]
    }
  }
  
  set.seed(seed)
  traces_subsets <- lapply(c(1:ncol(protein_combinations_noInteraction_max)), function(n){
    proteins <- protein_combinations_noInteraction_max[,n]
    proteins <- sort(proteins)
    peptides <- lapply(proteins, function(p){
      traces_sub <- subset(traces, p, trace_subset_type = "protein_id")
      all_peptides <- unique(traces_sub$trace_annotation$id)
      max_peptide_count <- ceiling(length(all_peptides)*peptide_ratio)
      n_peptides_to_select <- max(min_peptide_count,max_peptide_count)
      peptide_subset <- sample(all_peptides, n_peptides_to_select, replace = FALSE)
      return(peptide_subset)
    })
    traces_subset <- subset(traces, unlist(peptides))
    prot_names <- paste(proteins, collapse = "_")
    traces_subset$trace_annotation[,id := paste(paste(prot_names,id, sep = "_"), protein_id, sep="_")]
    traces_subset$traces[,id := traces_subset$trace_annotation$id]
    traces_subset$trace_annotation[,protein_id := prot_names]
    traces_subset$trace_annotation[,n_mixed_proteins := n_proteins]
    traces_subset$trace_annotation[,is_mixed:=TRUE]
    if ("genomic_coord" %in% names(traces_subset)) {
      traces_subset$genomic_coord <- NULL
    }
    return(traces_subset)
  })
  mixed_traces <- list(traces = do.call(rbind,lapply(traces_subsets, `[[`, "traces")),
                       trace_type = traces$trace_type,
                       trace_annotation = do.call(rbind,lapply(traces_subsets, `[[`, "trace_annotation")),
                       fraction_annotation = traces$fraction_annotation)
  class(mixed_traces) <- "traces"
  #.tracesTest(mixed_traces)
  return(mixed_traces)
}

evaluateProteoformClusteringOfMixtureTraces <- function(traces, score_cutoff){
  all_proteins <- unique(traces$trace_annotation$protein_id)
  mistakes <- lapply(all_proteins, function(x){
    sub <- subset(traces$trace_annotation, protein_id==x)
    
    sub_noOutliers <- subset(sub, cluster != 100)
    
    n_outlierPeptides <- nrow(subset(sub, cluster == 100))
    
    unique_proteins <- unique(sub$Entry_name)
    n_proteins <- length(unique_proteins)
    
    unique_proteins_noOutliers <- unique(sub_noOutliers$Entry_name)
    n_proteins_noOutliers <- length(unique_proteins_noOutliers)
    
    sub_noOutliers[, proteoform_id := ifelse(proteoform_score >= score_cutoff, paste(protein_id, cluster, sep ='_'), protein_id)]
    
    unique_proteoforms <- unique(sub_noOutliers$proteoform_id)
    n_proteoforms <- length(unique_proteoforms)
    
    n_mistakes <- sum(unlist(lapply(unique_proteoforms,function(p){
      subSub <- subset(sub_noOutliers,proteoform_id==p)
      nrow(subSub)-max(table(subSub$Entry_name))
    })))
    return(data.table(protein=x, n_outlierPeptides=n_outlierPeptides, n_proteins=n_proteins, n_proteins_noOutliers=n_proteins_noOutliers, n_proteoforms=n_proteoforms, n_mistakes=n_mistakes))
  })
  mistakes <- do.call(rbind,mistakes)
  return(mistakes)
}

source("/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")

#' ## Load data
tracesList_location <- readRDS("tracesList_location.rds")
design_matrix <- readRDS("design_matrix.rda")

removeOutlierPeptides <- function(traces){
  outlier_ids <- traces$trace_annotation[(n_proteoforms==2) & (cluster==100)]$id
  all_ids <- traces$trace_annotation$id
  non_outlier_ids <- all_ids[! all_ids %in% outlier_ids]
  sub_traces <- subset(traces, trace_subset_ids = non_outlier_ids)
  return(sub_traces)
}

filterProteoformsByThreshold <- function(traces, score_threshold=0.5){
  non_proteoform_proteins <- traces$trace_annotation[(n_proteoforms==2) & (proteoform_score<score_threshold)]$protein_id
  traces$trace_annotation[,n_proteoforms := ifelse(protein_id %in% non_proteoform_proteins, 0, n_proteoforms)]
  traces$trace_annotation[,proteoform_id := ifelse(n_proteoforms == 0, protein_id, proteoform_id)]
  return(traces)
}

#' ## Re-annotate tarces based on the selected  score threshold of 0.25
traces_list_pepClusters <- lapply(tracesList_location, function(x){
  n <- removeOutlierPeptides(x)
  m <- filterProteoformsByThreshold(n, score_threshold=0.25)
  return(m)
  })
class(traces_list_pepClusters) <- "tracesList"

#' ## Integrate traces across conditions
traces_sum_pepClusters <- integrateTraceIntensities(traces_list_pepClusters,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")

saveRDS(traces_sum_pepClusters,"traces_sum_pepClusters.rds")

#' ## Proteoform-specific feature finding
proteinFeatures  <- findProteinFeatures(traces=traces_sum_pepClusters, 
                                        corr_cutoff=0.9,
                                        window_size=7,
                                        parallelized=T,
                                        n_cores=3,
                                        collapse_method="apex_only",
                                        perturb_cutoff= "5%",
                                        rt_height=1,
                                        smoothing_length=7,
                                        useRandomDecoyModel=TRUE,
                                        quantLevel = "protein_id")

saveRDS(proteinFeatures, "proteinFeatures.rda")

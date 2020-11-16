source("/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")

#' ## Load data
traces_list_pepClusters <- readRDS("tracesList_location.rds")
traces_sum_pepClusters <- readRDS("traces_sum_pepClusters.rds")
proteinFeatures <- readRDS("proteinFeatures.rda")
design_matrix <- readRDS("design_matrix.rda")
calibrationFunctions <- readRDS("calibration.rds")

#' ## Resolve protein features by clusters
proteoformFeaturesResolved <- resolveProteoformSpecificFeatures(
    features=proteinFeatures,
    traces=traces_sum_pepClusters,
    minProteoformIntensityRatio=0.1,
    perturb_cutoff="5%")

saveRDS(proteoformFeaturesResolved,"proteoformFeaturesResolved.rds")

#' ## Score protein features
filteredDataProteoformResolved <- scoreFeatures(
  proteoformFeaturesResolved,
  FDR=0.05, PDF=T,
  name=paste0("qvalueStats_proteoformFeaturesResolved"))

saveRDS(filteredDataProteoformResolved,"filteredDataProteoformResolved.rds")


summarizeFeatures(filteredDataProteoformResolved, PDF=T, name="filteredDataProteoformResolved_summary")



##############

proteins_proteoforms_025 <- fread("ProteinTables/proteoforms_score_0.25.txt", header = FALSE, col.names = 'ID')

filteredDataProteoformResolved_proteoforms <- subset(filteredDataProteoformResolved, protein_id %in% proteins_proteoforms_025$ID)

length(unique(filteredDataProteoformResolved_proteoforms$protein_id))

filteredDataProteoformResolved_proteoforms[,n_unique_proteoforms := length(unique(proteoform_ids)), by="protein_id"]

non_SEC_proteoforms <- unique(filteredDataProteoformResolved_proteoforms[n_unique_proteoforms==1]$protein_id)
assembly_specific_proteoforms <- unique(filteredDataProteoformResolved_proteoforms[n_unique_proteoforms>1]$protein_id)

write.table(assembly_specific_proteoforms,paste0("ProteinTables/assembly_specific_proteoforms_score_0.25.txt"), sep="\t", quote = F, col.names = F, row.names = F)

source('../CCprofilerAnalysis/thesis/traces_plotting.R')

pdf("assembly_specific_proteoforms.pdf",height=4, width=6)
for (p in assembly_specific_proteoforms){
  plotSub <- subset(
    traces_list_pepClusters, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()

pdf("non_SEC_proteoforms.pdf",height=4, width=6)
for (p in non_SEC_proteoforms){
  plotSub <- subset(
    traces_list_pepClusters, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()


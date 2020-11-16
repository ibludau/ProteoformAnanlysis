source("/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")

#' ## Load data
traces_list_reassignedProteoforms <- readRDS("tracesList_location.rds")
reassignedProteoformFeatures <- readRDS("filteredDataProteoformResolved.rds")
design_matrix <- readRDS("design_matrix.rda")

#' ## Extract feature values
proteoform_featureVals <- extractFeatureVals(traces = traces_list_reassignedProteoforms,
                                             features = reassignedProteoformFeatures,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F)
saveRDS(proteoform_featureVals, "proteoform_featureVals.rda")

#' ## Fill feature values
proteoform_featureValsFilled <- fillFeatureVals(featureVals = proteoform_featureVals,
                                                tracesList = traces_list_reassignedProteoforms,
                                                design_matrix = design_matrix)
saveRDS(proteoform_featureValsFilled, "proteoform_featureValsFilled.rda")

#' ## Perform peptide-level differential expression testing for all features
proteoform_DiffExprPep <- testDifferentialExpression(featureVals = proteoform_featureValsFilled,
                                                     compare_between = "Condition",
                                                     level = "peptide",
                                                     measuredOnly = FALSE)

saveRDS(proteoform_DiffExprPep, "proteoform_DiffExprPep.rds")

#' ## Aggregate differential expression results to the proteoform level
proteoform_DiffExprProteoform <- aggregatePeptideTestsToProteoform(proteoform_DiffExprPep)

#' ## Aggregate differential expression results to the protein level
proteoform_DiffExprProtein <- aggregatePeptideTests(proteoform_DiffExprPep)

#' Change sign of fold-change to accomodate the control as the reference (reverse to how test was calculated)
proteoform_DiffExprPep[,log2FC:=-log2FC]
proteoform_DiffExprPep[,global_log2FC_imp:=-global_log2FC_imp]
saveRDS(proteoform_DiffExprPep, "proteoform_DiffExprPep.rda")

proteoform_DiffExprProteoform[,medianLog2FC:=-medianLog2FC]
proteoform_DiffExprProteoform[,global_medianLog2FC:=-global_medianLog2FC]
proteoform_DiffExprProteoform[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(proteoform_DiffExprProteoform, "proteoform_DiffExprProteoform.rda")

proteoform_DiffExprProtein[,medianLog2FC:=-medianLog2FC]
proteoform_DiffExprProtein[,global_medianLog2FC:=-global_medianLog2FC]
proteoform_DiffExprProtein[,global_medianLog2FC_imp:=-global_medianLog2FC_imp]
saveRDS(proteoform_DiffExprProtein, "proteoform_DiffExprProtein.rda")

#' Make volcano plots
#' General volcano plots
plotVolcano(proteoform_DiffExprPep, PDF = T, name = "proteoform_DiffExprPep")
plotVolcano(proteoform_DiffExprProteoform, PDF = T, name = "proteoform_DiffExprProteoform")
plotVolcano(proteoform_DiffExprProtein, PDF = T, name = "proteoform_DiffExprProtein")

#' Volcanoplts highlighting different proteoforms
library(ggrepel)
plotVolcano(proteoform_DiffExprProtein, highlight=c("Q7Z3B4"), PDF = T, name = "prot_DiffExprProt_NUP54_Q7Z3B4")

#' Volcanoplts highlighting different proteoforms >> example
library(ggrepel)
plotVolcano(proteoform_DiffExprProteoform, highlight=c("Q8WWM7"), PDF = T, name = "prot_DiffExprProteoform_Q8WWM7")
proteoform_DiffExprProteoform[feature_id=="Q8WWM7"][pBHadj <= 0.05][abs(medianLog2FC)>1]


diff <- proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]
nrow(proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)])
length(unique(proteoform_DiffExprProtein[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]$feature_id))

diff <- proteoform_DiffExprProteoform[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]
nrow(proteoform_DiffExprProteoform[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)])
length(unique(proteoform_DiffExprProteoform[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]$feature_id))

#####

proteins_proteoforms_025 <- fread("ProteinTables/proteoforms_score_0.25.txt", header = FALSE, col.names = 'ID')

proteoform_DiffExprProteoform_proteoforms <- subset(proteoform_DiffExprProteoform, feature_id %in% proteins_proteoforms_025$ID)
proteoform_DiffExprProteoform_proteoforms_diff <- proteoform_DiffExprProteoform_proteoforms[(abs(medianLog2FC)>=1) & (pBHadj <= 0.05)]

proteoform_proteins_diff <- unique(proteoform_DiffExprProteoform_proteoforms_diff$feature_id)
length(proteoform_proteins_diff)

write.table(proteoform_proteins_diff,paste0("ProteinTables/cellcycle_regulated_proteoforms_score_0.25.txt"), sep="\t", quote = F, col.names = F, row.names = F)

source('../CCprofilerAnalysis/thesis/PaperAnalysis/traces_plotting.R')

traces_list_pepClusters <- readRDS("tracesList_location.rds")

pdf("proteoform_proteins_diff.pdf",height=4, width=6)
for (p in proteoform_proteins_diff){
  plotSub <- subset(
    traces_list_pepClusters, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()


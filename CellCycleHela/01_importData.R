source("/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")

#' ## Read input data and annotate with RNAseq information
#' Read TRIC output file:
data <- fread("../InputData/20170906213658025-1332538/E1709051521_feature_alignment.tsv")
data[, filename := gsub(
  ".mzXML\\.gz",
  "", filename)]
data[, filename := gsub(
  "/scratch/.*.tmpdir/",
  "", filename)]
setkey(data, "filename")

#' Read SEC fraction annotation table:
ann <- fread("../InputData/HeLaCCL2_SEC_annotation_full.csv")
ann[, filename := file]
ann[, sample := paste0(condition, replicate)]

#' Craete tracesList object:
traces_list <- importMultipleCondiionsFromOpenSWATH(data=data,
                                                    annotation_table=ann,
                                                    rm_requantified=FALSE,
                                                    rm_decoys = TRUE,
                                                    rm_nonProteotypic = TRUE,
                                                    MS1Quant=FALSE,
                                                    proteogenomicsWF=FALSE,
                                                    verbose=F)


#' Clean memory:
rm(data)
gc()

#' Inspect traces list:
summary(traces_list)

#' Annotate traces with biomaRt info:
traces_list <- annotateTraces(traces_list,
                              exampleTraceAnnotation,
                              traces_id_column = "protein_id",
                              trace_annotation_id_column = "Entry")

#' ## Molecular weight calibration
#' Perform molecular weight calibration:
calibration = calibrateMW(exampleCalibrationTable,
                          PDF=TRUE,
                          plot=TRUE)
#' Annotate fractions in tares list with according molecular weights:
traces_list <- annotateMolecularWeight(traces_list, calibration)

#' Update traces with additional metrics for each fraction:
traces_list <- updateTraces(traces_list)
#' Inspect traces list:
summary(traces_list)

#' ## Create design matrix
design_matrix <- data.table(Sample_name = c("Mitosis1","Mitosis2","Mitosis3","Interphase1","Interphase2","Interphase3"),
                            Condition = c("Mitosis","Mitosis","Mitosis","Interphase","Interphase","Interphase"),
                            Replicate = c(1,2,3,1,2,3))

#' ## Save objects
saveRDS(calibration,"calibration.rds")
saveRDS(design_matrix, "design_matrix.rda")
saveRDS(traces_list, "pepTracesRaw.rda")

#' ## QC plots
test_proteins <- c("P06493")
pepTest <- subset(traces_list, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = T, PDF = TRUE,
     name = paste0("RAW_PeptideTraces_CDK1_","P06493"), isoformAnnotation = F, plot = TRUE, highlight = NULL, highlight_col = NULL, colour_by = "Entry_name")

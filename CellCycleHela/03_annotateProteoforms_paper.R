source("/Users/isabell/Desktop/projects/ProteoformProject/CCprofilerAnalysis/thesis/CellCycle_header.R")

####################################
## Load data #######################
####################################

pepTracesList_filtered <- readRDS("pepTracesList_filtered.rda")
design_matrix <- readRDS("design_matrix.rda")

####################################
## Preprocess traces ###############
####################################

pepTracesC <- integrateReplicates(pepTracesList_filtered,
                                  design_matrix = design_matrix,
                                  integrate_within = "Condition",
                                  repPepCorr_cutoff = 0)

pepTracesC <- combineTracesMutiCond(pepTracesC)

####################################
## Clustering ######################
####################################

#' ## Estimate proteoform scores by correlation clustering
traces_corr <- calculateGeneCorrMatrices(pepTracesC)

traces_clustered <- clusterPeptides(
  traces_corr,
  method = "average", plot = T, PDF=T,
  name="ProteoformClusters")

traces_clusteredInN <- cutClustersInNreal(traces_clustered,
                                          clusterN = 2,
                                          min_peptides_per_cluster = 2)

traces_scored <- calculateProteoformScore(traces_clusteredInN)

traces_scored$trace_annotation[, proteoform_score := ifelse(is.na(proteoform_score), 0, proteoform_score)]
traces_scored$trace_annotation[, proteoform_score_pval_adj := ifelse(is.na(proteoform_score_pval_adj), 1, proteoform_score_pval_adj)]

plotProteoformPvalHist(traces_scored, name="traces_scored_pval_hist", PDF=T)

plotProteoformVolcano(traces_scored, name="traces_scored_volcano", PDF=T)

source("../CCprofilerAnalysis/thesis/PaperAnalysis/CellCycleHela/plotProteoformVolcanoLines.R")
plotProteoformVolcanoLines(traces_scored, name="traces_scored_volcano_lines", score_cutoff = 0.1, adj_pval_cutoff = 0.1, PDF=T)

#' ## Annotate traces with proteoform_ids. 
# Selected zero cutoff for downstream evaluation.
traces_proteoforms <- annotateTracesWithProteoforms(
  traces_scored, score_cutoff = 0, adj_pval_cutoff =  1)

#' ## Evaluate sequence location of the determined proteoforms
traces_location <- evaluateProteoformLocation(
  traces_proteoforms)

####################################
## Load external phospho data ######
####################################

pdata <- fread("/Users/isabell/Downloads/Supplementary Table 4-Table 1.tsv", quote = F)
setnames(pdata, colnames(pdata), as.character(unlist(pdata[2,])))
pdata <- pdata[-c(1,2),]

pdata_sig <- subset(pdata, `Mitosis/Interphase` == TRUE)
pdata_sig[,log2R := as.numeric(gsub(",",".",`Log2_ratios Mitosis/Interphase`))]
pdata_sig <- subset(pdata_sig, abs(log2R) > 0.5)
pdata_sig[, protein_id := pdata_sig$`Protein Group`]

significantly_regulated_proteins <- unique(pdata_sig$protein_id)
significantly_regulated_peptides <- unique(pdata_sig$Sequence)

#' ## Important functions
library("seqinr")
library("Biostrings")
getPepStartSite <- function(pep, prot){
  if (prot %in% names(fasta)){
    return(words.pos(pep, toString(fasta[prot]))[1])
  } else {
    print(prot)
    return(NA)
  }
}

fasta <- readAAStringSet("../InputData/human.fasta")
names(fasta) = gsub(".*\\|(.*?)\\|.*", "\\1", names(fasta))

#' ## Annotate start and end sites:
pdata_sig[,PeptidePositionStart := getPepStartSite(Sequence, protein_id), by=c("Sequence", "protein_id")]

pdata_sig[, ptm_site := gsub("\"","",`Phospho_Modification with PhosphoRS probability`)]
pdata_sig[, ptm_site := gsub(" ","",ptm_site)]

getPsite <- function(s){
  x=as.data.table(t(as.data.table(lapply(unlist(strsplit(s, split=";")), function(x){unlist(strsplit(x,split=":"))}))))
  colnames(x)=c("site","probability")
  x[,probability := as.numeric(probability)]
  x_sig <- subset(x, probability > 90)
  if(nrow(x_sig) > 0){
    p_site = as.numeric(gsub('.*\\(([0-9]+)\\)','\\1',x_sig$site))
    p_site = paste(p_site, collapse = ";")
    return(p_site)
  } else {
    return("NA")
  }
}

pdata_sig[, phosphosite := lapply(ptm_site, getPsite), by=c("Sequence","Modification","protein_id")]
pdata_sig[, phosphosite_protSite := PeptidePositionStart + as.numeric(phosphosite) -1, by=c("phosphosite","Sequence","Modification","protein_id")]

phosphosite_annotation_table <- subset(pdata_sig, select=c("protein_id","phosphosite_protSite"))

saveRDS(phosphosite_annotation_table, 'phosphosite_annotation_table.rds')
write.table(phosphosite_annotation_table, 'phosphosite_annotation_table.tsv', quote = F, sep = '\t', row.names = F, col.names = T)


# for app visualization
pdata[, protein_id := pdata$`Protein Group`]
pdata[,PeptidePositionStart := getPepStartSite(Sequence, protein_id), by=c("Sequence", "protein_id")]
pdata[, ptm_site := gsub("\"","",`Phospho_Modification with PhosphoRS probability`)]
pdata[, ptm_site := gsub(" ","",ptm_site)]
pdata[, phosphosite := lapply(ptm_site, getPsite), by=c("Sequence","Modification","protein_id")]
pdata[, phosphosite_protSite := PeptidePositionStart + as.numeric(phosphosite) -1, by=c("phosphosite","Sequence","Modification","protein_id")]
all_phosphosite_annotation_table <- subset(pdata, select=c("protein_id","phosphosite_protSite"))
saveRDS(all_phosphosite_annotation_table, 'all_phosphosite_annotation_table.rds')
write.table(all_phosphosite_annotation_table, 'all_phosphosite_annotation_table.tsv', quote = F, sep = '\t', row.names = F, col.names = T)


### Annotate proteoform annotation with phosphosite regulation info
proteoform_annotation <- traces_location$trace_annotation

proteoform_annotation[,phophoProtein := ifelse(protein_id %in% significantly_regulated_proteins, 1, 0)]


getPhosphositeAnn <- function(prot,start,end, ann_table){
  ann <- subset(ann_table, (protein_id==prot) & (phosphosite_protSite >= start) & (phosphosite_protSite <= end))
  if (nrow(ann) > 0){
    return(1)
  } else {
    return(0)
  }
}

getPhosphositeSite <- function(prot,start,end, ann_table){
  ann <- subset(ann_table, (protein_id==prot) & (phosphosite_protSite >= start) & (phosphosite_protSite <= end))
  if (nrow(ann) > 0){
    return(paste0(ann$phosphosite_protSite, collapse = ';'))
  } else {
    return(0)
  }
}

proteoform_annotation[,phophoPeptide := getPhosphositeAnn(protein_id,
                                                          PeptidePositionStart,
                                                          PeptidePositionEnd, 
                                                          phosphosite_annotation_table), 
                      by=c("id","protein_id")]

proteoform_annotation[,phophoSite := getPhosphositeSite(protein_id,
                                                        PeptidePositionStart,
                                                        PeptidePositionEnd, 
                                                        phosphosite_annotation_table), 
                      by=c("id","protein_id")]



### Fisher's exact test on regulated phosphosites
proteoform_annotation_proteinInfo <- unique(subset(proteoform_annotation, select=c("protein_id","proteoform_score","proteoform_score_pval_adj","phophoProtein","n_proteoforms")))

uniprot_annotation <- fread("/Users/isabell/Desktop/projects/ProteoformProject/InputData/uniprotSequenceAnnotation_humanSwissprot_130820.tab", quote = FALSE)
uniprot_annotation <- subset(uniprot_annotation, select=c('Entry','Alternative products (isoforms)'))
setnames(uniprot_annotation, 'Alternative products (isoforms)', 'isoform_annotation')
uniprot_annotation[, has_isoforms := ifelse(length(grep("ALTERNATIVE PRODUCTS",isoform_annotation)) > 0, 1, 0), by="Entry"]
uniprot_annotation <- subset(uniprot_annotation, select=c('Entry','has_isoforms'))

proteoform_annotation_proteinInfo <- merge(proteoform_annotation_proteinInfo, uniprot_annotation, by.x="protein_id", by.y="Entry", all.x=T, all.y=F, sort = F)

fisher_res <- data.table()
for (i in seq(0,1,0.05)){
  for (j in seq(0,1,0.05)){
    proteoform_sigP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==1) & (proteoform_score >= i) & (proteoform_score_pval_adj <= i)])
    proteoform_noP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==0) & (proteoform_score >= i) & (proteoform_score_pval_adj <= i)])
    noProteoform_sigP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==1) & ((proteoform_score < i) | (proteoform_score_pval_adj <= j))])
    noProteoform_noP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==0) & ((proteoform_score < i) | (proteoform_score_pval_adj <= j))])
    #noProteoform_sigP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==1) & ((proteoform_score < i) )])
    #noProteoform_noP <- nrow(proteoform_annotation_proteinInfo[(phophoProtein==0) & ((proteoform_score < i) )])
    
    contingencyTable <- matrix(c(proteoform_sigP, noProteoform_sigP, proteoform_noP, noProteoform_noP),
                               nrow = 2)
    
    fisherProteoformPhosphoEnrichment = fisher.test(contingencyTable, alternative = "greater")
    fisherProteoformPhosphoEnrichment = data.table(score_threshold=i,
                                                   qval_threshold=j,
                                                   test='fisherProteoformPhosphoEnrichment', 
                                                   pvalue=fisherProteoformPhosphoEnrichment$p.value,
                                                   odds=round(fisherProteoformPhosphoEnrichment$estimate, digits = 3))
    
    fisher_res <- rbind(fisher_res, fisherProteoformPhosphoEnrichment)
  }
}

pdf("fisherProteoformPhosphoEnrichment.pdf", width=4, height=3)
ggplot(fisher_res, aes(x=score_threshold, y=odds, color=-log10(pvalue), group=qval_threshold)) + 
  geom_point() + 
  geom_line() + 
  theme_classic() + 
  ylim(0,8)+
  #ggtitle('fisherProteoformPhosphoEnrichment') +
  scale_color_gradient2(low='grey', midpoint=-log10(0.01), mid='blue',high='red', space = "rgb", guide = "colourbar") 
dev.off()

write.table(fisher_res, "fisherProteoformPhosphoEnrichment.txt", quote=F, sep="\t", row.names = F, col.names = T)

###############

fisher_res <- data.table()
for (i in seq(0,1,0.05)){
  for (j in seq(0,1,0.05)){
    proteoform_sigP <- nrow(proteoform_annotation_proteinInfo[(has_isoforms==1) & (proteoform_score >= i) & (proteoform_score_pval_adj <= j)])
    proteoform_noP <- nrow(proteoform_annotation_proteinInfo[(has_isoforms==0) & (proteoform_score >= i) & (proteoform_score_pval_adj <= j)])
    noProteoform_sigP <- nrow(proteoform_annotation_proteinInfo[(has_isoforms==1) & ((proteoform_score < i) | (proteoform_score_pval_adj <= j))])
    noProteoform_noP <- nrow(proteoform_annotation_proteinInfo[(has_isoforms==0) & ((proteoform_score < i) | (proteoform_score_pval_adj <= j))])
    
    contingencyTable <- matrix(c(proteoform_sigP, noProteoform_sigP, proteoform_noP, noProteoform_noP),
                               nrow = 2)
    
    fisherProteoformIsoformEnrichment = fisher.test(contingencyTable, alternative = "greater")
    fisherProteoformIsoformEnrichment = data.table(score_threshold=i,
                                                   qval_threshold=j,
                                                   test='fisherProteoformIsoformEnrichment', 
                                                   pvalue=fisherProteoformIsoformEnrichment$p.value,
                                                   odds=round(fisherProteoformIsoformEnrichment$estimate, digits = 3))
    
    fisher_res <- rbind(fisher_res, fisherProteoformIsoformEnrichment)
  }
}

pdf("fisherProteoformIsoformEnrichment.pdf", width=4, height=3)
ggplot(fisher_res, aes(x=score_threshold, y=odds, color=-log10(pvalue),group=qval_threshold)) + 
  geom_point() + 
  geom_line() + 
  theme_classic() + 
  ylim(0,8)+
  scale_color_gradient2(low='grey', midpoint=-log10(0.01), mid='blue',high='red', space = "rgb", guide = "colourbar") 
dev.off()

write.table(fisher_res, "fisherProteoformIsoformEnrichment.txt", quote=F, sep="\t", row.names = F, col.names = T)


#########

### Check clusters for regulated phosphosites.
traces_location$trace_annotation[,phophoProtein := ifelse(protein_id %in% significantly_regulated_proteins, 1, 0)]
traces_location$trace_annotation[,phophoPeptide := getPhosphositeAnn(protein_id,PeptidePositionStart,PeptidePositionEnd,phosphosite_annotation_table), by=c("id","protein_id")]

traces_location$trace_annotation[,phospho_fisher_pval := 1.00]
traces_location$trace_annotation[,phospho_fisher_odds := 1.00]
#protOfInterest <- traces_location$trace_annotation[(proteoform_pval_adj <= 0.05) & (n_proteoforms > 1) & (phophoPeptide ==1)]$protein_id
protOfInterest <- traces_location$trace_annotation[phophoPeptide ==1]$protein_id

for (p in protOfInterest){
  ann <- subset(traces_location$trace_annotation, protein_id==p)
  proteoforms <- unique(ann$proteoform_id)
  for (q in proteoforms){
    ann_q <- subset(ann, proteoform_id==q)
    proteoform_sigP <- nrow(ann_q[phophoPeptide==1])
    proteoform_noP <- nrow(ann_q[phophoPeptide==0])
    ann_nq <- subset(ann, proteoform_id!=q)
    noProteoform_sigP <- nrow(ann_nq[phophoPeptide==1])
    noProteoform_noP <- nrow(ann_nq[phophoPeptide==0])
    contingencyTable <- matrix(c(proteoform_sigP, noProteoform_sigP, proteoform_noP, noProteoform_noP),
                               nrow = 2)
    fisher_out = fisher.test(contingencyTable, alternative = "greater")
    traces_location$trace_annotation[(protein_id==p) & (proteoform_id==q)]$phospho_fisher_pval = fisher_out$p.value
    traces_location$trace_annotation[(protein_id==p) & (proteoform_id==q)]$phospho_fisher_odds = fisher_out$estimate
  }
}

traces_location$trace_annotation[,phospho_fisher_pval_min := min(phospho_fisher_pval), by="protein_id"]

#' ## Save traces 
saveRDS(traces_location,"traces_location.rds")

#' ## Extract global stats
#' ## Extract global stats
getProteoformStats <- function(traces, score_threshold, qval_threshold, localization_threshold=0.05, phospho_threshold=0.2){
  tab <- traces$trace_annotation
  
  proteins <- unique(tab$protein_id)
  n_proteins <- length(proteins)
  write.table(proteins,"ProteinTables/proteins_total.txt", sep="\t", quote = F, col.names = F, row.names = F)
  
  tab_sub <- subset(tab, (proteoform_score>=score_threshold) & (proteoform_score_pval_adj<=qval_threshold))
  
  proteins_proteoforms <- unique(tab_sub$protein_id)
  n_proteins_proteoforms <- length(proteins_proteoforms)
  write.table(proteins_proteoforms,paste0("ProteinTables/proteoforms_score_",score_threshold,"qval_",qval_threshold,".txt"), sep="\t", quote = F, col.names = F, row.names = F)
  
  proteins_notExtreme <- unique(tab_sub[(genomLocation_pval_lim_min > localization_threshold) & (genomLocation_pval_min > localization_threshold)]$protein_id)
  proteins_extreme <- unique(tab_sub[(genomLocation_pval_lim_min <= localization_threshold) | (genomLocation_pval_min <= localization_threshold)]$protein_id)
  proteins_moreExtreme <- unique(tab_sub[(genomLocation_pval_lim_min <= localization_threshold) & (genomLocation_pval_min > localization_threshold)]$protein_id)
  proteins_equallyExtreme <- unique(tab_sub[(genomLocation_pval_min <= localization_threshold)]$protein_id)
  proteins_proteoformPhosphoEnriched <- unique(tab_sub[(phophoProtein==1) & (phospho_fisher_pval<=phospho_threshold) & (phospho_fisher_odds >= 1.5)]$protein_id)
  
  proteins_extreme_and_phospho <- proteins_extreme[proteins_extreme %in% proteins_proteoformPhosphoEnriched]
  proteins_extreme_not_phospho <- proteins_extreme[! proteins_extreme %in% proteins_proteoformPhosphoEnriched]
  proteins_phospho_not_extreme <- proteins_proteoformPhosphoEnriched[! proteins_proteoformPhosphoEnriched %in% proteins_extreme]
  proteins_not_phospho_not_extreme <- proteins_proteoforms[! proteins_proteoforms %in% proteins_proteoformPhosphoEnriched]
  proteins_not_phospho_not_extreme <- proteins_not_phospho_not_extreme[! proteins_not_phospho_not_extreme %in% proteins_extreme]
  
  n_proteins_notExtreme <- length(proteins_notExtreme)
  n_proteins_extreme <- length(proteins_extreme)
  n_proteins_moreExtreme <- length(proteins_moreExtreme)
  n_proteins_equallyExtreme <- length(proteins_equallyExtreme)
  n_proteins_proteoformPhosphoEnriched <- length(proteins_proteoformPhosphoEnriched)
  
  n_proteins_extreme_and_phospho <- length(proteins_extreme_and_phospho)
  n_proteins_extreme_not_phospho <- length(proteins_extreme_not_phospho)
  n_proteins_phospho_not_extreme <- length(proteins_phospho_not_extreme)
  n_proteins_not_phospho_not_extreme <- length(proteins_not_phospho_not_extreme)
  
  dt <- data.table(score_threshold = score_threshold,
                   qval_threshold = qval_threshold,
                   localization_threshold=localization_threshold,
                   phospho_threshold=phospho_threshold,
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
                   p_proteins_equallyExtreme=n_proteins_equallyExtreme/n_proteins_proteoforms,
                   n_proteins_proteoformPhosphoEnriched=n_proteins_proteoformPhosphoEnriched,
                   p_proteins_proteoformPhosphoEnriched=n_proteins_proteoformPhosphoEnriched/n_proteins_proteoforms,
                   n_proteins_extreme_and_phospho=n_proteins_extreme_and_phospho,
                   n_proteins_extreme_not_phospho=n_proteins_extreme_not_phospho,
                   n_proteins_phospho_not_extreme=n_proteins_phospho_not_extreme,
                   n_proteins_not_phospho_not_extreme=n_proteins_not_phospho_not_extreme
  )
  return(dt)
}


res <- data.table()
for (k in seq(0,1,0.05)){
  for (i in seq(0,1,0.05)) {
    for (j in seq(0, 0.5, 0.05)){
      res_i <- getProteoformStats(traces_location, score_threshold=i, qval_threshold=k, localization_threshold=j, phospho_threshold=j)
      res <- rbind(res, res_i)
    }
  }
}

protein_stats = subset(res, select=c("score_threshold", "qval_threshold",
                                     "localization_threshold", "phospho_threshold",
                                     "n_proteins",
                                     "n_proteins_proteoforms","n_proteins_extreme",
                                     "n_proteins_equallyExtreme","n_proteins_moreExtreme", 
                                     "n_proteins_proteoformPhosphoEnriched",
                                     "n_proteins_extreme_and_phospho",
                                     "n_proteins_extreme_not_phospho",
                                     "n_proteins_phospho_not_extreme",
                                     "n_proteins_not_phospho_not_extreme"))

knitr::kable(protein_stats)
write.table(protein_stats, "protein_proteoform_stats.txt", 
            quote=F, sep="\t", row.names = F, col.names = T)

res[,'pval_threshold':=localization_threshold]
protein_stats_plot = subset(res, select=c("score_threshold", "qval_threshold",
                                          "pval_threshold",
                                          "n_proteins_proteoforms","p_proteins_proteoforms",
                                          "n_proteins_extreme","p_proteins_extreme",
                                          "n_proteins_proteoformPhosphoEnriched",
                                          "p_proteins_proteoformPhosphoEnriched"))
protein_stats_plot <- melt(protein_stats_plot, id.vars = c('score_threshold','qval_threshold','pval_threshold'))

pdf("protein_proteoform_stats_qval01.pdf", width=8, height=6)
ggplot(protein_stats_plot[qval_threshold==0.1], aes(x=score_threshold,y=value, group=pval_threshold, color=factor(pval_threshold))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ variable, scales='free', nrow=3) +
  theme_classic() +
  scale_color_discrete()
dev.off()

pdf("protein_proteoform_stats_score01.pdf", width=8, height=6)
ggplot(protein_stats_plot[score_threshold==0.1], aes(x=qval_threshold,y=value, group=pval_threshold, color=factor(pval_threshold))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ variable, scales='free', nrow=3) +
  theme_classic() +
  scale_color_discrete()
dev.off()


pdf("protein_proteoform_stats_percentProteoformIncrease_qval01.pdf", width=5, height=3)
ggplot(protein_stats_plot[variable=="p_proteins_proteoforms"], 
       aes(x=score_threshold,y=value, group=qval_threshold, colour=qval_threshold)) + 
  geom_point() + 
  geom_line() + 
  ylab("fraction of proteins with proteoforms") +
  theme_classic() 
ggplot(protein_stats_plot[qval_threshold==0.1][variable=="p_proteins_proteoforms"], aes(x=score_threshold,y=value, group=pval_threshold)) + 
  geom_point() + 
  geom_line() + 
  #facet_wrap(~ variable, scales='free', nrow=3) +
  ylab("fraction of proteins with proteoforms") +
  theme_classic() +
  scale_color_discrete()
dev.off()

pdf("protein_proteoform_stats_percentProteoformIncrease_score01.pdf", width=5, height=3)
ggplot(protein_stats_plot[variable=="p_proteins_proteoforms"], 
       aes(x=qval_threshold,y=value, group=score_threshold, colour=score_threshold)) + 
  geom_point() + 
  geom_line() + 
  ylab("fraction of proteins with proteoforms") +
  theme_classic() 
ggplot(protein_stats_plot[score_threshold==0.1][variable=="p_proteins_proteoforms"], aes(x=qval_threshold,y=value, group=pval_threshold)) + 
  geom_point() + 
  geom_line() + 
  #facet_wrap(~ variable, scales='free', nrow=3) +
  ylab("fraction of proteins with proteoforms") +
  theme_classic() +
  scale_color_discrete()
dev.off()


# Annotate proteoforms across traces list

final_peptide_set <- unique(traces_location$trace_annotation$id)
pepTracesList_filtered_sub <- subset(pepTracesList_filtered, trace_subset_ids = final_peptide_set)

pepTracesList_filtered_sub_ann <- lapply(pepTracesList_filtered_sub, function(x){
  x$trace_annotation <- subset(x$trace_annotation, select = c("protein_id", "id"))
  x$trace_annotation <- merge(x$trace_annotation,
                              traces_location$trace_annotation,by=c("protein_id", "id"),sort=F,all.x=T,all.y=F)
  x
})
class(pepTracesList_filtered_sub_ann) <- "tracesList"

tracesList_location <- copy(pepTracesList_filtered_sub_ann)
saveRDS(tracesList_location, "tracesList_location.rds")


#' ## Plot all proteoform profiles
sigProteins <- subset(traces_location$trace_annotation, (proteoform_score>=0.1) & (proteoform_score_pval_adj<=0.1))
sigProteins <- unique(sigProteins[order(proteoform_score, decreasing = T)]$protein_id)

source('../CCprofilerAnalysis/thesis/traces_plotting.R')

pdf("allProteoformClusters_score01_qval01.pdf",height=4, width=6)
for (p in sigProteins){
  plotPeptideCluster(
    traces_location,p, closeGaps=T, PDF=F)
  plotSub <- subset(
    pepTracesList_filtered_sub_ann, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()

## Plot all proteoform profiles with high sequence proximity
intCases <- subset(traces_location$trace_annotation, (proteoform_score>=0.1) & (proteoform_score_pval_adj<=0.1) & genomLocation_pval_lim_min<=0.1)
intCases <- unique(intCases[order(proteoform_score, decreasing = T)]$protein_id)

pdf("allProteoformClusters_score01_qval01_pval_lim_min_01.pdf",height=4, width=6)
for (p in intCases){
  plotPeptideCluster(
    traces_location, p, closeGaps=T, PDF=F)
  plotSub <- subset(
    pepTracesList_filtered_sub_ann, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3),
                              "; pval=", 
                              round(min(plotSub$trace_annotation$genomLocation_pval_lim_min), digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()

## Plot all proteoform profiles with phosphoEnrichment
intCases <- subset(traces_location$trace_annotation, (proteoform_score>=0.1) & (proteoform_score_pval_adj<=0.1) & phospho_fisher_pval<=0.1)
intCases <- unique(intCases[order(proteoform_score, decreasing = T)]$protein_id)

pdf("allProteoformClusters_score01_qval01_phospho_fisher_pval_01.pdf",height=4, width=6)
for (p in intCases){
  plotPeptideCluster(
    traces_location, p, closeGaps=T, PDF=F)
  plotSub <- subset(
    pepTracesList_filtered_sub_ann, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3),
                              "; pval=", 
                              round(min(plotSub$trace_annotation$phospho_fisher_pval), digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()

source("../CCprofilerAnalysis/thesis/proteoformClusterPhosphosite.R")
library('ggpubr')
pdf("allProteoformClusters_score01_qval01_phospho_fisher_pval_01_plot.pdf", width=5, height=3)
for (p in intCases){
  plotPeptideClusterSigPhosphosites(traces_location, p, closeGaps = TRUE)
}
dev.off()







setwd("/Users/isabell/Desktop/projects/ProteoformProject/Results/")
library(data.table)
library(ggplot2)

# proteoforms at score threshold 0.25
proteoforms_score01_qval01 <- fread("ProteinTables/proteoforms_score_0.1qval_0.1.txt", header = FALSE, col.names = 'ID')
proteoforms <- proteoforms_score01_qval01$ID

# proteoform stats
proteoform_stats <- fread("protein_proteoform_stats.txt")
#proteoform_stats_score01_qval01 <- subset(proteoform_stats, score_threshold==0.25)
proteoform_stats_score01_qval01 <- subset(proteoform_stats, (score_threshold==0.1) & (qval_threshold==0.1))
proteoform_stats_score01_qval01_pval10 <- subset(proteoform_stats_score01_qval01, localization_threshold==0.10)

#########################
# General proteoform pie
#########################

pie_dt <- data.table(id=c('proteoforms','qno proteofroms'), 
                     n=c(proteoform_stats_score01_qval01_pval10$n_proteins_proteoforms,
                         proteoform_stats_score01_qval01_pval10$n_proteins-proteoform_stats_score01_qval01_pval10$n_proteins_proteoforms))

pdf('pie_proteoforms_count.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=id)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=8.5, direction = 1) +
  theme_void() +
  scale_fill_manual(values=c("green4","#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=id)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4","#999999")) +
  theme_classic()
print(q)
dev.off()


#########################
# Uniprot isoform annotationpie
#########################

uniprot_annotation <- fread("/Users/isabell/Desktop/projects/ProteoformProject/InputData/uniprotSequenceAnnotation_humanSwissprot_130820.tab", quote = FALSE)

uniprot_annotation <- subset(uniprot_annotation, select=c('Entry','Alternative products (isoforms)'))
setnames(uniprot_annotation, 'Alternative products (isoforms)', 'isoform_annotation')

uniprot_annotation[, has_isoforms := ifelse(length(grep("ALTERNATIVE PRODUCTS",isoform_annotation)) > 0, 1, 0), by="Entry"]
uniprot_annotation_iso <- subset(uniprot_annotation, has_isoforms==1)

n_proteoforms_with_iso = length(which(proteoforms %in% uniprot_annotation_iso$Entry))
n_proteoforms_without_iso = length(which(! proteoforms %in% uniprot_annotation_iso$Entry))

####### Enrichment compared to full uniprot
contingencyTable <- matrix(c(n_proteoforms_with_iso, nrow(uniprot_annotation[has_isoforms==1]), n_proteoforms_without_iso, nrow(uniprot_annotation[has_isoforms==0])),
                           nrow = 2)

fisherEnrichment = fisher.test(contingencyTable, alternative = "greater")

fisherEnrichment = data.table(test='fisherProteoformPhosphoEnrichment',
                              pvalue=fisherEnrichment$p.value,
                              odds=round(fisherEnrichment$estimate, digits = 3))

print(fisherEnrichment)

####### Enrichment compared to full SEC dataset
proteoforms_score0 <- fread("ProteinTables/proteoforms_score_0qval_1.txt", header = FALSE, col.names = 'ID')
all_proteins <- proteoforms_score0$ID
n_proteins_with_iso = length(which(all_proteins %in% uniprot_annotation_iso$Entry))
n_proteins_without_iso = length(which(! all_proteins %in% uniprot_annotation_iso$Entry))

contingencyTable <- matrix(c(n_proteoforms_with_iso, n_proteins_with_iso, n_proteoforms_without_iso, n_proteins_without_iso),
                           nrow = 2)

fisherEnrichment = fisher.test(contingencyTable, alternative = "greater")

fisherEnrichment = data.table(test='fisherProteoformPhosphoEnrichment',
                              pvalue=fisherEnrichment$p.value,
                              odds=round(fisherEnrichment$estimate, digits = 3))

print(fisherEnrichment)

#######

pie_dt <- data.table(uniprot=c('annotated isoforms','no annotated isoforms'), n=c(n_proteoforms_with_iso,n_proteoforms_without_iso))

pdf('pie_proteoforms_annotatedUniprotIsoforms.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=uniprot)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4","#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=uniprot)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4","#999999")) +
  theme_classic()
print(q)
dev.off()


#####

proteoforms_score0 <- fread("ProteinTables/proteoforms_score_0qval_1.txt", header = FALSE, col.names = 'ID')

proteoforms0 <- proteoforms_score0$ID

n_proteins_with_iso = length(which(proteoforms0 %in% uniprot_annotation_iso$Entry))
n_proteins_without_iso = length(which(! proteoforms0 %in% uniprot_annotation_iso$Entry))

pie_dt <- data.table(uniprot=c('annotated isoforms','no annotated isoforms'), n=c(n_proteins_with_iso,n_proteins_without_iso))

#########################
# Sequence proximity pie
#########################

pie_dt <- data.table(stat=c('<=',
                               '<',
                               '>'), 
                     n=c(proteoform_stats_score01_qval01_pval10$n_proteins_equallyExtreme,
                         proteoform_stats_score01_qval01_pval10$n_proteins_moreExtreme,
                         proteoform_stats_score01_qval01_pval10$n_proteins_proteoforms - proteoform_stats_score01_qval01_pval10$n_proteins_extreme))

pdf('pie_proteoforms_sequenceProximity.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4",'blue',"#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4",'blue',"#999999")) +
  theme_classic()
print(q)
dev.off()

#########################
# Sequence proximity pie
#########################


pie_dt <- data.table(stat=c('phospho & extreme',
                            'extreme',
                            'phospho',
                            'none'), 
                     n=c(proteoform_stats_score01_qval01_pval10$n_proteins_extreme_and_phospho,
                         proteoform_stats_score01_qval01_pval10$n_proteins_extreme_not_phospho,
                         proteoform_stats_score01_qval01_pval10$n_proteins_phospho_not_extreme,
                         proteoform_stats_score01_qval01_pval10$n_proteins_not_phospho_not_extreme))

pie_dt$stat <- factor(pie_dt$stat, level=c('phospho',
                                           'phospho & extreme',
                                           'extreme',
                                           'none')) 

pdf('pie_proteoforms_phospho.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) +
  theme_classic()
print(q)
dev.off()

#########################
# Assembly specific pie
#########################

assembly_specific_proteoforms_score01_qval01 <- fread("ProteinTables/assembly_specific_proteoforms_score_0.1qval_0.1.txt", header = FALSE, col.names = 'ID')
assembly_specific_proteoforms_score01_qval01 <- assembly_specific_proteoforms_score01_qval01$ID

pie_dt <- data.table(stat=c('assembly specific',
                            'not assembly specific'), 
                     n=c(length(assembly_specific_proteoforms_score01_qval01),
                         proteoform_stats_score01_qval01_pval10$n_proteins_proteoforms-length(assembly_specific_proteoforms_score01_qval01)))

pdf('pie_proteoforms_assembly_specific.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4","#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4","#999999")) +
  theme_classic()
print(q)
dev.off()

#########################
# Cell Cycle regulated pie
#########################

cellcycle_regulated_proteoforms_score01_qval01 <- fread("ProteinTables/cellcycle_regulated_proteoforms__score01_qval01.txt", header = FALSE, col.names = 'ID')
cellcycle_regulated_proteoforms_score01_qval01 <- cellcycle_regulated_proteoforms_score01_qval01$ID

pie_dt <- data.table(stat=c('assembly specific',
                            'not assembly specific'), 
                     n=c(length(cellcycle_regulated_proteoforms_score01_qval01),
                         proteoform_stats_score01_qval01_pval10$n_proteins_proteoforms-length(cellcycle_regulated_proteoforms_score01_qval01)))

pdf('pie_proteoforms_cellcycle_regulated.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4","#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4","#999999")) +
  theme_classic()
print(q)
dev.off()

#########################
# Cell Cycle & assembly pie
#########################

cellcycle_and_assembly_proteoforms_score01_qval01 <- cellcycle_regulated_proteoforms_score01_qval01[cellcycle_regulated_proteoforms_score01_qval01 %in% assembly_specific_proteoforms_score01_qval01]
cellcycle_not_assembly_proteoforms_score01_qval01 <- cellcycle_regulated_proteoforms_score01_qval01[! cellcycle_regulated_proteoforms_score01_qval01 %in% assembly_specific_proteoforms_score01_qval01]
assembly_not_cellcycle_proteoforms_score01_qval01 <- assembly_specific_proteoforms_score01_qval01[! assembly_specific_proteoforms_score01_qval01 %in% cellcycle_regulated_proteoforms_score01_qval01]
proteins_not_assembly_not_cellcycle <- proteoforms[! proteoforms %in% assembly_specific_proteoforms_score01_qval01]
proteins_not_assembly_not_cellcycle <- proteins_not_assembly_not_cellcycle[! proteins_not_assembly_not_cellcycle %in% cellcycle_not_assembly_proteoforms_score01_qval01]

pie_dt <- data.table(stat=c('cellcycle_and_assembly',
                            'cellcycle_not_assembly',
                            'assembly_not_cellcycle',
                            'not_assembly_not_cellcycle'), 
                     n=c(length(cellcycle_and_assembly_proteoforms_score01_qval01),
                         length(cellcycle_not_assembly_proteoforms_score01_qval01),
                         length(assembly_not_cellcycle_proteoforms_score01_qval01),
                         length(proteins_not_assembly_not_cellcycle)))

pdf('pie_proteoforms_cellcycle_and_assembly.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) +
  theme_classic()
print(q)
dev.off()

#########################
# Phospho & cell cycle pie
#########################
traces_list_pepClusters <- readRDS("tracesList_location.rds")

phospho_specific_proteoforms_score01_qval01 <- lapply(traces_list_pepClusters, function(x){unique(x$trace_annotation[(proteoform_score >= 0.1) & (proteoform_score_pval_adj <= 0.1) & (phospho_fisher_pval <= 0.1)]$protein_id)})
phospho_specific_proteoforms_score01_qval01 <- unlist(phospho_specific_proteoforms_score01_qval01)
phospho_specific_proteoforms_score01_qval01 <- unique(phospho_specific_proteoforms_score01_qval01)

cellcycle_and_phospho_proteoforms_score01_qval01 <- cellcycle_regulated_proteoforms_score01_qval01[cellcycle_regulated_proteoforms_score01_qval01 %in% phospho_specific_proteoforms_score01_qval01]
cellcycle_not_phospho_proteoforms_score01_qval01 <- cellcycle_regulated_proteoforms_score01_qval01[! cellcycle_regulated_proteoforms_score01_qval01 %in% phospho_specific_proteoforms_score01_qval01]
phospho_not_cellcycle_proteoforms_score01_qval01 <- phospho_specific_proteoforms_score01_qval01[! phospho_specific_proteoforms_score01_qval01 %in% cellcycle_regulated_proteoforms_score01_qval01]
proteins_not_phospho_not_cellcycle <- proteoforms[! proteoforms %in% phospho_specific_proteoforms_score01_qval01]
proteins_not_phospho_not_cellcycle <- proteins_not_phospho_not_cellcycle[! proteins_not_phospho_not_cellcycle %in% cellcycle_not_phospho_proteoforms_score01_qval01]

pie_dt <- data.table(stat=c('cellcycle_and_phospho',
                            'cellcycle_not_phospho',
                            'phospho_not_cellcycle',
                            'not_phospho_not_cellcycle'), 
                     n=c(length(cellcycle_and_phospho_proteoforms_score01_qval01),
                         length(cellcycle_not_phospho_proteoforms_score01_qval01),
                         length(phospho_not_cellcycle_proteoforms_score01_qval01),
                         length(proteins_not_phospho_not_cellcycle)))

pie_dt$stat <- factor(pie_dt$stat, level=c('phospho_not_cellcycle',
                                           'cellcycle_and_phospho',
                                           'cellcycle_not_phospho',
                                           'not_phospho_not_cellcycle')) 

pdf('pie_proteoforms_cellcycle_and_phospho.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
  theme_void() +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=stat)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("green4",'darkgreen','green',"#999999")) +
  theme_classic()
print(q)
dev.off()


######
######
######

source('../CCprofilerAnalysis/thesis/PaperAnalysis/CellCycleHela/traces_plotting.R')

design_matrix <- readRDS("design_matrix.rda")

pdf("proteins_not_assembly_not_cellcycle.pdf",height=4, width=6)
for (p in proteins_not_assembly_not_cellcycle){
  plotSub <- subset(
    traces_list_pepClusters, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plot.tracesList(plotSub, legend = T, aggregateReplicates=TRUE, 
                  design_matrix=design_matrix, error_bar=FALSE,
                  name=paste0(p, "; score=", round(plotSub$Interphase1$trace_annotation$proteoform_score[1], digits = 3)),
                  colour_by = "proteoform_id")
}
dev.off()



setwd("/Users/isabell/Desktop/projects/ProteoformProject/Mouse/")
library(data.table)
library(ggplot2)

# proteoforms at adj_pval_cutoff = 10%
proteoforms_score <- fread("ProteinTables/proteoforms_adj_pval_cutoff_0.1.txt", header = FALSE, col.names = 'ID')
proteoforms <- proteoforms_score$ID

# proteoform stats
proteoform_stats <- fread("protein_proteoform_stats.txt")
proteoform_stats_score <- subset(proteoform_stats, adj_pval_cutoff==0.1)
proteoform_stats_score[,localization_threshold := seq(0, 0.5, 0.05)]
proteoform_stats_score_pval10 <- subset(proteoform_stats_score, localization_threshold==0.10)

#########################
# General proteoform pie
#########################

pie_dt <- data.table(id=c('proteoforms','no proteofroms'), 
                     n=c(proteoform_stats_score_pval10$n_proteins_proteoforms,
                         proteoform_stats_score_pval10$n_proteins-proteoform_stats_score_pval10$n_proteins_proteoforms))

pdf('pie_proteoforms_count.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=id)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=2, direction = 1) +
  theme_void() +
  scale_fill_manual(values=c("#999999","green4")) 
print(p)
q <- ggplot(pie_dt, aes(x="", y=n, fill=id)) +
  geom_bar(stat="identity", position='dodge',alpha=0.4, color="white") +
  geom_text(aes(y = n,label = n)) +
  scale_fill_manual(values=c("#999999","green4")) +
  theme_classic()
print(q)
dev.off()


#########################
# Uniprot isoform annotationpie
#########################

uniprot_annotation <- fread("/Users/isabell/Desktop/projects/ProteoformProject/InputData/uniprotSequenceAnnotation_mouseSwissprot_170820.tab", quote = FALSE)

uniprot_annotation <- subset(uniprot_annotation, select=c('Entry','Alternative products (isoforms)'))
setnames(uniprot_annotation, 'Alternative products (isoforms)', 'isoform_annotation')

uniprot_annotation[, has_isoforms := ifelse(length(grep("ALTERNATIVE PRODUCTS",isoform_annotation)) > 0, 1, 0), by="Entry"]
uniprot_annotation_iso <- subset(uniprot_annotation, has_isoforms==1)

n_proteoforms_with_iso = length(which(proteoforms %in% uniprot_annotation_iso$Entry))
n_proteoforms_without_iso = length(which(! proteoforms %in% uniprot_annotation_iso$Entry))

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


#########################
# Sequence proximity pie
#########################


pie_dt <- data.table(stat=c('<=',
                               '<',
                               '>'), 
                     n=c(proteoform_stats_score_pval10$n_proteins_equallyExtreme,
                         proteoform_stats_score_pval10$n_proteins_moreExtreme,
                         proteoform_stats_score_pval10$n_proteins_proteoforms - proteoform_stats_score_pval10$n_proteins_extreme))

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




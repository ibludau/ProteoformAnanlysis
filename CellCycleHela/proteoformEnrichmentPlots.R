setwd("/Users/isabell/Desktop/projects/ProteoformProject/Results/")
library(data.table)
library(ggplot2)

# enrichment analysis of proteoforms at score threshold 0.25
# https://david.ncifcrf.gov/ 13.08.2020

dt_enrichment <- fread("ProteinTables/proteoform_vs_total_score_25.txt")
dt_enrichment <- subset(dt_enrichment, select=c("Category","Term","Count","PValue","Fold Enrichment"))
setnames(dt_enrichment, "Fold Enrichment", "FoldEnrichment")

dt_enrichment <- subset(dt_enrichment, PValue<=0.01)
dt_enrichment <- subset(dt_enrichment, Count>=5)
dt_enrichment <- subset(dt_enrichment, FoldEnrichment>=1.2)

dt_enrichment[, Category := ifelse(Category=="UP_KEYWORDS", "UniProt keywords", "UniProt sequence feature")]

dt_enrichment <- dt_enrichment[order(rank(Count))]

dt_enrichment$Term <- factor(dt_enrichment$Term, levels = dt_enrichment$Term)

pdf("ProteoformEnrichmentPlot.pdf", width=8, height=3)
ggplot(dt_enrichment, aes(x=Term, y=Count, fill=FoldEnrichment)) + 
  geom_bar(stat='identity', alpha=0.8) +
  theme_classic() +
  facet_wrap(~ Category, ncol=2, scales = "free_y") +
  coord_flip() +
  scale_fill_gradient(low = "#C0D8C0", high = "green4")
dev.off()

########

# TMEM108B

dt_enrichment <- fread("../../proteoformViews/correlated_Q9NUM4_david_enrichment.txt")
dt_enrichment <- subset(dt_enrichment, select=c("Category","Term","Count","PValue","Fold Enrichment"))
setnames(dt_enrichment, "Fold Enrichment", "FoldEnrichment")

dt_enrichment <- subset(dt_enrichment, PValue<=0.01)
dt_enrichment <- subset(dt_enrichment, Count>=5)
dt_enrichment <- subset(dt_enrichment, FoldEnrichment>=1.2)

#dt_enrichment[, Category := "GO cellular component")]

dt_enrichment <- dt_enrichment[order(rank(Count))]

dt_enrichment$Term <- factor(dt_enrichment$Term, levels = dt_enrichment$Term)

pdf("Q9NUM4_MembraneEnrichmentPlot.pdf", width=8, height=3)
ggplot(dt_enrichment, aes(x=Term, y=Count, fill=FoldEnrichment)) + 
  geom_bar(stat='identity', alpha=0.8) +
  theme_classic() +
  facet_wrap(~ Category, ncol=2, scales = "free_y") +
  coord_flip() +
  scale_fill_gradient(low = "#C0D8C0", high = "green4")
dev.off()



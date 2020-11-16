#' ## Load CCprofiler and setup the work environment
library(devtools)
options(warn=-1)

install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
library('CCprofiler')
library('data.table')
library('ggplot2')

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Mouse/")
knitr::opts_knit$set(
  root.dir = '/Users/isabell/Desktop/projects/ProteoformProject/Mouse/')

options(knitr.table.format = function() {
  if (knitr::is_latex_output()) 
    "latex" else "pandoc"
})

traces_location <- readRDS("traces_location.rds")

#' ## Proteoform annotation
traces_proteoforms <- annotateTracesWithProteoforms(
  traces_location, score_cutoff = 0.4)

removeOutlierPeptides <- function(traces){
  outlier_ids <- traces$trace_annotation[(n_proteoforms==2) & (cluster==100)]$id
  all_ids <- traces$trace_annotation$id
  non_outlier_ids <- all_ids[! all_ids %in% outlier_ids]
  sub_traces <- subset(traces, trace_subset_ids = non_outlier_ids)
  return(sub_traces)
}

traces_proteoforms_noO <- removeOutlierPeptides(traces_proteoforms)

proteoform_traces <- proteinQuantification(traces_proteoforms_noO,
                                                quantLevel="proteoform_id",
                                                topN = 1000,
                                                keep_less = TRUE,
                                                rm_decoys = TRUE,
                                                use_sibPepCorr = FALSE,
                                                use_repPepCorr = FALSE,
                                                full_intersect_only = FALSE,
                                                verbose = FALSE)

traces_plain <- melt(proteoform_traces$traces)
traces_plain[,tissue:="brain"]
traces_plain[,tissue:=ifelse(variable %in% seq(9,16,1), "BAT", tissue)]
traces_plain[,tissue:=ifelse(variable %in% seq(17,24,1), "heart", tissue)]
traces_plain[,tissue:=ifelse(variable %in% seq(25,32,1), "liver", tissue)]
traces_plain[,tissue:=ifelse(variable %in% seq(33,40,1), "quadriceps", tissue)]

traces_plain[,protein:=gsub("_.+","",id)]

traces_plain$tissue <- factor(traces_plain$tissue, levels=c("brain","BAT","heart","liver","quadriceps"))
traces_plain$id <- factor(traces_plain$id)

traces_plain[,log_int:=log2(value)]
traces_plain[,log_int:=ifelse(log_int=="-Inf",NA, log_int)]

traces_plain_proteofroms <- traces_plain[grep("_.+", id)]
all_proteins <- unique(traces_plain_proteofroms$protein)
anova_out <- data.table(id=all_proteins, anova_pval_tissue=1, anova_pval_proteoform=1, anova_pval_mixed=1)
for (p in all_proteins){
  aov.model <- aov(traces_plain_proteofroms[protein==p]$log_int ~ traces_plain_proteofroms[protein==p]$tissue*traces_plain_proteofroms[protein==p]$id)
  aov.summary <- summary(aov.model)[[1]]
  anova_out[id==p]$anova_pval_tissue <- aov.summary[["Pr(>F)"]][1]
  anova_out[id==p]$anova_pval_proteoform <- aov.summary[["Pr(>F)"]][2]
  anova_out[id==p]$anova_pval_mixed <- aov.summary[["Pr(>F)"]][3]
}

anova_out[,anova_pval_tissue_adjBF := p.adjust(anova_pval_tissue, "bonferroni")]
anova_out[,anova_pval_proteoform_adjBF := p.adjust(anova_pval_proteoform, "bonferroni")]
anova_out[,anova_pval_mixed_adjBF := p.adjust(anova_pval_mixed, "bonferroni")]

anova_out <- anova_out[order(anova_pval_mixed_adjBF)]
pdf("anova_out_sorted.pdf", width=8, height=5)
for (p in anova_out$id){
  plotSub <- subset(
    traces_location, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plotSub <- plotSub[c("traces","trace_type","trace_annotation","fraction_annotation")]
  class(plotSub) <- 'traces'
  plot(plotSub, legend = T, colour_by = "proteoform_id", name = paste0(p," pval=",anova_out[id==p]$anova_pval_mixed_adjBF))
}
dev.off()


pdf("anova_out_sorted_proteoformQuant.pdf", width=8, height=5)
for (p in anova_out$id){
  proteoform_traces$trace_annotation[,protein_id := gsub("_.+","",protein_id)]
  plotSub <- subset(
    proteoform_traces, trace_subset_ids = p, 
    trace_subset_type = "protein_id")
  plotSub <- plotSub[c("traces","trace_type","trace_annotation","fraction_annotation")]
  class(plotSub) <- 'traces'
  plot(plotSub, legend = T, colour_by = "proteoform_id", name = paste0(p," pval=",anova_out[id==p]$anova_pval_mixed_adjBF))
}
dev.off()

####

length(unique(anova_out[anova_pval_mixed_adjBF <= 0.01]$id))
length(unique(anova_out[anova_pval_mixed_adjBF > 0.01]$id))

#########################
# Tissue specific pie
#########################

pie_dt <- data.table(id=c('tissue specific','not tissue specific'), 
                     n=c(length(unique(anova_out[anova_pval_mixed_adjBF <= 0.01]$id)),
                         length(unique(anova_out[anova_pval_mixed_adjBF > 0.01]$id))))

pie_dt$id <- factor(pie_dt$id, levels=c('tissue specific','not tissue specific'))

pdf('pie_tissue_specific_proteoforms_count.pdf', height=2, width=5)
p <- ggplot(pie_dt, aes(x="", y=n, fill=id)) +
  geom_bar(stat="identity", width=1, alpha=0.4, color="white") +
  coord_polar("y", start=0, direction = -1) +
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



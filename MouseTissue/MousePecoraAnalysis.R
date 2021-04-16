library(devtools)
library(data.table)
install_github("jessegmeyerlab/PeCorA")
library(PeCorA)

setwd("/Users/isabell/Desktop/projects/ProteoformProject/Revisions/")

benchmark_traces <- readRDS("../Mouse/traces_location.rds")

benchmark_q <- benchmark_traces$traces
benchmark_q <- melt(benchmark_q, id.vars = "id", variable.name = "sample", value.name = "quant")

benchmark_qa <- merge(benchmark_q, benchmark_traces$trace_annotation, by = "id")

benchmark_qa[,"Peptide":=id]
benchmark_qa[,"Protein":=protein_id]
benchmark_qa[,"Peptide.Modified.Sequence":=id]

benchmark_qa[,"Condition":=ifelse(sample %in% c(1:8), "brain", "no")]
benchmark_qa[,"Condition":=ifelse(sample %in% c(9:16), "BAT", Condition)]
benchmark_qa[,"Condition":=ifelse(sample %in% c(17:24), "heart", Condition)]
benchmark_qa[,"Condition":=ifelse(sample %in% c(25:32), "liver", Condition)]
benchmark_qa[,"Condition":=ifelse(sample %in% c(33:40), "quad", Condition)]

rep_df = data.table(sample=as.factor(c(1:40)), BioReplicate = rep(c(1:8),5))
benchmark_qa = merge(benchmark_qa,rep_df,by='sample')

benchmark_qa[,"Normalized.Area":=quant]

sel_cols = names(benchmark_qa)[names(benchmark_qa) %in% names(t)]
benchmark_final <- subset(benchmark_qa, select=sel_cols)

benchmark_df <- setDF(benchmark_final)

scaled_peptides <- PeCorA_preprocessing(benchmark_df,
                                        area_column_name=6,
                                        threshold_to_filter=2,
                                        control_name="brain")

disagree_peptides <- PeCorA(scaled_peptides)

pdf("pecora_mouse_examples.pdf")
for (i in c(1:10)){
  PeCorA_plotting_plot<-PeCorA_plotting(disagree_peptides,disagree_peptides[i,],scaled_peptides)
  print(PeCorA_plotting_plot)
}
dev.off()

res <- data.table(disagree_peptides)
res[,adj_adj_pval  := p.adjust(adj_pval, "fdr")]

write.table(res, "pecora_mouse_res.tsv", sep="\t", quote = F, row.names = F)

hist(res$adj_pval)
hist(res$adj_adj_pval)


pecora_10 <- res[adj_pval<=0.1]
pecora_adj10 <- res[adj_adj_pval<=0.1]
pecora_10e10 <- res[adj_pval<=10e-10]
pecora_10e20 <- res[adj_pval<=10e-20]

pecora_10[,n_pep := .N, by="protein"]
pecora_adj10[,n_pep := .N, by="protein"]
pecora_10e10[,n_pep := .N, by="protein"]
pecora_10e20[,n_pep := .N, by="protein"]

length(unique(pecora_10$protein))/length(unique(res$protein))
length(unique(pecora_adj10$protein))/length(unique(res$protein))
length(unique(pecora_10e10$protein))/length(unique(res$protein))
length(unique(pecora_10e20$protein))/length(unique(res$protein))

pecora_10_2pep <- pecora_10[n_pep>=2]
length(unique(pecora_10_2pep$protein))/length(unique(res$protein))

pecora_10_prots <- unique(pecora_10$protein)
pecora_adj10_prots <- unique(pecora_adj10$protein)
pecora_10e10_prots <- unique(pecora_10e10$protein)
pecora_10e20_prots <- unique(pecora_10e20$protein)
pecora_10_2pep_prots <- unique(pecora_10_2pep$protein)

res_copf <- benchmark_traces$trace_annotation
res_copf_sig <- res_copf[(proteoform_score_pval_adj <= 0.1) & (proteoform_score>=0.1)]
res_prots_copf <- unique(res_copf_sig$protein_id)

all_proteins <- unique(benchmark_traces$trace_annotation$protein_id)

library(Vennerable)
pdf("pecora_venn.pdf")
venn <- Venn(list("PeCorA"=pecora_10_prots,"COPF"=res_prots_copf))
plot(venn, doWeights = TRUE)
venn <- Venn(list("PeCorA"=pecora_10_prots,"PeCorA adj10"=pecora_adj10_prots))
plot(venn, doWeights = TRUE)
venn <- Venn(list("PeCorA"=pecora_10_prots,"PeCorA 2pep"=pecora_10_2pep_prots))
plot(venn, doWeights = TRUE)
venn <- Venn(list("all"=all_proteins,"PeCorA"=pecora_10_prots,"PeCorA 10e-10"=pecora_10e10_prots,"COPF"=res_prots_copf))
plot(venn, doWeights = TRUE)
venn <- Venn(list("all"=all_proteins,"PeCorA"=pecora_10_prots,"COPF"=res_prots_copf))
plot(venn, doWeights = TRUE)
dev.off()


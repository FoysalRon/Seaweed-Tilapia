########Differential abundance analysis#######
#Lefse
library(microbiomeMarker)
#Load data
OTU<-read.csv("ASVtable.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxa.csv", header=TRUE, row.names=1)
META<-read.csv("metadata.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
pseq<-t(physeq)
library(MicrobiotaProcess)
#Rarefy the data
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps
# contrast must be specified for two groups comparison
mm_edger <- run_edger(
  ps,
  group = "Group",
  pvalue_cutoff = 0.05,
  p_adjust = "fdr"
)
mm_edger
# bar plot&cladogram
da1 <- plot_ef_bar(mm_edger) + theme (axis.text.y = element_text(size = 12)) + scale_fill_manual(values=c("#D55E00", "#009E73")) + theme(axis.text.y = element_text(hjust = 1, size=12, face="italic"))
da1
#LEfSe
mm_lefse <- run_lefse(
  ps,
  wilcoxon_cutoff = 0.05,
  group = "Group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)
# bar plot&cladogram
da5 <- plot_ef_bar(mm_lefse) + theme (axis.text.y = element_text(size = 12)) + scale_fill_manual(values=c("#D55E00", "#009E73")) + theme(axis.text.y = element_text(hjust = 1, size=12, face="italic"))
da5
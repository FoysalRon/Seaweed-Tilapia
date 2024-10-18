#############Co-occurrence network###############
#Load data
OTU<-read.csv("asv_table.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxa_info.csv", header=TRUE, row.names=1)
META<-read.csv("sample_info.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
physeq<-t(physeq)
#Load package
library("MicrobiotaProcess");packageVersion("MicrobiotaProcess")
#Rarefy the data
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps

#Transform data for MicroEco (plot to abundant phyla)
library('file2meco')
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
meco_dataset <- phyloseq2meco(ps)
# The parameter cor_method in trans_network is used to select correlation calculation method.
# default pearson or spearman correlation invoke R base cor.test, a little slow
t1 <- trans_network$new(dataset = meco_dataset, cor_method = "spearman", filter_thres = 0.001)
# return t1$res_cor_p list, containing two tables: correlation coefficient table and p value table
# require WGCNA package
if(!require("WGCNA")) install.packages("WGCNA", repos = BiocManager::repositories())
t1 <- trans_network$new(dataset = meco_dataset, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)
# construct network; require igraph package
t1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# invoke igraph cluster_fast_greedy function for this undirected network 
t1$cal_module(method = "cluster_fast_greedy")
# require rgexf package to be installed
t1$save_network(filepath = "network.gexf")

# use_col is used to select a column of t1$res_node_table
tmp <- t1$trans_comm(use_col = "module", abundance = FALSE)
tmp
tmp$otu_table[tmp$otu_table > 0] <- 1
tmp$tidy_dataset()
tmp$cal_abund()
tmp2 <- trans_abund$new(tmp, taxrank = "Phylum", ntaxa = 10)
tmp2$data_abund$Sample %<>% factor(., levels = rownames(tmp$sample_table))
tmp2$plot_line(xtext_angle = 30, color_values = RColorBrewer::brewer.pal(12, "Paired")) + ylab("OTUs ratio (%)")

t1$cal_sum_links(taxa_level = "Phylum")
# interactive visualization; require chorddiag package; see https://github.com/mattflor/chorddiag
t1$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
# From v1.2.0, method = "circlize" is available for conveniently saving the static plot
if(!require("circlize")) install.packages("circlize")
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))

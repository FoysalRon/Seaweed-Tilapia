###########Seaweed-Tilapia:Analysis############
setwd("path_to_directory")
getwd()
##Load library
library("microbiomeSeq"); packageVersion("microbiomeSeq")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("ggtree"); packageVersion("ggtree")
library("microbiome"); packageVersion("microbiome")
#library("hrbrthemes"); packageVersion("hrbrthemes")
#library("gcookbook"); packageVersion("gcookbook")
library("tidyverse"); packageVersion("tidyverse")
library("dplyr"); packageVersion("dplyr")
#library("coin"); packageVersion("coin")
library("VennDiagram"); packageVersion("VennDiagram")
#library("UpSetR"); packageVersion("UpSetR")
#library("patchwork"); packageVersion("patchwork")
library("RColorBrewer"); packageVersion("RColorBrewer")
library("ggpubr"); packageVersion("ggpubr")
library('microeco'); packageVersion('microeco')
library('ggtext'); packageVersion('ggtext')
library('paletteer'); packageVersion('paletteer')

#Set colour
MyPalette <- c ("grey30", "#FB8072", "#FF7F00", "#C74B1D", "#77CCCC", "#999999", "#117777", "#8B2323", "#777711",  "#117744", "#774411", "#FB8072", "#d45087", "#665191")
MyPal <- c ("grey30", "#4F7CBA", "#7873C0", "#A26DC2", "#CE69BE", "#EB73B3", "#FC719E", "#F64971", "#E03426", "#F06719", "#F89217", "#F8B620", "#D5BB21", "#A2B627", "#57A337", "#33A65C", "#21B087", "#30BCAD", "#2CB5C0", "#999999", "#1BA3C6", "#117777")
MyPalC <- c ("grey30", "#6F63BB", "#8A60B0", "#BA43B4", "#C7519C", "#D63A3A", "#FF7F0E", "#FFAA0E", "#FFBF50", "#BCBD22", "#78A641", "#2CA030", "#1BA3C6", "#999999", "#12A2A8", "#1F83B4", "#117777")

#Load data (all data)
otu_table_16S<-read.csv("asv_file.csv", header=TRUE, row.names=1)
taxonomy_table_16S<-read.csv("taxa_file.csv", header=TRUE, row.names=1)
sample_info_16S<-read.csv("sample_info_file.csv", header=TRUE, row.names=1)
env_data_16S<-read.csv("env_data_file.csv", header=TRUE, row.names=1)
# generally add the metadata
dataset <- microtable$new(otu_table = otu_table_16S, sample_table = sample_info_16S)
dataset
# Create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S)
dataset
dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset$tidy_dataset()
print(dataset)
dataset$sample_sums() %>% range
# As an example, use 10000 sequences in each sample
dataset$rarefy_samples(sample.size = 12211)
dataset$sample_sums() %>% range
# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]
dataset$save_table(dirpath = "basic_files", sep = ",")
# use default parameters
dataset$cal_abund()
# return dataset$taxa_abund
class(dataset$taxa_abund)
# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]
dataset$save_abund(dirpath = "taxa_abund")
# tab-delimited, i.e. mpa format
dataset$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)
# remove those unclassified
dataset$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)
# return dataset$alpha_diversity
class(dataset$alpha_diversity)
# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")
test <- dataset$merge_taxa(taxa = "Genus")
test
test <- dataset$merge_samples(use_group = "Group")
test

# select top 10 abundant Phyla.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)
# select top 10 abundant Phyla.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
phy1 <- t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)
phy2 <- phy1 + scale_fill_manual(values=MyPalette) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(size=12)) + theme(axis.title.y = element_text(size = 12))
phy2

#Set colour
MyPie <- c ("#774411", "#117744", "#117777", "#8B2323", "#777711", "#FB8072", "#FF7F00", "#C74B1D", "#d45087", "#77CCCC", "grey30", "#999999", "#FB8072", "#665191")

# Pie-chart - make a ggplot2 object
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)
pie1 <- t1$plot_pie(facet_nrow = 1, add_label = TRUE)
pie1 <- pie1 + scale_fill_manual(values=MyPie) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) 
pie1

# Pie-chart - make a ggplot2 object
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15, groupmean = "Group")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)
pie2 <- t1$plot_pie(facet_nrow = 1, add_label = TRUE)
pie2 <- pie2 + scale_fill_manual(values=MyPie) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold"))
pie2

MyPal <- c ("grey30", "#4F7CBA", "#7873C0", "#A26DC2", "#CE69BE", "#EB73B3", "#FC719E", "#F64971", "#E03426", "#F06719", "#F89217", "#F8B620", "#D5BB21", "#A2B627", "#57A337", "#33A65C", "grey50", "#30BCAD", "#2CB5C0", "#1BA3C6", "#117777")

# The facet parameter can be used to obtain replicate-specific barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20)
g1 <- t1$plot_bar(others_color = "grey70", facet = "Group", legend_text_italic = FALSE, color_values = paletteer_d("ggthemes::Hue_Circle"))
g1 + theme_classic() + theme(axis.title.y = element_text(size = 14))
g1 <- g1 + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic")) + theme(axis.title.y = element_text(size = 14)) + theme(strip.text.x = element_text(size = 10, colour = "grey15"))
g1
g2 <- g1 + scale_fill_manual(values = MyPal)
g2

#Load data
OTU<-read.csv("asv.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxa.csv", header=TRUE, row.names=1)
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

## Plot alpha-diversity
(p = plot_richness(ps, x = "Group", measures=c("Observed", "Shannon", "Simpson", "Chao")))
p2 <- p + geom_boxplot(data = p$data, aes(x = Group, y = value, fill = Group), alpha = 0.8, lwd =0.8, outlier.shape = NA) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face="bold")) + scale_fill_manual(values=c("grey30", "brown4")) + theme(legend.position = "NULL") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + geom_pwc(aes(group = Group), method = "t_test", label = "{p.adj.format}{p.adj.signif}", hide.ns = TRUE, vjust = -0.1)
p2$layers <- p2$layers[-1]
p2

#Venn diagram (group)
library('MicEco'); packageVersion('MicEco')
et_vn1 <- ps_venn(
  ps,
  group="Group",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels = list(cex = 1.2),
  size=18,
  hjust=1.0,
  fill=c("grey30", "brown4"))
et_vn1

## Beta-diversity
# Random tree
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
physeq1 = merge_phyloseq(physeq, META, random_tree)

#Group - CNT-GE
ordu.unwt.uni <- ordinate(physeq1, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(physeq1, ordu.unwt.uni, color="Group", shape="Group") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac")
unwt.unifrac <- unwt.unifrac + theme_bw()
unw_ord1 <- unwt.unifrac + geom_point(size = 4, alpha = 0.9) + theme_bw(base_size = 12, base_line_size = 1.5) + theme(text = element_text(size = 12, face = "bold"))+ stat_ellipse(aes(group = Group), geom="polygon",level=0.8,alpha=0.2) + labs("Group") + scale_color_manual(values = c("grey30", "brown4"), name = "Group") + scale_shape_manual(values = c(19, 15), name = "Group") + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white"))
unw_ord1

#Group - CNT-GE
ordu.wt.uni <- ordinate(physeq1, "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, ordu.wt.uni, color="Group", shape="Group") 
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac")
wt.unifrac <- wt.unifrac + theme_bw()
wt_ord1 <- wt.unifrac + geom_point(size = 4, alpha = 0.9) + theme_bw(base_size = 20, base_line_size = 1.5) + theme(text = element_text(size = 12, face = "bold"))+ stat_ellipse(aes(group = Group), geom="polygon",level=0.8,alpha=0.2) + labs("Group") + scale_color_manual(values = c("grey30", "brown4"), name = "Group") + scale_shape_manual(values = c(19, 15), name = "Group") + theme(legend.position = "right") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white"))
wt_ord1

#PERMANOVA_ANODIS
library(vegan)
metadf <- data.frame(sample_data(physeq1))
unifrac.dist <- UniFrac(physeq1, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE)
permanova <- adonis2(unifrac.dist ~ Group, data = metadf)
permanova

#PERMANOVA_ANODIS
library(vegan)
metadf <- data.frame(sample_data(physeq1))
unifrac.dist <- UniFrac(physeq1, weighted = FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE)
permanova <- adonis2(unifrac.dist ~ Group, data = metadf)
permanova


####Bray-Custis distance####
library('reshape2')
#Load data
otumat<-read.csv("otutable.csv", header=TRUE, row.names=1)
taxmat<-read.csv("taxtable.csv", header=TRUE, row.names=1)
metadata<-read.csv("metatable.csv", header=TRUE, row.names=1)
setdiff(rownames(otumat),taxmat$Genus)
rownames(taxmat) <- taxmat$Genus
taxmat =  as.matrix(taxmat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
sampledata = sample_data(metadata)
physeq1 = phyloseq(OTU, TAX, sampledata)
physeq1
sample_data(physeq1)$Group <- factor((sample_data(physeq1)$Group), levels=c("CNT","GE"))
relab_genera = transform_sample_counts(physeq1, function(x) x / sum(x) * 100)
head(otu_table(relab_genera)[,1:6])
abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]
sub_dist <- list()
groups_all <- sample_data(relab_genera)$Group
groups_all
for (group in levels(groups_all)) {
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}
braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))
head(df.bray)
p4 <- ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() +
  geom_boxplot(alpha=0.6) +
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("grey30", "brown4")) + scale_color_manual(values = c("grey30", "brown4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + geom_pwc(aes(group = group), method = "t_test", label = "{p.adj.format}{p.adj.signif}", hide.ns = TRUE, vjust = 0.5)
p4
# Arrange plots
b_row <- plot_grid(unw_ord1, wt_ord1, p4, labels = c('c', 'd', 'e'), rel_widths = c(1.85, 1.95, 1.3), nrow = 1, align = "hv")
b_row
#Save plot
ggsave("Fig 6e.tiff", units = c("in"), width=10.5, height=4.5, dpi=300, compression="lzw")
ord = ordinate(relab_genera, method="PCoA", distance="bray")
plot_ordination(relab_genera, ord, color="Group") + geom_point(size=4) + stat_ellipse(aes(group=Group))
## Rare curve
# for reproducibly random number
set.seed(1024)
rareres <- get_rarecurve(obj=ps, chunks=400)
prare2 <- ggrarecurve(obj=rareres,
                      factorNames="Group",
                      shadow=FALSE,
                      indexNames=c("Observe","Chao1", "Shannon", "Simpson")
) +
  scale_color_manual(values=c("grey30", "brown4"))+
  theme_bw()+
  theme(axis.text=element_text(size=10), axis.text.x=element_text(size=8, angle=45, hjust = 1), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
prare2

## Another rarefac
library("ranacapa")
p <- ggrare(ps, step = 1000, color = "Group", se = FALSE)

p <- p + geom_line(size=1) +
  scale_color_manual(values=c("grey30", "brown4")) +
  theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face="bold")) + scale_fill_manual(values=c("grey30", "brown4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + ylab("Species Richness") + xlab("Number of sequence") + theme(legend.position="bottom")
p

# Arrange plots
t_row <- plot_grid(p, et_vn1, p2, labels = c('a', 'b', 'c'), rel_widths = c(1.5, 1.0, 2), nrow = 1, align = "v")
t_row
# Arrange plots
b_row <- plot_grid(unw_ord1, wt_ord1, p4, labels = c('c', 'd', 'e'), rel_widths = c(1.85, 2, 1), nrow = 1, align = "hv")
b_row
plot_grid(t_row, b_row, ncol = 1, rel_heights = c(2, 2.0))
#Save plot
ggsave("Fig 6.tiff", units = c("in"), width=10.5, height=8.5, dpi=300, compression="lzw")
#Save venn to rotate
et_vn1
#saving the plot (medium definition)
dev.copy(tiff,'plot_venn_gr.tiff', width=4, height=2.0, units="in", res=300)
dev.off()

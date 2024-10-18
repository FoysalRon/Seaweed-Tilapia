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

#Set colour
MyPie <- c ("#774411", "#117744", "#117777", "#8B2323", "#777711", "#FB8072", "#FF7F00", "#C74B1D", "#d45087", "#77CCCC", "grey30", "#999999", "#FB8072", "#665191")

# Pie-chart - phylum
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)
pie1 <- t1$plot_pie(facet_nrow = 1, add_label = TRUE)
pie1 <- pie1 + scale_fill_manual(values=MyPie) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) 
pie1

MyPal <- c ("grey30", "#4F7CBA", "#7873C0", "#A26DC2", "#CE69BE", "#EB73B3", "#FC719E", "#F64971", "#E03426", "#F06719", "#F89217", "#F8B620", "#D5BB21", "#A2B627", "#57A337", "#33A65C", "grey50", "#30BCAD", "#2CB5C0", "#1BA3C6", "#117777")

# Pie-chart - genus
# The facet parameter can be used to obtain replicate-specific barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20)
g1 <- t1$plot_bar(others_color = "grey70", facet = "Group", legend_text_italic = FALSE, color_values = paletteer_d("ggthemes::Hue_Circle"))
g1 + theme_classic() + theme(axis.title.y = element_text(size = 14))
g1 <- g1 + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic")) + theme(axis.title.y = element_text(size = 14)) + theme(strip.text.x = element_text(size = 10, colour = "grey15"))
g1
g2 <- g1 + scale_fill_manual(values = MyPal)
g2
#! /bin/bash

qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trim-left-f 10 --p-trim-left-r 10 --p-trunc-len-f 260 --p-trunc-len-r 220 --o-table table.qza --o-representative-sequences rep-seqs.qza --p-n-threads 8 --o-denoising-stats denoising-stats.qza
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats
qiime feature-table filter-features --i-table table.qza --p-min-frequency 10 --o-filtered-table singletons-filtered-table.qza
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table singletons-filtered-table.qza --o-filtered-data singletons-filtered-rep-seqs.qza
qiime feature-table tabulate-seqs --i-data singletons-filtered-rep-seqs.qza --o-visualization singletons-filtered-rep-seqs
qiime alignment mafft --i-sequences singletons-filtered-rep-seqs.qza --o-alignment aligned-rep-seqs
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree
qiime feature-classifier classify-consensus-blast --i-query singletons-filtered-rep-seqs.qza --i-reference-taxonomy /home/faisalron/Documents/Bioinformatics/SILVA138/silva_138_tax_341_806.qza --i-reference-reads /home/faisalron/Documents/Bioinformatics/SILVA138/silva_138_seq_341_806.qza --o-search-results BLAST --o-classification taxonomy-B --p-perc-identity 0.80 --p-maxaccepts 1 --verbose
qiime metadata tabulate --m-input-file taxonomy-B.qza --o-visualization taxonomy.qzv
qiime taxa collapse --i-table singletons-filtered-table.qza --i-taxonomy taxonomy-B.qza --p-level 7 --output-dir blast-taxtable
qiime tools export --input-path blast-taxtable/collapsed_table.qza --output-path blast-taxtable2
biom convert -i blast-taxtable2/feature-table.biom -o blast-taxtable2/feature-table-collapsed.txt --to-tsv
qiime tools export --input-path  singletons-filtered-table.qza --output-path bl-exported-feature-table
qiime tools export --input-path taxonomy-B.qza --output-path bl_exported-feature-table
cp bl_exported-feature-table/taxonomy.tsv bl_biom-taxonomy.tsv
sed -i 's/Feature ID/#OTUID/g' bl_biom-taxonomy.tsv
sed -i 's/Taxon/taxonomy/g' bl_biom-taxonomy.tsv
sed -i 's/Consensus/confidence/g' bl_biom-taxonomy.tsv
biom add-metadata -i bl-exported-feature-table/feature-table.biom -o bl-table-taxonomy.biom --observation-metadata-fp bl_biom-taxonomy.tsv --sc-separated taxonomy
biom convert -i bl-table-taxonomy.biom -o bl-table-with-taxonomy.txt --to-tsv --header-key taxonomy

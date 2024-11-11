# Plant Microbiome
# Mairui Gao (mrgao@umd.edu)

conda activate qiime2-2023.5

qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path config.txt  --output-path  mjsample.qza  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize --i-data mjsample.qza --o-visualization mjsample.qzv

qiime tools export --input-path mjsample.qzv --output-path mjsample_statistic

qiime dada2 denoise-paired --i-demultiplexed-seqs mjsample.qza p-trunc-len-f 260 --p-trunc-len-r 250 --p-trim-left-f 10 --p-trim-left-r 0 --p-n-threads 8 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats stats-dada2.qza

qiime metadata tabulate --m-input-file stats-dada2.qza   --o-visualization stats-dada2.qzv

qiime tools export --input-path stats-dada2.qzv --output-path stats2

qiime feature-table summarize   --i-table table.qza   --o-visualization table.qzv   --m-sample-metadata-file sample.tsv

qiime tools export --input-path table.qzv --output-path table_stat

qiime tools export --input-path rep-seqs.qza --output-path rep-seqs 

qiime tools export --input-path table.qza --output-path table
biom convert -i table/feature-table.biom -o asv_table.txt  --table-type "OTU table" --to-tsv

qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path silva138_16s.tax --output-path ref-taxonomy.qza

qiime tools import --type 'FeatureData[Sequence]' --input-path silva138_16s.fasta --output-path ref-seqs.qza

qiime feature-classifier classify-consensus-blast --i-query rep-seqs.qza --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --p-maxaccepts 1 --o-classification taxonomy_blast.qza

qiime metadata tabulate --m-input-file taxonomy_blast.qza --o-visualization taxonomy_blast.qzv

qiime taxa barplot  --i-table table.qza  --i-taxonomy taxonomy_blast.qza  --m-metadata-file sample.txt --o-visualization taxa-bar-plots.qzv

qiime tools export --input-path  taxa-bar-plots.qzv --output-path  taxa-bar-plots
 
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy_blast.qza --p-exclude mitochondria,chloroplast --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime taxa filter-table   --i-table table-no-mitochondria-no-chloroplast.qza   --i-taxonomy taxonomy_blast.qza   --p-include p__ --o-filtered-table table-no-m-no-c-with-phyla.qza
qiime taxa barplot  --i-table table-no-m-no-c-with-phyla.qza  --i-taxonomy taxonomy_blast.qza  --m-metadata-file sample.txt --o-visualization table-no-m-no-c-with-phyla.qzv

qiime tools export --input-path table-no-m-no-c-with-phyla.qzv --output-path taxa-bar-plots-phyla-no-mitochondria-no-chloroplast 

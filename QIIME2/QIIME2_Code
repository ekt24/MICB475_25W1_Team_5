## Determine ASVs with DADA2
### Run Denoising and Clustering (Single-End)

```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 240 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
```

### Run Denoising and Clustering (Paired-End)

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
```


## Visualize ASVs stats

```bash
qiime feature-table summarize \
  --i-table table-final2.qza \
  --o-visualization final_table2.qzv \
  --m-sample-metadata-file ~/opioid_metadata_revised.txt
```

```bash
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-final2.qza \
  --o-visualization rep-seqs-final2.qzv
```

## Transfer .qzv file from server to local computer

```bash
scp root@10.19.139.118:~/rep-seqs-final2.qzv .
```
# Generate a tree for phylogenetic diversity analyses

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-final2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

## Alpha-rarefaction

```bash
qiime diversity alpha-rarefaction \
  --i-table table-final2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 42735 \
  --m-metadata-file Final_Altered_metadata.txt \
  --o-visualization alpha-rarefaction-final.qzv
```

## Transfer stats file and rarefaction qzv file for visualization

```bash
scp final_stats2.qza root@10.19.139.155:~
scp root@10.19.139.118:~/alpha-rarefaction-final.qzv .
```

## Taxonomic analysis
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```

```bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

## Taxonomy barplots

```bash
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /datasets/project_1/moving_pictures/sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

## Filter table to exclude mitochondria and chloroplast sequences

```bash
  qiime taxa filter-table \
  --i-table table-final2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza
```

## Generate a tree for phylogenetic diversity analyses

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-final2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

# Alpha-rarefaction with maximum sequencing depth as 42735

```bash
qiime diversity alpha-rarefaction \
  --i-table table-final2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 42735 \
  --m-metadata-file Final_Altered_metadata.txt \
  --o-visualization alpha-rarefaction-final.qzv
```

# Transfer stats file and rarefaction qzv file for visualization

```bash
scp final_stats2.qza root@10.19.139.155:~
scp root@10.19.139.118:~/alpha-rarefaction-final.qzv .
```

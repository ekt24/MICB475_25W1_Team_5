## Determine ASVs with DADA2
### Run Denoising and Clustering (Single-End)

```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 240 \
  --o-representative-sequences rep-seqs-final2.qza \
  --o-table table-final2.qza \
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

## Generate a Taxonomy Table using trained classifier
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-v34-classifier.qza \
  --i-reads rep-seqs-final2.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table-final2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Final_Altered_metadata.txt \
  --o-visualization taxa-bar-plots.qzv
```

## Generate a Tree for phylogenetic diversity analyses

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-final2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

## Alpha-rarefaction with maximum sequencing depth as 42735

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


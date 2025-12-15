## Train Classifier

### Extract the V3–V4 region from SILVA reference sequences
```bash
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 240 \
  --o-reads ref-seqs-v34.qza
```

### Train the Naive Bayes classifier
```bash
  qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-v34.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-v34-classifier.qza
```
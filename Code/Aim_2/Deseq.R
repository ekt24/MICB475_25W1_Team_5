# Loading libraries required for differential abundance analysis 
library(DESeq2)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(ggplot2)

# Load phyloseq object
load('data/object_final.RData')

# Add +1 to feature counts to avoid zeros that aren't compatible with DESeq2 when estimating size factors
phyloseq_object_plus1 <- transform_sample_counts(object_final, function(x) x+1)

# Converting phyloseq object to DESeq2 object using opioid status as predictor
deseq_object <- phyloseq_to_deseq2(phyloseq_object_plus1, ~Opioid.Substance)

# Run DESeq2 differential abundance analysis 
deseq_output <- DESeq(deseq_object)
res <- results(deseq_output, tidy=TRUE)
 
# Convert DESeq2 results to data frame for volcano plot 
deseq_df <- as.data.frame(res)

# Separate taxon IDs as separate column
deseq_df$taxa <- rownames(deseq_df)

# Remove NA values
deseq_df <- na.omit(deseq_df)

# Compute -log10 adjusted p-values
deseq_df$neg_log10_padj <- -log10(deseq_df$padj)

# Volcano plot
volcano_plot <- ggplot(deseq_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.8) +
  scale_color_manual(values = c("grey70", "red")) +
  labs(
    title = "Volcano Plot of Differential Abundance (DESeq2)",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Significant (padj < 0.05)"
  ) +
  theme_bw(base_size = 14)

#ggsave("plots/volcano_plot.png", volcano_plot, width = 8, height = 10)

# Changing volcano plot labels from taxa ID to genus and species
# Extract taxonomy table from original phyloseq object and convert to data frame 
tax <- as.data.frame(as(tax_table(object_final), "matrix"))
write.csv(tax, "data/taxonomy_table.csv", row.names = TRUE)

# Merge differential abundance results with taxonomy information using taxon IDs as key
res_sig_tax <- deseq_df %>%
  left_join(
    tax %>% tibble::rownames_to_column("taxa"),
    by = "taxa"
  )

# Create an initial label combining genus and species
res_sig_tax$label <- paste(res_sig_tax$Genus, res_sig_tax$Species, sep = " ")

# If species missing, use genus only for the label
res_sig_tax$label[res_sig_tax$Species == "" | is.na(res_sig_tax$Species)] <- 
  res_sig_tax$Genus[res_sig_tax$Species == "" | is.na(res_sig_tax$Species)]

# If label is still NA, change back to taxon ID
res_sig_tax$label[is.na(res_sig_tax$label)] <- res_sig_tax$taxa[is.na(res_sig_tax$label)]

# Filter for significant taxa and sort using fold change
sig <- res_sig_tax %>% 
  filter(padj < 0.05) %>% 
  arrange(log2FoldChange)

# Keep top 20 taxa for bar plot
sig <- sig %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice(1:20)

# Create an intermediate label with genus and species prefixes
sig$clean_label <- paste0("g_", sig$Genus, " s_", sig$Species)

# If species is missing, change back to genus only or taxon ID
sig$clean_label[is.na(sig$Species) | sig$Species == ""] <- sig$Genus[is.na(sig$Species) | sig$Species == ""]
sig$clean_label[is.na(sig$clean_label)] <- sig$taxa[is.na(sig$clean_label)]

# Remove "g__" prefix from genus names
sig$Genus <- gsub("g__", "", sig$Genus)

# Remove "s__" prefix from species names and underscores
sig$Species <- gsub("s__", "", sig$Species)
sig$Species <- gsub("_", " ", sig$Species)   

# If a species name is missing or not informative, replace with genus
sig$Species[is.na(sig$Species) | sig$Species == "" | grepl("uncultured", sig$Species, ignore.case = TRUE)] <- ""

# Create final clean label
sig$label <- trimws(paste(sig$Genus, sig$Species))

# Remove duplicated genus names 
sig$Species <- gsub("^" %+% sig$Genus %+% "\\s*", "", sig$Species)

# Remove underscores
# Remove repeated genus prefixes
sig$Species <- mapply(
  function(genus, species) {
    species <- gsub("_", " ", species)             # convert underscores
    species <- gsub(paste0("^", genus, "\\s*"), "", species)  # remove genus prefix
    species
  },
  sig$Genus,
  sig$Species
)

# Reconstruct final labels after cleaning
sig$label <- trimws(paste(sig$Genus, sig$Species))

# Reorder taxa in descending order wrt to log2 fold change (x-axis)
sig <- sig %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(label = factor(label, levels = unique(label)))

# Keep strongest effect per label
sig <- sig %>%
  group_by(label) %>%
  slice_max(abs(log2FoldChange), n = 1) %>%  # keep strongest effect per genus
  ungroup()

#Bar plot for top 20 significant taxa 
sig_taxa <- ggplot(sig, aes(x = log2FoldChange, y = label, fill = log2FoldChange)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Top 20 Significant Differential Taxa",
       x = "Log2 Fold Change",
       y = "Taxa (Genus species)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))
labs(
  title = "Top 20 Significant Differential Taxa",
  x = "Taxa (Genus, species)",
  y = "Log2 Fold Change",
  fill = "log2FoldChange"
) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 7))

ggsave("plots/sig_taxa.png", sig_taxa, width = 10, height = 7)

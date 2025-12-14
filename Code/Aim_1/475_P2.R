library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

metaFP <- "Final_Altered_metadata.txt" 
meta <- read.delim(file=metaFP)

otuFP <- "feature-table.txt" 
otu <- read.delim(file=otuFP, skip=1, row.names = 1)

taxFP <- "taxonomy.tsv"
tax <- read.delim(file=taxFP)

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#format

otu_mat <- as.matrix(otu[,-1])
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'X.SampleID'
SAMP <- sample_data(samp_df)
class(SAMP)

tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature.ID`
TAX <- tax_table(tax_mat)
class(TAX)

# build phyloseq object
object <- phyloseq(OTU, SAMP, TAX, phylotree)

#filter and finalize object

object_filt <- subset_taxa(object,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
object_filt_nolow <- filter_taxa(object_filt, function(x) sum(x)>5, prune = TRUE)
object_filt_nolow_samps <- prune_samples(sample_sums(object_filt_nolow)>100, object_filt_nolow)
object_final <- subset_samples(object_filt_nolow_samps, !is.na(month) & Opioid.Substance %in% c("Yes", "No"))

#object_final is the phyloseq object to be used for Aims 1/2
save(object_final, file = "object_final.RData")

#sample size based on Min.
summary(sample_sums(object_final))

#rarefied to 2 lengths
rarecurve(t(as.data.frame(otu_table(object_final))), cex=0.1)
object_rare1 <- rarefy_even_depth(object_final, rngseed = 1, sample.size = 1000)

#rarecurve(t(as.data.frame(otu_table(object_final))), cex=0.1)
#object_rare2 <- rarefy_even_depth(object_final, rngseed = 1, sample.size = 267)

#Shannon on Rare 1 (1000) & Rare 2 (267) 

shannon_df <- data.frame(
  Shannon = estimate_richness(object_rare1, measures = "Shannon")$Shannon,
  Opioid.Substance = sample_data(object_rare1)$Opioid.Substance
)

plot1 <- ggplot(shannon_df, aes(x = Opioid.Substance, y = Shannon)) +
  labs(y = "Shannon Diversity Index", x = "Opioid Use") +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +   # optional points
  theme_bw()

print(plot1)
ggsave("Shannon_boxplot.png", plot = plot1, width = 6, height = 4, dpi = 300)

#shannon1 <- estimate_richness(object_rare1, measures = "Shannon")$Shannon
#plot1 <- plot_richness(object_rare1, measures = c("Shannon"), x = "Opioid.Substance")
#print(plot1)
#ggsave(filename = "Shannon_plot.png", plot = plot1, 
#       width = 6, height = 4, units = "in", dpi = 300)

#Wilcox Rank Sum Test on Rare 1 (1000) & Rare 2 (267) 

wilcox_result1 <- wilcox.test(Shannon ~ Opioid.Substance, data = shannon_df)
wilcox_result1

save(wilcox_result1, file = "wilcox_result.RData")

#v1 had W=1029, p-value=0.3742, not significant

#shannon_df2 <- data.frame(
#  Shannon = shannon2,
#  Opioid.Substance = sample_data(object_rare2)$Opioid.Substance
#)
#shannon_df_final2 <- shannon_df2 %>%
#  filter(Opioid.Substance %in% c("Yes", "No"))
#wilcox_result2 <- wilcox.test(Shannon ~ Opioid.Substance, data = shannon_df_final1)

#v2 had W=1156, p-value=0.1899, not significant

#Bray & PERMANOVA, on Rare 1 (1000) & Rare 2 (267)
bray_dist1 <- phyloseq::distance(object_rare1, method = "bray")
save(bray_dist1, file = "bray_dist1.RData")

pn_1 <- adonis2(bray_dist1 ~ Opioid.Substance, data = data.frame(sample_data(object_rare1)), permutations = 999)
pn_1
save(pn_1, file = "permanova_result1.RData")

pcoa_bc1 <- ordinate(object_rare1, method="PCoA", distance=bray_dist1)
gg_pcoa1 <- plot_ordination(object_rare1, pcoa_bc1, color = "Opioid.Substance") + stat_ellipse()
print(gg_pcoa1)
ggsave(filename = "PCoA_plot.png", plot = gg_pcoa1, 
      width = 6, height = 5, units = "in", dpi = 300)

#v1 p=0.044, significant 

#bray_dist2 <- phyloseq::distance(object_rare2, method = "bray")
#pn_2 <- adonis2(bray_dist2 ~ Opioid.Substance, data = data.frame(sample_data(object_rare2)), permutations = 999)
#pn_2
#pcoa_bc2 <- ordinate(object_rare2, method="PCoA", distance=bray_dist2)
#gg_pcoa2 <- plot_ordination(object_rare2, pcoa_bc2, color = "Opioid.Substance") + stat_ellipse()
#print(gg_pcoa2)

#v2 p=0.069, not significant


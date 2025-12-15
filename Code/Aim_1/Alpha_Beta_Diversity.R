library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

#loading all necessary files from qiime

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

#number of samples remaning after filtering
nsamples(object_final)

#sample size based on Min.
summary(sample_sums(object_final))

#rarefied to 2 lengths
rarecurve(t(as.data.frame(otu_table(object_final))), cex=0.1)
object_rare1 <- rarefy_even_depth(object_final, rngseed = 1, sample.size = 19434)

#number of sample remainng after rarefraction
nsamples(object_rare1)

#Shannon Diversity
shannon_df <- data.frame(
  Shannon = estimate_richness(object_rare1, measures = "Shannon")$Shannon,
  Opioid.Substance = sample_data(object_rare1)$Opioid.Substance
)

#Plotting Shannon on a Boxplot
plot1 <- ggplot(shannon_df, aes(x = Opioid.Substance, y = Shannon)) +
  labs(y = "Shannon Diversity Index", x = "Opioid Use") +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  theme_bw()

print(plot1)
ggsave("Shannon_boxplot.png", plot = plot1, width = 6, height = 4, dpi = 300)

#Wilcox Rank Sum Test 

wilcox_result1 <- wilcox.test(Shannon ~ Opioid.Substance, data = shannon_df)
wilcox_result1

save(wilcox_result1, file = "wilcox_result.RData")

#W=479, p-value=0.4527, not significant

#Bray & PERMANOVA
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

#p=0.54, not significant 




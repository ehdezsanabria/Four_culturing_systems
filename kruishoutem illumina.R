# source("https://bioconductor.org/biocLite.R")
# 
# biocLite("phyloseq")
# biocLite("Biobase")

devtools::install_github("karthik/wesanderson")

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("SpiecEasi")
library("dplyr")
library("Phenoflow")
library("grid")
library("vegan")
library("igraph")
library("DESeq2")
source("functions.R")

#### Convert the OTU table from csv to phyloseq object, to make the abundance plots ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('Metadata Kruishoutem.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- interaction(metadata$Timepoint, metadata$Substrate, metadata$Fertilizer)

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T0" = "T0",
                                        "T1" = "T1",
                                        "T2" = "T2",
                                        "T3" = "T3",
                                        "T4" = "T4",
                                        "T5" = "T5",
                                        "T6" = "T6",
                                        "T7" = "T7",
                                        "T8" = "T8"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T0" = "T0",
                                                                          "T1" = "T1",
                                                                          "T2" = "T2",
                                                                          "T3" = "T3",
                                                                          "T4" = "T4",
                                                                          "T5" = "T5",
                                                                          "T6" = "T6",
                                                                          "T7" = "T7",
                                                                          "T8" = "T8"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_Kruishoutem.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus) 
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver004:Oliver144)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxa1 <- df_phy_genus_scale_pruned_remerged %>% 
  dplyr::filter(Substrate == "Soil") %>% 
  ggplot(aes(Fertilizer, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  facet_grid(Substrate~Timepoint)+
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxa1)

plot.taxa2 <- df_phy_genus_scale_pruned_remerged %>% 
  dplyr::filter(Substrate == "GB") %>% 
  ggplot(aes(Fertilizer, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  facet_grid(Substrate~Timepoint)+
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxa2)


tiff(file = "Relative abundance Kruishoutem2.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.taxa1, plot.taxa2, ncol = 1)
dev.off()

png(file = "Relative abundance Kruishoutem2.png", width = 12, height = 12, units = "in", res = 300)
cowplot::plot_grid(plot.taxa1, plot.taxa2, ncol = 1)
dev.off()

####RANDOM TRIAL FOR EACH SUBSTRATE####

# Add interaction factor levels to pool replicates
metadata <- read.csv('GB_organic.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- metadata$Timepoint

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T0" = "T0",
                                        "T1" = "T1",
                                        "T2" = "T2",
                                        "T3" = "T3",
                                        "T4" = "T4",
                                        "T5" = "T5",
                                        "T6" = "T6",
                                        "T7" = "T7",
                                        "T8" = "T8"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T0" = "T0",
                                                                          "T1" = "T1",
                                                                          "T2" = "T2",
                                                                          "T3" = "T3",
                                                                          "T4" = "T4",
                                                                          "T5" = "T5",
                                                                          "T6" = "T6",
                                                                          "T7" = "T7",
                                                                          "T8" = "T8"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_GB_organic.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus) 
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver004:Oliver144)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxaGBor <- df_phy_genus_scale_pruned_remerged %>% 
  ggplot(aes(Timepoint, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  ylab("Relative Abundance (%)")+ xlab("")+
  labs(title = "Organic Fertilizer")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxaGBor)

tiff(file = "Relative abundance_GB_organic.tiff", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaGBor)
dev.off()

png(file = "Relative abundance_GB_organic.png", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaGBor)
dev.off()

#### GB FISH ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('GB_fish.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- metadata$Timepoint

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T0" = "T0",
                                        "T1" = "T1",
                                        "T2" = "T2",
                                        "T3" = "T3",
                                        "T4" = "T4",
                                        "T5" = "T5",
                                        "T6" = "T6",
                                        "T7" = "T7",
                                        "T8" = "T8"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T0" = "T0",
                                                                          "T1" = "T1",
                                                                          "T2" = "T2",
                                                                          "T3" = "T3",
                                                                          "T4" = "T4",
                                                                          "T5" = "T5",
                                                                          "T6" = "T6",
                                                                          "T7" = "T7",
                                                                          "T8" = "T8"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_GB_fish.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus) 
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver008:Oliver120)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxaGBF <- df_phy_genus_scale_pruned_remerged %>% 
  ggplot(aes(Timepoint, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  ylab("Relative Abundance (%)")+ xlab("")+
  labs(title = "Fish-derived Fertilizer")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxaGBF)

tiff(file = "Relative abundance_GB_fish.tiff", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaGBF)
dev.off()

png(file = "Relative abundance_GB_fish.png", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaGBF)
dev.off()


### Soil Animal ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('Soil_animal.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- metadata$Timepoint

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T0" = "T0",
                                        "T1" = "T1",
                                        "T2" = "T2",
                                        "T3" = "T3",
                                        "T4" = "T4",
                                        "T5" = "T5",
                                        "T6" = "T6",
                                        "T7" = "T7",
                                        "T8" = "T8"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T0" = "T0",
                                                                          "T1" = "T1",
                                                                          "T2" = "T2",
                                                                          "T3" = "T3",
                                                                          "T4" = "T4",
                                                                          "T5" = "T5",
                                                                          "T6" = "T6",
                                                                          "T7" = "T7",
                                                                          "T8" = "T8"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_Soil_animal.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus) 
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver032:Oliver124)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxaSoilA <- df_phy_genus_scale_pruned_remerged %>% 
  ggplot(aes(Timepoint, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  ylab("Relative Abundance (%)")+ xlab("")+
  labs(title = "Animal manure Fertilizer")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxaSoilA)

tiff(file = "Relative abundance_Soil_animal.tiff", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaSoilA)
dev.off()

png(file = "Relative abundance_Soil_animal.png", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaSoilA)
dev.off()



### Soil Plant ####

# Add interaction factor levels to pool replicates
metadata <- read.csv('Soil_plant.csv')

metadata <- data.frame(metadata)
metadata$merge_factor <- metadata$Timepoint

metadata$Timepoint <- plyr::revalue(metadata$Timepoint, replace = 
                                      c("T0" = "T0",
                                        "T1" = "T1",
                                        "T2" = "T2",
                                        "T3" = "T3",
                                        "T4" = "T4",
                                        "T5" = "T5",
                                        "T6" = "T6",
                                        "T7" = "T7",
                                        "T8" = "T8"))
metadata$Timepoint <- factor(as.character(metadata$Timepoint), levels = c("T0" = "T0",
                                                                          "T1" = "T1",
                                                                          "T2" = "T2",
                                                                          "T3" = "T3",
                                                                          "T4" = "T4",
                                                                          "T5" = "T5",
                                                                          "T6" = "T6",
                                                                          "T7" = "T7",
                                                                          "T8" = "T8"))

# Create phyloseq object from csv file

emma_df <- read.csv("./Counts_Soil_plant.csv")

# Taxonomy table
tax_df <-emma_df %>%  dplyr::select(OTU:Genus) 
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- emma_df %>% dplyr::select(Oliver016:Oliver128)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-OTU)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum" ,"Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

# Merge replicate samples for visualization purposes
phy_genus <- merge_samples(physeq, "merge_factor")

# Check library sizes
hist(sample_sums(phy_genus), breaks = 100)
print(sample_sums(phy_genus))

# rescale to minimum sampling depth
phy_genus_scale <- scale_reads(phy_genus)
print(sample_sums(phy_genus_scale))

# Pool at genus level
phy_genus_scale = tax_glom(phy_genus_scale, "Genus")

# Calculate relative abundances
phy_genus_scale <- transform_sample_counts(phy_genus_scale, function(x) 100*x/sum(x))

# Select top 12 Genera
TopNGenera <- names(sort(taxa_sums(phy_genus_scale), TRUE)[1:12])
tax_table(phy_genus_scale)[!taxa_names(phy_genus_scale) %in% TopNGenera, "Genus"] <- "Other"

# psmelt into dataframe
df_phy_genus_scale_pruned <- psmelt(phy_genus_scale)

# Add the following part if the version of phyloseq is giving issues and does not add metadata
# correctly

# The following instractions are to replace metadata again 
# since the merging of samples removes this information

df_phy_genus_scale_pruned <- Filter(function(x)!all(is.na(x)), df_phy_genus_scale_pruned)
df_phy_genus_scale_pruned <- df_phy_genus_scale_pruned[, -c(4:11)]
df_phy_genus_scale_pruned <- dplyr::left_join(df_phy_genus_scale_pruned, metadata[-c(1)],
                                              by = c("Sample" = "merge_factor")) %>% distinct()

# adjust order Genus label so that "other" group is last and the other
# Genera are ordered according to their total abundance across data set
sum_table_genera <- df_phy_genus_scale_pruned %>% group_by(Genus) %>% summarize(sum_abund = sum(Abundance))

names_sorted_genera <- as.character(sum_table_genera$Genus)[order(sum_table_genera$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_genera <- c(names_sorted_genera[names_sorted_genera != "Other"],
                         names_sorted_genera[names_sorted_genera == "Other"]) # Put "other" group last

df_phy_genus_scale_pruned$Genus <- factor(as.character(df_phy_genus_scale_pruned$Genus),
                                          levels = names_sorted_genera)

# Finally we merge the Genera present in "Others"
df_phy_genus_scale_pruned_remerged <- df_phy_genus_scale_pruned %>% 
  group_by(Genus, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_phy_genus_scale_pruned_remerged <- left_join(df_phy_genus_scale_pruned_remerged,
                                                metadata[-c(1)], 
                                                by = c("Sample" = "merge_factor")) %>% distinct()

# Make barplot #
library (RColorBrewer)
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
#getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
#getPalette <- colorRampPalette(wes_palette(n=5, "Zissou1"))

# Generate a plot with the relative abundances pooled by location
library(ggplot2)
library(wesanderson)


plot.taxaSoilP <- df_phy_genus_scale_pruned_remerged %>% 
  ggplot(aes(Timepoint, y = sum_abund, 
             fill = Genus, color = Genus), data = .)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("gray50", 13))+
  ylab("Relative Abundance (%)")+ xlab("")+
  labs(title = "Malt sprouts (Plant-derived) Fertilizer")+
  scale_fill_manual(values = rev(getPalette(13)))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

print(plot.taxaSoilP)

tiff(file = "Relative abundance_Soil_plant.tiff", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaSoilP)
dev.off()

png(file = "Relative abundance_Soil_plant.png", width = 12, height = 9, units = "in", res = 300)
print(plot.taxaSoilP)
dev.off()

tiff(file = "Relative abundance_GB.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.taxaGBF, plot.taxaGBor, ncol = 1)
dev.off()

png(file = "Relative abundance_GB.png", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.taxaGBF, plot.taxaGBor, ncol = 1)
dev.off()

tiff(file = "Relative abundance_soil.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.taxaSoilA, plot.taxaSoilP, ncol = 1)
dev.off()

png(file = "Relative abundance_soil.png", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.taxaSoilA, plot.taxaSoilP, ncol = 1)
dev.off()


################################################################################
### Alpha diversity analysis
################################################################################

# Species richness and eveness indices
library(vegan)
library(ggplot2)
counts <- t(otu_df)
shannon <- diversity(counts)
totalspecies<-specnumber(counts)
Pielou <- diversity(counts)/log(totalspecies)

# Calculate diversity from rescaled OTU table
diversity_results <- Diversity_16S(scale_reads(physeq), brea = FALSE, thresh = 1000, R=100)
diversity_results <- data.frame(Sample = rownames(diversity_results), diversity_results)

# Merge diversity results with metadata
diversity_results <- left_join(diversity_results, metadata, by = c("Sample")) %>% 
  distinct()

# Plot diversity results
## Overall difference between locations - D2 (Inverse Simpson)
plot.D2_1 <- diversity_results %>% 
  dplyr::filter(Substrate == "Soil") %>% 
  ggplot(aes(Timepoint, y = D2, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Inverse Simpson")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_1)

plot.D2_2 <- diversity_results %>% 
  dplyr::filter(Substrate == "GB") %>% 
  ggplot(aes(Timepoint, y = D2, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Inverse Simpson")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Dark2"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_2)

tiff(file = "Alpha diversity Kruishoutem.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

png(file = "Alpha diversity Kruishoutem.png", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

## Observed diversity

plot.D2_1 <- diversity_results %>% 
  dplyr::filter(Substrate == "Soil") %>% 
  ggplot(aes(Timepoint, y = D0, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Richness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_1)

plot.D2_2 <- diversity_results %>% 
  dplyr::filter(Substrate == "GB") %>% 
  ggplot(aes(Timepoint, y = D0, 
             fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Richness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Dark2"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")


print(plot.D2_2)

tiff(file = "Observed diversity Kruishoutem.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

png(file = "Observed diversity Kruishoutem.png", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.D2_1, plot.D2_2, ncol = 1)
dev.off()

# Export diversity values with metadata
indices <- cbind(Pielou, shannon, totalspecies)

# You can use the values in these files for statistical analysis
# Try to run a mixed model in SAS or GraphPad or any other software
write.csv(file = "evenness_Kruishoutem.csv", indices)
write.csv(file = "diversity_Kruishoutem.csv", diversity_results)

## Plot Evenness

# Change factor names

evenness_results <- read.csv('evenness_Kruishoutem.csv', stringsAsFactors = TRUE)

evenness_results$Timepoint <- plyr::revalue(evenness_results$Timepoint, replace = 
                                              c("T0" = "T0",
                                                "T1" = "T1",
                                                "T2" = "T2",
                                                "T3" = "T3",
                                                "T4" = "T4",
                                                "T5" = "T5",
                                                "T6" = "T6",
                                                "T7" = "T7",
                                                "T8" = "T8"))

plot.ev1 <- evenness_results %>% 
  dplyr::filter(Substrate == "Soil") %>% 
  ggplot(aes(Timepoint, y = Pielou,fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Evenness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Accent"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")

print(plot.ev1)

plot.ev2 <- evenness_results %>% 
  dplyr::filter(Substrate == "GB") %>% 
  ggplot(aes(Timepoint, y = Pielou,fill = Fertilizer), data = .)+
  geom_jitter(shape = 21, size = 4, position=position_dodge(0.8),
              aes(shape=Fertilizer))+geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  ylab("Evenness")+ xlab(".")+
  scale_fill_manual(values=brewer.pal(n=8,"Dark2"))+
  facet_grid(Substrate~Fertilizer)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14,angle = 30, hjust = 1),
        strip.text = element_text(size = 14),
        legend.position = "right")

print(plot.ev2)


tiff(file = "Evenness Kruishoutem.tiff", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.ev1, plot.ev2, ncol = 1)
dev.off()

png(file = "Evenness Kruishoutem.png", width = 12, height = 9, units = "in", res = 300)
cowplot::plot_grid(plot.ev1, plot.ev2, ncol = 1)
dev.off()



################################################################################
### Beta diversity analysis
################################################################################

# Add metadata to phyloseq object

sample_data(physeq) <- metadata

# Rescale OTU table to account for library size differences
df_phy_scaled <- scale_reads(physeq)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Now lets transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled))

# And calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Lets run an exploratory permanova to see suggested effect sizes.
# Similar variances across treatments: 
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled)))
disper.seq_tpt <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Timepoint)

anova(disper.seq_tpt)
print(disper.seq_tpt)

disper.seq_fer <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Fertilizer)

anova(disper.seq_fer)
print(disper.seq_fer)


disper.seq_sub <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Substrate)

anova(disper.seq_sub)
print(disper.seq_sub)

# Permutations are constrained within each Location 
perm_results <- vegan::adonis(dist.seq ~ Timepoint + Fertilizer + Substrate, 
                              data = data.frame(sample_data(df_phy_scaled)))

# Add this information on plots
my_grob = grid::grobTree(textGrob(bquote(paste(r[Timepoint]^2 == 
                                                 .(round(100 * perm_results$aov.tab[1, 5], 1)), 
                                               "%")), x = 0.55, y = 0.40, 
                                  hjust = 0, gp = gpar(col = "black", 
                                                       fontsize = 12, fontface = "italic")))
my_grob2 = grid::grobTree(textGrob(bquote(paste(r[Fertilizer]^2 == 
                                                  .(format(round(100 * perm_results$aov.tab[2, 5],1), 
                                                           nsmall = 1)), "%")), x = 0.55, y = 0.35, 
                                   hjust = 0, gp = gpar(col = "black", fontsize = 12, 
                                                        fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Substrate]^2 == 
                                            .(round(100 * perm_results$aov.tab[3, 5], 1)), 
                                          "%")), x = 0.55, y = 0.30, hjust = 0, gp = gpar(col = "black", 
                                                                                          fontsize = 12, fontface = "italic")))
# Now we can plot the beta diversity plot
library(wesanderson)
 
p_beta_bulk1 <- pcoa.df %>%
  dplyr::filter(Substrate == "GB") %>% 
  ggplot(aes(x=Axis.1, y=Axis.2,shape=Fertilizer), data = .)+
  geom_point(alpha=0.7, size=5,aes(fill=Timepoint))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()+
  scale_fill_brewer(palette = "Dark2")+
  #scale_fill_brewer(palette = "Accent")+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  #facet_grid(Timepoint~.)+
  theme(strip.text.x = element_text(size=16))

print(p_beta_bulk1)

tiff(file = "Betadispersion Kruishoutem_GB.tiff", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk1)
dev.off()

png(file = "Betadispersion Kruishoutem_GB.png", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk1)
dev.off()


p_beta_bulk2 <- pcoa.df %>%
  dplyr::filter(Substrate == "Soil") %>% 
  ggplot(aes(x=Axis.1, y=Axis.2,shape=Fertilizer), data = .)+
  geom_point(alpha=0.7, size=5,aes(fill=Timepoint))+
  scale_shape_manual(values = c(22,25))+
  theme_bw()+
  scale_fill_brewer(palette = "Dark2")+
  #scale_fill_brewer(palette = "Accent")+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  #facet_grid(Timepoint~.)+
  theme(strip.text.x = element_text(size=16))

print(p_beta_bulk2)

tiff(file = "Betadispersion Kruishoutem_Soil.tiff", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk2)
dev.off()

png(file = "Betadispersion Kruishoutem_Soil.png", width = 12, height = 9, units = "in", res = 300)
print(p_beta_bulk2)
dev.off()
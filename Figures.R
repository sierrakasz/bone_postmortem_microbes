#Kaszubinski 2021

#set up directory
rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Lucas_grant_project/R_scripts/")

#packages
library(car)
library(ggpubr)
library(microbiome)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)
library(vegan)

#set.seed
set.seed(5678)

#import three files: otu table, taxonomy table, and phylo tree
otu <- read.csv("table.csv")
tax <- read.csv("tax_format.csv")
tree <- read_tree('tree.nwk')

#load in metadata
metadata=(read.csv("aquatic_bones_lucas_grant_metadata.csv",header=TRUE))
#change metadata to factors for sorting 
metadata$pig <- as.factor(metadata$pig)
metadata$Year <- as.factor(metadata$Year)
metadata$date_collected <- factor(metadata$date_collected, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                                      '4/23/2019', '7/26/2019', '10/26/2019',
                                                                      '1/24/2020'))
#format metadata into phyloseq format
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

#format otu table in phyloseq format
rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
#merge tree and otu tables
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

#format taxonomy table in phyloseq format
rownames(tax) <- tax$OTUID
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
#only pulling taxa to genus
tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)
taxa_names(TAX)=row.names(OTU)

#merge it all together into one phyloseq object
physeq <- merge_phyloseq(physeq_otu.tree, TAX)
physeq
#7879 taxa, 71 samples

#triming out taxa that are not representative of .01% of mean sequence number
physeq_trim <- prune_taxa(taxa_sums(physeq) > sum(otu) *.001 / 71, physeq)
physeq_trim
#2479 taxa

#rarefy
physeq_7000 <- rarefy_even_depth(physeq_trim, sample.size = 7000)
physeq_7000
#2475, 65 samples

#removes positive and negative control samples
physeq_in <- subset_samples(physeq_7000, sample_type == 'internal_microbiome')
physeq_ex <- subset_samples(physeq_7000, sample_type == 'external_microbiome')
physeq_npn <- merge_phyloseq(physeq_in, physeq_ex)

# Figure 2 --------------------------------------------------------------

#separate out seasons for separate panel figures
physeq_fall <- subset_samples(physeq_npn, Season == 'Fall')
physeq_spri <- subset_samples(physeq_npn, Season == 'Spring')
physeq_summ <- subset_samples(physeq_npn, Season == 'Summer')
physeq_wint <- subset_samples(physeq_npn, Season == 'Winter')

#calculate beta div with unifrac dissimilarity matrix
ord_f = ordinate(physeq_fall, method="PCoA", distance="unifrac")
#create plot
ordplot_f=plot_ordination(physeq_fall, ord_f, color="Year", shape = 'Season') +
  scale_color_manual(values = c('#8A6C64', '#92250B')) +
  scale_shape_manual(values=c(15)) + geom_point(size = 5) +
  theme(legend.position="bottom")

ordplot_f

#calculate beta div with unifrac dissimilarity matrix
ord_s = ordinate(physeq_summ, method="PCoA", distance="unifrac")
#create plot
ordplot_s=plot_ordination(physeq_summ, ord_s, color="Year", shape = 'Season') +
  scale_color_manual(values = c('#E0C08E', '#E49417')) +
  scale_shape_manual(values=c(16)) + geom_point(size = 5) +
  theme(legend.position="bottom")

ordplot_s

#calculate beta div with unifrac dissimilarity matrix
ord_w = ordinate(physeq_wint, method="PCoA", distance="unifrac")
#create plot
ordplot_w=plot_ordination(physeq_wint, ord_w, color="Year", shape = 'Season') +
  scale_color_manual(values = c('#8E91E0', '#141BDC')) +
  scale_shape_manual(values=c(17)) + geom_point(size = 5) +
  theme(legend.position="bottom")

ordplot_w

#calculate beta div with unifrac dissimilarity matrix
ord_t = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot_t=plot_ordination(physeq_npn, ord_t, color="Year", shape = 'Season') +
  scale_color_manual(values = c('#BBD8C2','#1F8334')) +
  scale_shape_manual(values=c(15,18,16,17)) + geom_point(size = 5) 

ordplot_t

theme_set(theme_classic(base_size = 14))
tiff("FIG2.TIF", width = 3000, height = 3000, res=300)
ggarrange(ordplot_t, ordplot_f, ordplot_s, ordplot_w, 
          labels = c("A.", "B.", 'C.', 'D.'),
          nrow = 2, ncol = 2)
dev.off()

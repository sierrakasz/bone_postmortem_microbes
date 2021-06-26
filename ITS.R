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
library(reshape2)
library(tidyverse)
library(vegan)

#set.seed
set.seed(5678)


biom <- import_biom('SFKBonesApr2021.biom', 'SFKBonesTree.nwk', parseFunction=parse_taxonomy_greengenes)

#load in metadata
metadata=(read.csv("aquatic_bones_lucas_grant_metadata_euks.csv",header=TRUE))
#change metadata to factors for sorting 
metadata$pig <- as.factor(metadata$pig)
metadata$Year <- as.factor(metadata$Year)
metadata$date_collected <- factor(metadata$date_collected, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                                      '4/23/2019', '7/26/2019', '10/26/2019',
                                                                      '1/24/2020'))
#format metadata into phyloseq format
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

physeq <- merge_phyloseq(biom, sampdat)
physeq
#3950 taxa

#triming out taxa that are not representative of .01% of mean sequence number
otu <- data.frame(otu_table(physeq))
physeq_trim <- prune_taxa(taxa_sums(physeq) > sum(otu) *.001 / 65, physeq)
physeq_trim
#2469 taxa

#determine what to rarefy to
#show rarefaction curves
tiff("SFIG1_euk.TIF", width = 2000, height = 1500, res=300)
rarecurve(t(otu_table(physeq_trim)), step=1000, cex=0.3)
dev.off()
#4000

#rarefy
physeq_4000 <- rarefy_even_depth(physeq_trim, sample.size = 4000)
physeq_4000
#2467, 65 samples

#removes positive and negative control samples
physeq_in <- subset_samples(physeq_4000, sample_type == 'internal_microbiome')
physeq_ex <- subset_samples(physeq_4000, sample_type == 'external_microbiome')

physeq_npn <- merge_phyloseq(physeq_in, physeq_ex)


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

### bacteria
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

### eukaryotes

biom <- import_biom('SFKBonesApr2021.biom', parseFunction=parse_taxonomy_greengenes)
#fix the tax table formatting 
tax_biom <- as.data.frame(tax_table(biom))
tax_biom <- tax_biom[,-c(2:7)]

otu_biom <- as.data.frame(otu_table(biom))

tree_euk <- read_tree('SFKBonesTree.nwk')

#format otu table in phyloseq format
OTU_biom=otu_table(otu_biom, taxa_are_rows=TRUE)

#merge tree and otu tables
physeq_otu.tree=phyloseq(OTU_biom, tree_euk)

#format taxonomy table in phyloseq format
empty_as_na <- function(x){
  ifelse(is.na(x), "Unassigned", x)
}

tax_biom <- apply(tax_biom, 2, empty_as_na)

tax_reduc <- merge(tax_biom, otu_biom, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
#only pulling taxa to genus
tax_f <- tax_reduc[,1:6]
tax_f <- as.matrix(tax_f)
TAX_biom=tax_table(tax_f)
taxa_names(TAX_biom)=row.names(OTU_biom)

#merge it all together into one phyloseq object
physeq_euk <- merge_phyloseq(physeq_otu.tree, TAX_biom)
physeq_euk
#3950 taxa, 65 samples

#load in metadata
metadata_euk=(read.csv("aquatic_bones_lucas_grant_metadata_euks.csv",header=TRUE))
#change metadata to factors for sorting 
metadata_euk$pig <- as.factor(metadata_euk$pig)
metadata_euk$Year <- as.factor(metadata_euk$Year)
metadata_euk$date_collected <- factor(metadata_euk$date_collected, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                                              '4/23/2019', '7/26/2019', '10/26/2019',
                                                                              '1/24/2020'))
#format metadata into phyloseq format
sampdat=sample_data(metadata_euk)
sample_names(sampdat)=metadata_euk$SampleID

physeq_euk <- merge_phyloseq(physeq_euk, sampdat)
physeq_euk
#3950 taxa

#triming out taxa that are not representative of .01% of mean sequence number
otu <- data.frame(otu_table(physeq_euk))
physeq_trim_euk <- prune_taxa(taxa_sums(physeq_euk) > sum(otu) *.001 / 65, physeq_euk)
physeq_trim_euk

# Figure S1 ---------------------------------------------------------------

#determine what to rarefy to
#show rarefaction curves
tiff("SFIG1.TIF", width = 2000, height = 1500, res=300)
rarecurve(t(otu_table(physeq_trim)), step=1000, cex=0.3)
dev.off()

#rarefy
physeq_7000 <- rarefy_even_depth(physeq_trim, sample.size = 7000)
physeq_7000
#2475, 65 samples


# Figure S2 ---------------------------------------------------------------

#determine what to rarefy to
#show rarefaction curves
tiff("SFIG2.TIF", width = 2000, height = 1500, res=300)
rarecurve(t(otu_table(physeq_trim_euk)), step=1000, cex=0.3)
dev.off()
#4000

#rarefy
physeq_4000 <- rarefy_even_depth(physeq_trim_euk, sample.size = 4000)
physeq_4000
#2467, 65 samples

#removes positive and negative control samples
physeq_in_euk <- subset_samples(physeq_4000, sample_type == 'internal_microbiome')
physeq_ex_euk <- subset_samples(physeq_4000, sample_type == 'external_microbiome')

physeq_npn_euk <- merge_phyloseq(physeq_in_euk, physeq_ex_euk)

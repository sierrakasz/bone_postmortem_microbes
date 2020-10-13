rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Lucas_grant_project/R_scripts/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

#import three files: otu table, taxonomy table, and phylo tree
otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

# to note: average ST of confidence for tax assignment was: 0.932 +/- 0.097

#take the OTU names, make them row names to combine with the tax file
rownames(otu) <- otu$OTUID
otu <- otu[,-1]

#format it into a phyloseq obj.
OTU=otu_table(otu, taxa_are_rows=TRUE)
physeq_all=phyloseq(OTU,tree)
physeq_all
#taxa: 7879; 71 samples

#get the tax table into a useful format
#before importing, erased confidence section (av and std above)
#text to column in excel using semi colon as separater
#add the column names
tax <- tax %>% 
  mutate(Kingdom = str_replace(Kingdom, "D_0__", "")) %>% 
  mutate(Phylum = str_replace(Phylum, "D_1__", "")) %>% 
  mutate(Class = str_replace(Class, "D_2__", "")) %>% 
  mutate(Order = str_replace(Order, "D_3__", "")) %>% 
  mutate(Family = str_replace(Family, "D_4__", "")) %>% 
  mutate(Genus = str_replace(Genus, "D_5__", "")) %>% 
  mutate(Species = str_replace(Species, "D_6__", ""))
#get rids of an empty cells
tax[tax==""]<-"Unassigned"

#save file in case
write.csv(tax, "tax_format.csv")

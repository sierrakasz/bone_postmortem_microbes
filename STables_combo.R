#Kaszubinski 2021

#set up directory
rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Lucas_grant_project/R_scripts/")

#packages
library(car)
library(FEAST)
library(ggpubr)
library(microbiome)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(RCurl)
library(rsample)
library(reshape2)
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

#rarefy
physeq_7000 <- rarefy_even_depth(physeq_trim, sample.size = 7000)
physeq_7000
#2475, 65 samples

#removes positive and negative control samples
physeq_in <- subset_samples(physeq_7000, sample_type == 'internal_microbiome')
physeq_ex <- subset_samples(physeq_7000, sample_type == 'external_microbiome')
physeq_npn <- merge_phyloseq(physeq_in, physeq_ex)

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

rarefy
physeq_4000 <- rarefy_even_depth(physeq_trim_euk, sample.size = 4000)
physeq_4000
#2464, 65 samples

#removes positive and negative control samples
physeq_in_euk <- subset_samples(physeq_4000, sample_type == 'internal_microbiome')
physeq_ex_euk <- subset_samples(physeq_4000, sample_type == 'external_microbiome')

physeq_npn_euk <- merge_phyloseq(physeq_in_euk, physeq_ex_euk)


# Table S1 ----------------------------------------------------------------
#alpha-diversity among microbiome types 
# A & B) Mean and standard deviation of alpha diversity metrics
# C & D) Kruskal-Wallis test among alpha diversity metrics

# A
#organize data
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(sample_type) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
#organize data
erich_euk <- estimate_richness(physeq_npn_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich_euk <- add_rownames(erich_euk, "SampleID")
#summary stats
erich_sums <- merge(erich_euk, metadata_euk)
erich_sums %>% group_by(sample_type) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# C
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

#list to hold alpha div metrics 
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ sample_type, data = kw_values[[i]]))
}

# D
# get data organized and formatted for stats testing
erich_euk <- erich_euk %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich_euk, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

#list to hold alpha div metrics 
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ sample_type, data = kw_values[[i]]))
}

# Table S2 ----------------------------------------------------------------

##beta-diversity among microbiome types (internal, external)
# A & B) summary stats
# C & D) PERMANOVA results for comparing beta-diversity and beta-dispersion for 999 permuations. 

###bacteria 
# A

#calculate unifrac distances 
dis_unifrac_in <- UniFrac(physeq_in, weighted=FALSE, 
                          normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_in  <- melt(as.matrix(dis_unifrac_in))

dis_unifrac_ex <- UniFrac(physeq_ex, weighted=FALSE, 
                          normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_ex  <- melt(as.matrix(dis_unifrac_ex))

#find the average for each sample
df_unifrac_in <- df_unifrac_in[,-1]
df_unifrac_in_sum <- df_unifrac_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_ex <- df_unifrac_ex[,-1]
df_unifrac_ex_sum <- df_unifrac_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_in_sum$value)
sd(df_unifrac_in_sum$value)

mean(df_unifrac_ex_sum$value)
sd(df_unifrac_ex_sum$value)

# B
#calculate unifrac distances 
dis_unifrac_in <- UniFrac(physeq_in_euk, weighted=FALSE, 
                          normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_in  <- melt(as.matrix(dis_unifrac_in))

dis_unifrac_ex <- UniFrac(physeq_ex_euk, weighted=FALSE, 
                          normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_ex  <- melt(as.matrix(dis_unifrac_ex))

#find the average for each sample
df_unifrac_in <- df_unifrac_in[,-1]
df_unifrac_in_sum <- df_unifrac_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_ex <- df_unifrac_ex[,-1]
df_unifrac_ex_sum <- df_unifrac_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_in_sum$value)
sd(df_unifrac_in_sum$value)

mean(df_unifrac_ex_sum$value)
sd(df_unifrac_ex_sum$value)

# C
beta_diversity_calc_inex <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ sample_type, data = sampledf)))
}

beta_dispersion_calc_inex <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$sample_type)
  print(return(permutest(beta)))
}

beta_diversity_calc_inex(physeq_npn)
beta_dispersion_calc_inex(physeq_npn)

#D)

beta_diversity_calc_inex(physeq_npn_euk)
beta_dispersion_calc_inex(physeq_npn_euk)

# Table S3 ----------------------------------------------------------------

# Analysis of the composition of microbiomes (ANCOM) pairwise comparisons for differentially abundant phyla among days. 
# Only significant results at a 0.9 detection level or higher were included.  

# ANCOM function 
# code from https://sites.google.com/site/siddharthamandal1985/research 
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

#create function to pull OTU table for ANCOM analysis, and put it into usable format
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

### bacteria
# add metadata into correct format for ANCOM function
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#blank list to collect results from ANCOM function in
ancom_results_inex <- list()

#run ancom 
otu_ancom <- otu_ancom_make(physeq_npn)
comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "sample_type",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_inex[[1]] <- ancom_sign_taxa


### Eukaryotes
# add metadata into correct format for ANCOM function
metadata_ancom_euk <- metadata_euk
colnames(metadata_ancom_euk)[1] <- 'Sample.ID'

#blank list to collect results from ANCOM function in
ancom_results_inex_euk <- list()

#run ancom 
otu_ancom_euk <- otu_ancom_make(physeq_npn_euk)
comparison_test <- ANCOM.main(OTUdat = otu_ancom_euk,
                              Vardat = metadata_ancom_euk,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "sample_type",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn_euk))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_inex_euk[[1]] <- ancom_sign_taxa

# Table S4 ----------------------------------------------------------------
#interaction among microbiome/mycobiome type (external, internal) and date collected 
# A. & B. alpha-diversity
# C. & D. beta-diversity

# A. 
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

kw_values <- prepare_samples_kw(rich)

# ANOVA

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[1]])
Anova(model,
      type = "II")

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[2]])
Anova(model,
      type = "II")

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[3]])
Anova(model,
      type = "II")

# B. 
erich <- estimate_richness(physeq_npn_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

kw_values <- prepare_samples_kw(rich)

# ANOVA

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[1]])
Anova(model,
      type = "II")

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[2]])
Anova(model,
      type = "II")

model = lm(Observation ~ date_collected:sample_type,
           data = kw_values[[3]])
Anova(model,
      type = "II")

#C.
#PERMANOVA
beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ date_collected * sample_type, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)

#D.
#PERMANOVA
beta_diversity_calc_int(physeq_npn_euk)

# Table S5 ----------------------------------------------------------------
#### core microbiome 
preparing_data_for_core <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata[,c('SampleID', 'sample_type')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$sample_type <- factor(mtotus$sample_type, 
                               levels = c('internal_microbiome', 'external_microbiome'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  new_df <- mtotus %>% group_by(sample_type) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}

#data ready for core analysis
all_otus <- preparing_data_for_core(physeq_npn)

#clean up the data
all_otus <- all_otus %>% filter(total != 0)
all_otus$internal_microbiome <- as.numeric(as.character(all_otus$internal_microbiome))
all_otus$external_microbiome <- as.numeric(as.character(all_otus$external_microbiome))

#separate out by microbiome type
int_otus <- all_otus %>% filter(external_microbiome == 0) %>% 
  filter(internal_microbiome != 0)
ex_otus <- all_otus %>% filter(external_microbiome != 0) %>% 
  filter(internal_microbiome == 0)
core_otus <- all_otus %>% filter(external_microbiome != 0) %>% 
  filter(internal_microbiome != 0)

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names
impTaxa <- core_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
core_phy <- pop_taxa(physeq_npn, impTaxa)
core_phy

#pull out just taxa names
impTaxa <- int_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
int_phy <- pop_taxa(physeq_npn, impTaxa)
int_phy

#pull out just taxa names
impTaxa <- ex_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
ex_phy <- pop_taxa(physeq_npn, impTaxa)
ex_phy

#show the numbers of otus present at certain times. 
find_if_core_increases <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata[,c('SampleID', 'date_collected')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$date_collected <- factor(mtotus$date_collected, 
                                  levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                             '4/23/2019', '7/26/2019', '10/26/2019',
                                             '1/24/2020'))
  new_df <- mtotus %>% group_by(date_collected) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  new_df <- dplyr::select(new_df, -c(OTU)) 
  df_bin <- data.frame(sapply(new_df, function(x) as.numeric(as.character(x))))
  df_bin <- data.frame(sapply(df_bin, function(x) ifelse(x > 0, 1, 0)))
  
  total <- as.vector(colSums(Filter(is.numeric, df_bin)))
  
  #normalize using alpha-div
  erich <- estimate_richness(core_phy, measures = c("Observed"))
  erich <- add_rownames(erich, "SampleID")
  erich_sums <- merge(erich, metadata)
  obs <- erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed'), funs(mean))
  
  date <- colnames(new_df)
  plot_data_df <- data.frame(date,total,obs$Observed)
  plot_data_df$obs.Observed
  plot_data_df$total_norm <- plot_data_df$total / plot_data_df$obs.Observed
  plot_data_df$date <- factor(plot_data_df$date, 
                              levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                         '4/23/2019', '7/26/2019', '10/26/2019',
                                         '1/24/2020'))
  
  return(plot_data_df)
}

#partition number of core OTUs over time
core_df <- find_if_core_increases(core_phy)

#change data type for linear/quadratic models
core_df$date_num <- as.numeric(core_df$date)

#quadratic normalized by alpha-diversity
m1 <- lm(date_num ~ poly(total_norm), data = core_df)
summary(m1)

#linear normalized by alpha-diversity
m2 <- lm(date_num ~ total_norm, data = core_df)
summary(m2)


####best model
#quadratic
m3 <- lm(date_num ~ poly(total), data = core_df)
summary(m3)
#look at slope, increasing at different rates? 
m3

a <- ggplot(core_df, aes(y=total, x=date_num)) +
  geom_point(color = '#3E8E18', size=4) + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab('Number of Core Taxa') + xlab('Collection Date') 
a

core_df$date_num

#linear
m4 <- lm(date_num ~ total, data = core_df)
summary(m4)

#partition number of external OTUs over time
ex_df <- find_if_core_increases(ex_phy)

#change data type for linear/quadratic models
ex_df$date_num <- as.numeric(ex_df$date)

#quadratic normalized by alpha-diversity
m1 <- lm(date_num ~ poly(total_norm), data = ex_df)
summary(m1)

#linear normalized by alpha-diversity
m2 <- lm(date_num ~ total_norm, data = ex_df)
summary(m2)

####best model
#quadratic
m3 <- lm(date_num ~ poly(total), data = ex_df)
summary(m3)
#look at slope, increasing at different rates? 
m3


#linear
#overfit
m4 <- lm(date_num ~ total, data = ex_df)
summary(m4)

#partition number of internal OTUs over time
in_df <- find_if_core_increases(int_phy)

#change data type for linear/quadratic models
in_df$date_num <- as.numeric(in_df$date)

#quadratic normalized by alpha-diversity
m1 <- lm(date_num ~ poly(total_norm), data = in_df)
summary(m1)

#linear normalized by alpha-diversity
m2 <- lm(date_num ~ total_norm, data = in_df)
summary(m2)

####best model
#quadratic
m3 <- lm(date_num ~ poly(total), data = in_df)
summary(m3)
m3



#linear
#overfit
m4 <- lm(date_num ~ total, data = in_df)
summary(m4)


### core mycobiome
preparing_data_for_core <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata_euk[,c('SampleID', 'sample_type')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$sample_type <- factor(mtotus$sample_type, 
                               levels = c('internal_microbiome', 'external_microbiome'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  new_df <- mtotus %>% group_by(sample_type) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}
#data ready for core analysis
all_otus <- preparing_data_for_core(physeq_npn_euk)

#clean up the data
all_otus <- all_otus %>% filter(total != 0)
all_otus$internal_microbiome <- as.numeric(as.character(all_otus$internal_microbiome))
all_otus$external_microbiome <- as.numeric(as.character(all_otus$external_microbiome))

#separate out by microbiome type
int_otus <- all_otus %>% filter(external_microbiome == 0) %>% 
  filter(internal_microbiome != 0)
ex_otus <- all_otus %>% filter(external_microbiome != 0) %>% 
  filter(internal_microbiome == 0)
core_otus <- all_otus %>% filter(external_microbiome != 0) %>% 
  filter(internal_microbiome != 0)

#pull out just taxa names
impTaxa <- core_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
core_phy <- pop_taxa(physeq_npn_euk, impTaxa)
core_phy

#pull out just taxa names
impTaxa <- int_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
int_phy <- pop_taxa(physeq_npn_euk, impTaxa)
int_phy

#pull out just taxa names
impTaxa <- ex_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
ex_phy <- pop_taxa(physeq_npn_euk, impTaxa)
ex_phy

#show the numbers of otus present at certain times. 
find_if_core_increases <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata_euk[,c('SampleID', 'date_collected')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$date_collected <- factor(mtotus$date_collected, 
                                  levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                             '4/23/2019', '7/26/2019', '10/26/2019',
                                             '1/24/2020'))
  new_df <- mtotus %>% group_by(date_collected) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  new_df <- dplyr::select(new_df, -c(OTU)) 
  df_bin <- data.frame(sapply(new_df, function(x) as.numeric(as.character(x))))
  df_bin <- data.frame(sapply(df_bin, function(x) ifelse(x > 0, 1, 0)))
  
  total <- as.vector(colSums(Filter(is.numeric, df_bin)))
  
  #normalize using alpha-div
  erich <- estimate_richness(core_phy, measures = c("Observed"))
  erich <- add_rownames(erich, "SampleID")
  erich_sums <- merge(erich, metadata_euk)
  obs <- erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed'), funs(mean))
  
  date <- colnames(new_df)
  plot_data_df <- data.frame(date,total,obs$Observed)
  plot_data_df$obs.Observed
  plot_data_df$total_norm <- plot_data_df$total / plot_data_df$obs.Observed
  plot_data_df$date <- factor(plot_data_df$date, 
                              levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                         '4/23/2019', '7/26/2019', '10/26/2019',
                                         '1/24/2020'))
  
  return(plot_data_df)
}

#partition number of core OTUs over time
core_df <- find_if_core_increases(core_phy)

#change data type for linear/quadratic models
core_df$date_num <- as.numeric(core_df$date)

#quadratic
m3 <- lm(date_num ~ poly(total), data = core_df)
summary(m3)
#look at slope, increasing at different rates? 
m3

a <- ggplot(ex_df, aes(y=total, x=date_num)) +
  geom_point(color = '#3E8E18', size=4) + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab('Number of Core Taxa') + xlab('Collection Date') 
a

#partition number of external OTUs over time
ex_df <- find_if_core_increases(ex_phy)

#change data type for linear/quadratic models
ex_df$date_num <- as.numeric(ex_df$date)

#quadratic
m3 <- lm(date_num ~ poly(total), data = ex_df)
summary(m3)
#look at slope, increasing at different rates? 
m3

#partition number of internal OTUs over time
in_df <- find_if_core_increases(int_phy)

#change data type for linear/quadratic models
in_df$date_num <- as.numeric(in_df$date)

#quadratic
m3 <- lm(date_num ~ poly(total), data = in_df)
summary(m3)
m3





# Table S6 ----------------------------------------------------------------

pig1 <- subset_samples(physeq_npn, pig == '1')
pig2 <- subset_samples(physeq_npn, pig == '2')
pig3 <- subset_samples(physeq_npn, pig == '3')
pig4 <- subset_samples(physeq_npn, pig == '4')
pig5 <- subset_samples(physeq_npn, pig == '5')

#change this for all of the pigs
#make sure you delete all the empty rows and columns, throws error
meta_mic_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_meta_FEAST_pig1.txt")

otu_mic_feast <- data.frame(otu_table(pig1))
write.table(otu_mic_feast, 'otu_mic_feast_pig5.txt')

#move over column names by one in excel 
otu_mic_feast <- Load_CountMatrix(CountMatrix_path = "otu_mic_feast_pig1.txt")

FEAST_output <- FEAST(C = otu_mic_feast, metadata = meta_mic_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="micro_FEAST_pig1")



PlotSourceContribution(SinkNames = rownames(FEAST_output)[c(5:8)],
                       SourceNames = colnames(FEAST_output), dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                       mixing_proportions = FEAST_output, Plot_title = "Micro_pig1",Same_sources_flag = 0, N = 4)

# Make Boxplots for Source Sink
FEAST_avg <- read.csv("FEAST_ExIn_percent.csv")
FEAST_avg$date <- factor(FEAST_avg$date, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                                      '4/23/2019', '7/26/2019', '10/26/2019',
                                                                      '1/24/2020'))

res.aov2 <- aov(Percent_Source ~ date * SourceSink, data = FEAST_avg)
summary(res.aov2)


##mycobiome
pig1 <- subset_samples(physeq_npn_euk, pig == '1')
pig2 <- subset_samples(physeq_npn_euk, pig == '2')
pig3 <- subset_samples(physeq_npn_euk, pig == '3')
pig4 <- subset_samples(physeq_npn_euk, pig == '4')
pig5 <- subset_samples(physeq_npn_euk, pig == '5')

#change this for all of the pigs
#make sure you delete all the empty rows and columns, throws error
meta_mic_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_metadata_euks_FEAST_pig4.txt")

otu_mic_feast <- data.frame(otu_table(pig4))
write.table(otu_mic_feast, 'otu_mic_feast_euks_pig4.txt')

#move over column names by one in excel 
otu_mic_feast <- Load_CountMatrix(CountMatrix_path = "otu_mic_feast_euks_pig4.txt")

FEAST_output <- FEAST(C = otu_mic_feast, metadata = meta_mic_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="micro_FEAST_euk_pig4")


# Make Boxplots for Source Sink
FEAST_avg <- read.csv("FEAST_ExIn_percent_euk.csv")
FEAST_avg$date <- factor(FEAST_avg$date, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                    '4/23/2019', '7/26/2019', '10/26/2019',
                                                    '1/24/2020'))

res.aov2 <- aov(Percent_Source ~ date * SourceSink, data = FEAST_avg)
summary(res.aov2)

# Table S7 ----------------------------------------------------------------

#A) Microbiome: alpha-diversity interaction among years, seasons, and sample type

erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

kw_values <- prepare_samples_kw(rich)

# ANOVA, including interaction effect
model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[1]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[2]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[3]])
Anova(model,
      type = "II")

#B) Mycobiome: alpha-diversity interaction among years, seasons, and sample type

erich <- estimate_richness(physeq_npn_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

kw_values <- prepare_samples_kw(rich)

# ANOVA, including interaction effect
model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[1]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[2]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year:sample_type,
           data = kw_values[[3]])
Anova(model,
      type = "II")

#C) Microbiome: beta-diversity interaction effect among years, seasons, and sample type
beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Season * Year * sample_type, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)

#D) Mycobiome: beta-diversity interaction effect among years, seasons, and sample type
beta_diversity_calc_int(physeq_npn_euk)


# Table S8 ----------------------------------------------------------------

#microbiomes
##alpha-diversity among years
# A) Mean and standard deviation 
# C) Kruskal-Wallis test

# A
#internal
#calculate diversity metrics
erich <- estimate_richness(physeq_in, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(Year) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Year, data = kw_values[[i]]))
}

# external
#calculate diversity metrics
erich <- estimate_richness(physeq_ex, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(Year) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Year, data = kw_values[[i]]))
}

#mycobiomes
##alpha-diversity among years
# A) Mean and standard deviation 
# C) Kruskal-Wallis test

# A
#internal
#calculate diversity metrics
erich <- estimate_richness(physeq_in_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata_euk)
erich_sums %>% group_by(Year) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Year, data = kw_values[[i]]))
}

# external
#calculate diversity metrics
erich <- estimate_richness(physeq_ex_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata_euk)
erich_sums %>% group_by(Year) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Year, data = kw_values[[i]]))
}



# Table S9 ----------------------------------------------------------------

##beta-diversity among years
# A) Beta-diversity mean/ sd among years
# B) PERMANOVA results for comparing beta-diversity and beta-dispersion 999 permuations. 

# A microbiome
#separate out years to get average beta-diversity 
physeq_year1_in <- subset_samples(physeq_in, Year == '1')
physeq_year2_in <- subset_samples(physeq_in, Year == '2')
physeq_year1_ex <- subset_samples(physeq_ex, Year == '1')
physeq_year2_ex <- subset_samples(physeq_ex, Year == '2')

#calculate unifrac distances 
dis_unifrac_year1_in <- UniFrac(physeq_year1_in, weighted=FALSE, 
                           normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year1_in  <- melt(as.matrix(dis_unifrac_year1_in))

dis_unifrac_year2_in <- UniFrac(physeq_year2_in, weighted=FALSE, 
                             normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year2_in  <- melt(as.matrix(dis_unifrac_year2_in))

dis_unifrac_year1_ex <- UniFrac(physeq_year1_ex, weighted=FALSE, 
                             normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year1_ex  <- melt(as.matrix(dis_unifrac_year1_ex))

dis_unifrac_year2_ex <- UniFrac(physeq_year2_ex, weighted=FALSE, 
                             normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year2_ex  <- melt(as.matrix(dis_unifrac_year2_ex))

#find the average for each sample
df_unifrac_year1_in <- df_unifrac_year1_in[,-1]
df_unifrac_year1_in_sum <- df_unifrac_year1_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year2_in <- df_unifrac_year2_in[,-1]
df_unifrac_year2_in_sum <- df_unifrac_year2_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year1_ex <- df_unifrac_year1_ex[,-1]
df_unifrac_year1_ex_sum <- df_unifrac_year1_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year2_ex <- df_unifrac_year2_ex[,-1]
df_unifrac_year2_ex_sum <- df_unifrac_year2_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

#summarize beta-div. 
mean(df_unifrac_year1_in_sum$value)
sd(df_unifrac_year1_in_sum$value)

mean(df_unifrac_year2_in_sum$value)
sd(df_unifrac_year2_in_sum$value)

mean(df_unifrac_year1_ex_sum$value)
sd(df_unifrac_year1_ex_sum$value)

mean(df_unifrac_year2_ex_sum$value)
sd(df_unifrac_year2_ex_sum$value)

# B mycobiome
#separate out years to get average beta-diversity 
physeq_year1_in <- subset_samples(physeq_in_euk, Year == '1')
physeq_year2_in <- subset_samples(physeq_in_euk, Year == '2')
physeq_year1_ex <- subset_samples(physeq_ex_euk, Year == '1')
physeq_year2_ex <- subset_samples(physeq_ex_euk, Year == '2')

#calculate unifrac distances 
dis_unifrac_year1_in <- UniFrac(physeq_year1_in, weighted=FALSE, 
                                normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year1_in  <- melt(as.matrix(dis_unifrac_year1_in))

dis_unifrac_year2_in <- UniFrac(physeq_year2_in, weighted=FALSE, 
                                normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year2_in  <- melt(as.matrix(dis_unifrac_year2_in))

dis_unifrac_year1_ex <- UniFrac(physeq_year1_ex, weighted=FALSE, 
                                normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year1_ex  <- melt(as.matrix(dis_unifrac_year1_ex))

dis_unifrac_year2_ex <- UniFrac(physeq_year2_ex, weighted=FALSE, 
                                normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year2_ex  <- melt(as.matrix(dis_unifrac_year2_ex))

#find the average for each sample
df_unifrac_year1_in <- df_unifrac_year1_in[,-1]
df_unifrac_year1_in_sum <- df_unifrac_year1_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year2_in <- df_unifrac_year2_in[,-1]
df_unifrac_year2_in_sum <- df_unifrac_year2_in %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year1_ex <- df_unifrac_year1_ex[,-1]
df_unifrac_year1_ex_sum <- df_unifrac_year1_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year2_ex <- df_unifrac_year2_ex[,-1]
df_unifrac_year2_ex_sum <- df_unifrac_year2_ex %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

#summarize beta-div. 
mean(df_unifrac_year1_in_sum$value)
sd(df_unifrac_year1_in_sum$value)

mean(df_unifrac_year2_in_sum$value)
sd(df_unifrac_year2_in_sum$value)

mean(df_unifrac_year1_ex_sum$value)
sd(df_unifrac_year1_ex_sum$value)

mean(df_unifrac_year2_ex_sum$value)
sd(df_unifrac_year2_ex_sum$value)

#C microbiome
beta_diversity_calc_year <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Year, data = sampledf)))
}

beta_dispersion_calc_year <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Year)
  print(return(permutest(beta)))
}

beta_diversity_calc_year(physeq_in)
beta_dispersion_calc_year(physeq_in)

beta_diversity_calc_year(physeq_ex)
beta_dispersion_calc_year(physeq_ex)

#D mycobiome
beta_diversity_calc_year(physeq_in_euk)
beta_dispersion_calc_year(physeq_in_euk)

beta_diversity_calc_year(physeq_ex_euk)
beta_dispersion_calc_year(physeq_ex_euk)

# Table S10 ----------------------------------------------------------------
### bacteria
# add metadata into correct format for ANCOM function
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#create function to pull OTU table for ANCOM analysis, and put it into usable format
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

#blank list to collect results from ANCOM function in
ancom_results_year <- list()

#run ANCOM to find differentially abundant taxa among year
otu_ancom <- otu_ancom_make(physeq_in)
comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "Year",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_year[[1]] <- ancom_sign_taxa

otu_ancom <- otu_ancom_make(physeq_ex)
comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "Year",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_year[[1]] <- ancom_sign_taxa


### Eukaryotes
# add metadata into correct format for ANCOM function
metadata_ancom_euk <- metadata_euk
colnames(metadata_ancom_euk)[1] <- 'Sample.ID'

#blank list to collect results from ANCOM function in
ancom_results_year <- list()

#run ANCOM to find differentially abundant taxa among year
otu_ancom <- otu_ancom_make(physeq_in_euk)
comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom_euk,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "Year",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn_euk))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_year[[1]] <- ancom_sign_taxa

#blank list to collect results from ANCOM function in
ancom_results_year <- list()

#run ANCOM to find differentially abundant taxa among year
otu_ancom <- otu_ancom_make(physeq_ex_euk)
comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom_euk,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "Year",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
w_values$otu.names <- gsub('X', '', w_values$otu.names)
tax <- data.frame(tax_table(physeq_npn_euk))
tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_results_year[[1]] <- ancom_sign_taxa

# Table S11 ---------------------------------------------------------------

#microbiomes
##alpha-diversity among years
# A) Mean and standard deviation 
# C) Kruskal-Wallis test

# A
#internal
#calculate diversity metrics
erich <- estimate_richness(physeq_in, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ date_collected, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$date_collected, 
                                      dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

# external
#calculate diversity metrics
erich <- estimate_richness(physeq_ex, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ date_collected, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$date_collected, 
                                      dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

#mycobiomes
##alpha-diversity among years
# A) Mean and standard deviation 
# C) Kruskal-Wallis test

# A
#internal
#calculate diversity metrics
erich <- estimate_richness(physeq_in_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata_euk)
erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ date_collected, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$date_collected, 
                                      dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

# external
#calculate diversity metrics
erich <- estimate_richness(physeq_ex_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata_euk)
erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}
#list to hold results
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ date_collected, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$date_collected, 
                                      dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}


# Table S12 ---------------------------------------------------------------

##beta-diversity over time (date collected)
# A & B) summary stats
# C & D) PERMANOVA results for comparing beta-diversity and beta-dispersion for 999 permuations. 
# E & F) Pairwise PERMANOVA results for comparing beta-diversity 

#A
#internal
#separate out physeq objects
physeq_time1 <- subset_samples(physeq_in, date_collected == '8/20/2018')
physeq_time2 <- subset_samples(physeq_in, date_collected == '10/25/2018')
physeq_time3 <- subset_samples(physeq_in, date_collected == '1/22/2019')
physeq_time4 <- subset_samples(physeq_in, date_collected == '4/23/2019')
physeq_time5 <- subset_samples(physeq_in, date_collected == '7/26/2019')
physeq_time6 <- subset_samples(physeq_in, date_collected == '10/26/2019')
physeq_time7 <- subset_samples(physeq_in, date_collected == '1/24/2020')

#calculate unifrac distances 
dis_unifrac_1 <- UniFrac(physeq_time1, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_1  <- melt(as.matrix(dis_unifrac_1))

dis_unifrac_2 <- UniFrac(physeq_time2, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_2  <- melt(as.matrix(dis_unifrac_2))

dis_unifrac_3 <- UniFrac(physeq_time3, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_3  <- melt(as.matrix(dis_unifrac_3))

dis_unifrac_4 <- UniFrac(physeq_time4, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_4  <- melt(as.matrix(dis_unifrac_4))

dis_unifrac_5 <- UniFrac(physeq_time5, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_5  <- melt(as.matrix(dis_unifrac_5))

dis_unifrac_6 <- UniFrac(physeq_time6, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_6  <- melt(as.matrix(dis_unifrac_6))

dis_unifrac_7 <- UniFrac(physeq_time7, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_7  <- melt(as.matrix(dis_unifrac_7))

#find the average for each sample
df_unifrac_1 <- df_unifrac_1[,-1]
df_unifrac_1_sum <- df_unifrac_1 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_2 <- df_unifrac_2[,-1]
df_unifrac_2_sum <- df_unifrac_2 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_3 <- df_unifrac_3[,-1]
df_unifrac_3_sum <- df_unifrac_3 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_4 <- df_unifrac_4[,-1]
df_unifrac_4_sum <- df_unifrac_4 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_5 <- df_unifrac_5[,-1]
df_unifrac_5_sum <- df_unifrac_5 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_6 <- df_unifrac_6[,-1]
df_unifrac_6_sum <- df_unifrac_6 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_7 <- df_unifrac_7[,-1]
df_unifrac_7_sum <- df_unifrac_7 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_1_sum$value)
sd(df_unifrac_1_sum$value)

mean(df_unifrac_2_sum$value)
sd(df_unifrac_2_sum$value)

mean(df_unifrac_3_sum$value)
sd(df_unifrac_3_sum$value)

mean(df_unifrac_4_sum$value)
sd(df_unifrac_4_sum$value)

mean(df_unifrac_5_sum$value)
sd(df_unifrac_5_sum$value)

mean(df_unifrac_6_sum$value)
sd(df_unifrac_6_sum$value)

mean(df_unifrac_7_sum$value)
sd(df_unifrac_7_sum$value)

# C 
# internal
beta_diversity_calc_time <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ date_collected, data = sampledf)))
}

beta_dispersion_calc_time <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$date_collected)
  print(return(permutest(beta)))
}

beta_diversity_calc_time(physeq_in)
beta_dispersion_calc_time(physeq_in)

#A 
#external
#separate out physeq objects
physeq_time1 <- subset_samples(physeq_ex, date_collected == '8/20/2018')
physeq_time2 <- subset_samples(physeq_ex, date_collected == '10/25/2018')
physeq_time3 <- subset_samples(physeq_ex, date_collected == '1/22/2019')
physeq_time4 <- subset_samples(physeq_ex, date_collected == '4/23/2019')
physeq_time5 <- subset_samples(physeq_ex, date_collected == '7/26/2019')
physeq_time6 <- subset_samples(physeq_ex, date_collected == '10/26/2019')
physeq_time7 <- subset_samples(physeq_ex, date_collected == '1/24/2020')

#calculate unifrac distances 
dis_unifrac_1 <- UniFrac(physeq_time1, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_1  <- melt(as.matrix(dis_unifrac_1))

dis_unifrac_2 <- UniFrac(physeq_time2, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_2  <- melt(as.matrix(dis_unifrac_2))

dis_unifrac_3 <- UniFrac(physeq_time3, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_3  <- melt(as.matrix(dis_unifrac_3))

dis_unifrac_4 <- UniFrac(physeq_time4, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_4  <- melt(as.matrix(dis_unifrac_4))

dis_unifrac_5 <- UniFrac(physeq_time5, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_5  <- melt(as.matrix(dis_unifrac_5))

dis_unifrac_6 <- UniFrac(physeq_time6, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_6  <- melt(as.matrix(dis_unifrac_6))

dis_unifrac_7 <- UniFrac(physeq_time7, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_7  <- melt(as.matrix(dis_unifrac_7))

#find the average for each sample
df_unifrac_1 <- df_unifrac_1[,-1]
df_unifrac_1_sum <- df_unifrac_1 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_2 <- df_unifrac_2[,-1]
df_unifrac_2_sum <- df_unifrac_2 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_3 <- df_unifrac_3[,-1]
df_unifrac_3_sum <- df_unifrac_3 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_4 <- df_unifrac_4[,-1]
df_unifrac_4_sum <- df_unifrac_4 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_5 <- df_unifrac_5[,-1]
df_unifrac_5_sum <- df_unifrac_5 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_6 <- df_unifrac_6[,-1]
df_unifrac_6_sum <- df_unifrac_6 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_7 <- df_unifrac_7[,-1]
df_unifrac_7_sum <- df_unifrac_7 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_1_sum$value)
sd(df_unifrac_1_sum$value)

mean(df_unifrac_2_sum$value)
sd(df_unifrac_2_sum$value)

mean(df_unifrac_3_sum$value)
sd(df_unifrac_3_sum$value)

mean(df_unifrac_4_sum$value)
sd(df_unifrac_4_sum$value)

mean(df_unifrac_5_sum$value)
sd(df_unifrac_5_sum$value)

mean(df_unifrac_6_sum$value)
sd(df_unifrac_6_sum$value)

mean(df_unifrac_7_sum$value)
sd(df_unifrac_7_sum$value)

# C
#external
beta_diversity_calc_time(physeq_ex)
beta_dispersion_calc_time(physeq_ex)

# B
#Internal
#separate out physeq objects
physeq_time1 <- subset_samples(physeq_in_euk, date_collected == '8/20/2018')
physeq_time2 <- subset_samples(physeq_in_euk, date_collected == '10/25/2018')
physeq_time3 <- subset_samples(physeq_in_euk, date_collected == '1/22/2019')
physeq_time4 <- subset_samples(physeq_in_euk, date_collected == '4/23/2019')
physeq_time5 <- subset_samples(physeq_in_euk, date_collected == '7/26/2019')
physeq_time6 <- subset_samples(physeq_in_euk, date_collected == '10/26/2019')
physeq_time7 <- subset_samples(physeq_in_euk, date_collected == '1/24/2020')

#calculate unifrac distances 
dis_unifrac_1 <- UniFrac(physeq_time1, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_1  <- melt(as.matrix(dis_unifrac_1))

dis_unifrac_2 <- UniFrac(physeq_time2, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_2  <- melt(as.matrix(dis_unifrac_2))

dis_unifrac_3 <- UniFrac(physeq_time3, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_3  <- melt(as.matrix(dis_unifrac_3))

dis_unifrac_4 <- UniFrac(physeq_time4, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_4  <- melt(as.matrix(dis_unifrac_4))

dis_unifrac_5 <- UniFrac(physeq_time5, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_5  <- melt(as.matrix(dis_unifrac_5))

dis_unifrac_6 <- UniFrac(physeq_time6, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_6  <- melt(as.matrix(dis_unifrac_6))

dis_unifrac_7 <- UniFrac(physeq_time7, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_7  <- melt(as.matrix(dis_unifrac_7))

#find the average for each sample
df_unifrac_1 <- df_unifrac_1[,-1]
df_unifrac_1_sum <- df_unifrac_1 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_2 <- df_unifrac_2[,-1]
df_unifrac_2_sum <- df_unifrac_2 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_3 <- df_unifrac_3[,-1]
df_unifrac_3_sum <- df_unifrac_3 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_4 <- df_unifrac_4[,-1]
df_unifrac_4_sum <- df_unifrac_4 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_5 <- df_unifrac_5[,-1]
df_unifrac_5_sum <- df_unifrac_5 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_6 <- df_unifrac_6[,-1]
df_unifrac_6_sum <- df_unifrac_6 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_7 <- df_unifrac_7[,-1]
df_unifrac_7_sum <- df_unifrac_7 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_1_sum$value)
sd(df_unifrac_1_sum$value)

mean(df_unifrac_2_sum$value)
sd(df_unifrac_2_sum$value)

mean(df_unifrac_3_sum$value)
sd(df_unifrac_3_sum$value)

mean(df_unifrac_4_sum$value)
sd(df_unifrac_4_sum$value)

mean(df_unifrac_5_sum$value)
sd(df_unifrac_5_sum$value)

mean(df_unifrac_6_sum$value)
sd(df_unifrac_6_sum$value)

mean(df_unifrac_7_sum$value)
sd(df_unifrac_7_sum$value)

# D
#internal
beta_diversity_calc_time(physeq_in_euk)
beta_dispersion_calc_time(physeq_in_euk)

# B
#External
#separate out physeq objects
physeq_time1 <- subset_samples(physeq_ex_euk, date_collected == '8/20/2018')
physeq_time2 <- subset_samples(physeq_ex_euk, date_collected == '10/25/2018')
physeq_time3 <- subset_samples(physeq_ex_euk, date_collected == '1/22/2019')
physeq_time4 <- subset_samples(physeq_ex_euk, date_collected == '4/23/2019')
physeq_time5 <- subset_samples(physeq_ex_euk, date_collected == '7/26/2019')
physeq_time6 <- subset_samples(physeq_ex_euk, date_collected == '10/26/2019')
physeq_time7 <- subset_samples(physeq_ex_euk, date_collected == '1/24/2020')

#calculate unifrac distances 
dis_unifrac_1 <- UniFrac(physeq_time1, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_1  <- melt(as.matrix(dis_unifrac_1))

dis_unifrac_2 <- UniFrac(physeq_time2, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_2  <- melt(as.matrix(dis_unifrac_2))

dis_unifrac_3 <- UniFrac(physeq_time3, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_3  <- melt(as.matrix(dis_unifrac_3))

dis_unifrac_4 <- UniFrac(physeq_time4, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_4  <- melt(as.matrix(dis_unifrac_4))

dis_unifrac_5 <- UniFrac(physeq_time5, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_5  <- melt(as.matrix(dis_unifrac_5))

dis_unifrac_6 <- UniFrac(physeq_time6, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_6  <- melt(as.matrix(dis_unifrac_6))

dis_unifrac_7 <- UniFrac(physeq_time7, weighted=FALSE, 
                         normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_7  <- melt(as.matrix(dis_unifrac_7))

#find the average for each sample
df_unifrac_1 <- df_unifrac_1[,-1]
df_unifrac_1_sum <- df_unifrac_1 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_2 <- df_unifrac_2[,-1]
df_unifrac_2_sum <- df_unifrac_2 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_3 <- df_unifrac_3[,-1]
df_unifrac_3_sum <- df_unifrac_3 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_4 <- df_unifrac_4[,-1]
df_unifrac_4_sum <- df_unifrac_4 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_5 <- df_unifrac_5[,-1]
df_unifrac_5_sum <- df_unifrac_5 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_6 <- df_unifrac_6[,-1]
df_unifrac_6_sum <- df_unifrac_6 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_7 <- df_unifrac_7[,-1]
df_unifrac_7_sum <- df_unifrac_7 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))


#summarize beta-div. 
mean(df_unifrac_1_sum$value)
sd(df_unifrac_1_sum$value)

mean(df_unifrac_2_sum$value)
sd(df_unifrac_2_sum$value)

mean(df_unifrac_3_sum$value)
sd(df_unifrac_3_sum$value)

mean(df_unifrac_4_sum$value)
sd(df_unifrac_4_sum$value)

mean(df_unifrac_5_sum$value)
sd(df_unifrac_5_sum$value)

mean(df_unifrac_6_sum$value)
sd(df_unifrac_6_sum$value)

mean(df_unifrac_7_sum$value)
sd(df_unifrac_7_sum$value)

# D
#External
beta_diversity_calc_time(physeq_ex_euk)
beta_dispersion_calc_time(physeq_ex_euk)

# E & F
date_pairwise <- function(physeq) {
  physeq_1 <- subset_samples(physeq, date_collected == '8/20/2018')
  physeq_2 <- subset_samples(physeq, date_collected == '10/25/2018')
  physeq_3 <- subset_samples(physeq, date_collected == '1/22/2019')
  physeq_4 <- subset_samples(physeq, date_collected == '4/23/2019')
  physeq_5 <- subset_samples(physeq, date_collected == '7/26/2019')
  physeq_6 <- subset_samples(physeq, date_collected == '10/26/2019')
  physeq_7 <- subset_samples(physeq, date_collected == '1/24/2020')
  return(list(physeq_12 <- merge_phyloseq(physeq_1, physeq_2),
              physeq_13 <- merge_phyloseq(physeq_1, physeq_3),
              physeq_14 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_15 <- merge_phyloseq(physeq_1, physeq_5),
              physeq_16 <- merge_phyloseq(physeq_1, physeq_6),
              physeq_17 <- merge_phyloseq(physeq_1, physeq_7),
              physeq_23 <- merge_phyloseq(physeq_2, physeq_3),
              physeq_24 <- merge_phyloseq(physeq_2, physeq_4),
              physeq_25 <- merge_phyloseq(physeq_2, physeq_5),
              physeq_26 <- merge_phyloseq(physeq_2, physeq_6),
              physeq_27 <- merge_phyloseq(physeq_2, physeq_7),
              physeq_34 <- merge_phyloseq(physeq_3, physeq_4),
              physeq_35 <- merge_phyloseq(physeq_3, physeq_5),
              physeq_36 <- merge_phyloseq(physeq_3, physeq_6),
              physeq_37 <- merge_phyloseq(physeq_3, physeq_7),
              physeq_45 <- merge_phyloseq(physeq_4, physeq_5),
              physeq_46 <- merge_phyloseq(physeq_4, physeq_6),
              physeq_47 <- merge_phyloseq(physeq_4, physeq_7),
              physeq_56 <- merge_phyloseq(physeq_5, physeq_6),
              physeq_57 <- merge_phyloseq(physeq_5, physeq_7),
              physeq_67 <- merge_phyloseq(physeq_6, physeq_7)))
}

# E
# Internal
date_list <- date_pairwise(physeq_in)

for(i in 1:length(date_list)) {
  print(beta_diversity_calc_time(date_list[[i]]))
}

# External
date_list <- date_pairwise(physeq_ex)

for(i in 1:length(date_list)) {
  print(beta_diversity_calc_time(date_list[[i]]))
}

# F
# Internal
date_list <- date_pairwise(physeq_in_euk)

for(i in 1:length(date_list)) {
  print(beta_diversity_calc_time(date_list[[i]]))
}

# External
date_list <- date_pairwise(physeq_ex_euk)

for(i in 1:length(date_list)) {
  print(beta_diversity_calc_time(date_list[[i]]))
}

# Table S13 ---------------------------------------------------------------

## Microbiome
# Internal
#change this for all of the date combos
#make sure you delete all the empty rows and columns, throws error
meta_mic_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_meta_FEAST_timeseries6_int.txt")

#change number for corresponding date combination
otu_mic_feast <- data.frame(otu_table(physeq_in))
write.table(otu_mic_feast, 'otu_mic_feast_time1_in.txt')

#move over column names by one in excel 
otu_mic_feast <- Load_CountMatrix(CountMatrix_path = "otu_mic_feast_time1_in.txt")

FEAST_output <- FEAST(C = otu_mic_feast, metadata = meta_mic_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="micro_FEAST_time6_in")

#External
#change this for all of the date combos
#make sure you delete all the empty rows and columns, throws error
meta_mic_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_meta_FEAST_timeseries6_ex.txt")

#change number for corresponding date combination
otu_mic_feast <- data.frame(otu_table(physeq_ex))
write.table(otu_mic_feast, 'otu_mic_feast_time1_ex.txt')

#move over column names by one in excel 
otu_mic_feast <- Load_CountMatrix(CountMatrix_path = "otu_mic_feast_time1_ex.txt")

FEAST_output <- FEAST(C = otu_mic_feast, metadata = meta_mic_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="micro_FEAST_time6_ex")

## Mycobiome
#Internal
#change this for all of the date combos
#make sure you delete all the empty rows and columns, throws error
meta_myc_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_meta_FEAST_timeseries6_int_euk.txt")

#change number for corresponding date combination
otu_myc_feast <- data.frame(otu_table(physeq_in_euk))
write.table(otu_myc_feast, 'otu_myc_feast_time1_in.txt')

#move over column names by one in excel 
otu_myc_feast <- Load_CountMatrix(CountMatrix_path = "otu_myc_feast_time1_in.txt")

FEAST_output <- FEAST(C = otu_myc_feast, metadata = meta_myc_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="myco_FEAST_time6_in")

#External
#change this for all of the date combos
#make sure you delete all the empty rows and columns, throws error
meta_myc_feast <-Load_metadata(metadata_path = "aquatic_bones_lucas_grant_meta_FEAST_timeseries6_ex_euk.txt")

#change number for corresponding date combination
otu_myc_feast <- data.frame(otu_table(physeq_ex_euk))
write.table(otu_myc_feast, 'otu_myc_feast_time1_ex.txt')

#move over column names by one in excel 
otu_myc_feast <- Load_CountMatrix(CountMatrix_path = "otu_myc_feast_time1_ex.txt")

FEAST_output <- FEAST(C = otu_myc_feast, metadata = meta_myc_feast, 
                      different_sources_flag = 0,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="myco_FEAST_time6_ex")


## reformatted outputs in excel
feast_output <- read.csv('timeseries_FEAST.csv')

feast_op_micro <- feast_output %>% filter(Community == 'Microbiome')

feast_op_micro$Compare <- factor(feast_op_micro$Compare, levels = c('Comp1', 'Other1',
                                                                    'Comp2', 'Other2',
                                                                    'Comp3', 'Other3',
                                                                    'Comp4', 'Other4',
                                                                    'Comp5', 'Other5',
                                                                    'Comp6', 'Other6'
                                                                    
))


feast_op_micro %>% group_by(Compare) %>% summarise_at(c('Contribution'), funs(mean, sd))


feast_op_myco <- feast_output %>% filter(Community == 'Mycobiome')

feast_op_myco$Compare <- factor(feast_op_myco$Compare, levels = c('Comp1', 'Other1',
                                                                  'Comp2', 'Other2',
                                                                  'Comp3', 'Other3',
                                                                  'Comp4', 'Other4',
                                                                  'Comp5', 'Other5',
                                                                  'Comp6', 'Other6'
                                                                  
))


feast_op_myco %>% group_by(Compare) %>% summarise_at(c('Contribution'), funs(mean, sd))


# Table S14 ---------------------------------------------------------------
otu <- as.data.frame(t(otu_table(physeq_npn)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, date_collected)
meta_sa$date_collected <- as.numeric(meta_sa$date_collected)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = date_collected ~ .,
  data    = otu,
  ntree= 500
)

m1

#let's see what the contributing taxa are
imp <- importance(m1)
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(IncNodePurity))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.100 <- imp.sort[1:100, ]
imp.100$predictors<- gsub('X', '', imp.100$predictors)
tax <- data.frame(tax_table(physeq))
tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
tax$predictors <- rownames(tax)
imp.100 <- merge(imp.100, tax)


# Table S15 ---------------------------------------------------------------
# Random forest regression of add using out-of-bag (OOB) classification and 500 decision trees.
# Variation explained, mean square residuals, and mean square error are reported. 

### Microbiome
#first model - all samples OOB
otu <- as.data.frame(t(otu_table(physeq_npn)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, add)
meta_sa$add <- as.numeric(meta_sa$add)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = add ~ .,
  data    = otu,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

#second model - all samples test set error

valid_split <- initial_split(otu, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "add")]
y_test <- otu_valid$add

rf_oob_comp <- randomForest(
  formula = add~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 500
)

oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

rf_oob_comp
sqrt(m1$mse[which.min(rf_oob_comp$mse)])

#decided on OOB over test/train

#function to run random forest
random_forest_date <- function(physeq) {
  otu <- as.data.frame(t(otu_table(physeq)))
  otu$SampleID <- rownames(otu)
  meta_sa <- metadata %>% select(SampleID, add)
  meta_sa$add <- as.numeric(meta_sa$add)
  otu <- merge(meta_sa, otu, by = 'SampleID')
  otu <- otu[,-1]
  names(otu) <- make.names(names(otu))
  
  m1 <- randomForest(
    formula = add ~ .,
    data    = otu,
    ntree= 500
  )
  
  print(m1)
  print(sqrt(m1$mse[which.min(m1$mse)]))
}

#model 3 - only internal microbiome
random_forest_date(physeq_in)

#model 4 - only external microbiome
random_forest_date(physeq_ex)

### Mycobiome
#first model - all samples OOB
otu <- as.data.frame(t(otu_table(physeq_npn_euk)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata_euk %>% select(SampleID, add)
meta_sa$add <- as.numeric(meta_sa$add)
otu <- merge(meta_sa, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = add ~ .,
  data    = otu,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

#second model - all samples test set error

valid_split <- initial_split(otu, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "add")]
y_test <- otu_valid$add

rf_oob_comp <- randomForest(
  formula = add~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 500
)

oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

# compare error rates
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

rf_oob_comp
sqrt(m1$mse[which.min(rf_oob_comp$mse)])

#decided on OOB over test/train

#function to run random forest
random_forest_date <- function(physeq) {
  otu <- as.data.frame(t(otu_table(physeq)))
  otu$SampleID <- rownames(otu)
  meta_sa <- metadata_euk %>% select(SampleID, add)
  meta_sa$add <- as.numeric(meta_sa$add)
  otu <- merge(meta_sa, otu, by = 'SampleID')
  otu <- otu[,-1]
  names(otu) <- make.names(names(otu))
  
  m1 <- randomForest(
    formula = add ~ .,
    data    = otu,
    ntree= 500
  )
  
  print(m1)
  print(sqrt(m1$mse[which.min(m1$mse)]))
}

#model 3 - only internal mycobiome
random_forest_date(physeq_in_euk)

#model 4 - only external mycobiome
random_forest_date(physeq_ex_euk)

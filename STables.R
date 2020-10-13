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
library(reshape2)
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


# Table S1 ----------------------------------------------------------------

##alpha-diversity among pig replicates
# A) Mean and standard deviation 
# B) Kruskal-Wallis test
# C) Pair-wise posthoc Nemenyi test

# A
#calculate diversity metrics
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(pig) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B and C
# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

#each list is a different alpha-div metric
kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ pig, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$pig, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}


# Table S2 ----------------------------------------------------------------

##beta-diversity among pig replicates
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion (999 permuations)
# B) Pair-wise PERMANOVA results for comparing beta-diversity (999 permuations) 

# A
beta_diversity_calc_pig <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ pig, data = sampledf)))
}

beta_dispersion_calc_pig <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$pig)
  print(return(permutest(beta)))
}

beta_diversity_calc_pig(physeq_npn)
beta_dispersion_calc_pig(physeq_npn)

# B
#construct pairwise comparisons 
pig_pairwise <- function(physeq) {
  physeq_1 <- subset_samples(physeq, pig == '1')
  physeq_2 <- subset_samples(physeq, pig == '2')
  physeq_3 <- subset_samples(physeq, pig == '3')
  physeq_4 <- subset_samples(physeq, pig == '4')
  physeq_5 <- subset_samples(physeq, pig == '5')
  return(list(physeq_12 <- merge_phyloseq(physeq_1, physeq_2),
              physeq_13 <- merge_phyloseq(physeq_1, physeq_3),
              physeq_14 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_15 <- merge_phyloseq(physeq_1, physeq_5),
              physeq_23 <- merge_phyloseq(physeq_2, physeq_3),
              physeq_24 <- merge_phyloseq(physeq_2, physeq_4),
              physeq_25 <- merge_phyloseq(physeq_2, physeq_5),
              physeq_34 <- merge_phyloseq(physeq_3, physeq_4),
              physeq_35 <- merge_phyloseq(physeq_3, physeq_5),
              physeq_45 <- merge_phyloseq(physeq_4, physeq_5)))
}

pig_list <- pig_pairwise(physeq_npn)

for(i in 1:length(pig_list)) {
  print(beta_diversity_calc_pig(pig_list[[i]]))
}


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

otu_ancom <- otu_ancom_make(physeq_npn)

# add metadata into correct format for ANCOM function
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#blank list to collect results from ANCOM function in
ancom_results_pig <- list()

#tests ANCOM for each pairwise comparison 
for(i in 1:length(pig_list)) {
  otu_ancom <- otu_ancom_make(pig_list[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "pig",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  w_values$otu.names <- gsub('X', '', w_values$otu.names)
  w_values
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_pig[[i]] <- ancom_sign_taxa
}


# Table S4 ----------------------------------------------------------------

#random forest classification of pig
#using OOB error rate and 1000 decision trees
random_foresting_pig <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$pig)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

m1 <- random_foresting_pig(physeq_npn)
plot(m1)
m1


# Table S5 ----------------------------------------------------------------

##alpha-diversity among years
# A) Mean and standard deviation 
# B) Kruskal-Wallis test

# A
#calculate diversity metrics
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
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

# Table S6 ----------------------------------------------------------------

##beta-diversity among years
# A) Beta-diversity mean/ sd among years
# B) PERMANOVA results for comparing beta-diversity and beta-dispersion 999 permuations. 

# A
#separate out years to get average beta-diversity 
physeq_year1 <- subset_samples(physeq_npn, Year == '1')
physeq_year2 <- subset_samples(physeq_npn, Year == '2')

#calculate unifrac distances 
dis_unifrac_year1 <- UniFrac(physeq_year1, weighted=FALSE, 
                           normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year1  <- melt(as.matrix(dis_unifrac_year1))

dis_unifrac_year2 <- UniFrac(physeq_year2, weighted=FALSE, 
                             normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_year2  <- melt(as.matrix(dis_unifrac_year2))

#find the average for each sample
df_unifrac_year1 <- df_unifrac_year1[,-1]
df_unifrac_year1_sum <- df_unifrac_year1 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_year2 <- df_unifrac_year2[,-1]
df_unifrac_year2_sum <- df_unifrac_year2 %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

#summarize beta-div. 
mean(df_unifrac_year1_sum$value)
sd(df_unifrac_year1_sum$value)

mean(df_unifrac_year2_sum$value)
sd(df_unifrac_year2_sum$value)

# B
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

beta_diversity_calc_year(physeq_npn)
beta_dispersion_calc_year(physeq_npn)



# Table S7 ----------------------------------------------------------------

#blank list to collect results from ANCOM function in
ancom_results_year <- list()

#run ANCOM to find differentially abundant taxa amoung year
otu_ancom <- otu_ancom_make(physeq_npn)
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


# Table S8 ----------------------------------------------------------------

#random forest classification of pig
#using OOB error rate and 1000 decision trees
random_foresting_year <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Year)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

m1 <- random_foresting_year(physeq_npn)
plot(m1)
m1

# Table S9 ----------------------------------------------------------------

##alpha-diversity among seasons (spring, summer, fall, winter)
# A) Mean and standard deviation of alpha diversity metrics
# B) Kruskal-Wallis test among alpha diversity metrics
# C) Pair-wise posthoc Nemenyi test

# A
#organize data
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#calculate summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(Season) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

# B and C

# get data organized and formatted for stats testing
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
              rich_cha <- rich %>% filter(Index == 'Chao1'),
              rich_sha <- rich %>% filter(Index == 'Shannon')))
}

#list containing different metrics 
kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Season, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Season, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

# Table S10 ---------------------------------------------------------------

##beta-diversity among seasons
# A) summary stats
# B) PERMANOVA results for comparing beta-diversity and beta-dispersion for 999 permuations. 
# C) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion for 999 permuations. 

# A
#separate out seasons to get average beta-diversity 
physeq_fall <- subset_samples(physeq_npn, Season == 'Fall')
physeq_spri <- subset_samples(physeq_npn, Season == 'Spring')
physeq_summ <- subset_samples(physeq_npn, Season == 'Summer')
physeq_wint <- subset_samples(physeq_npn, Season == 'Winter')

#calculate unifrac distances 
dis_unifrac_fall <- UniFrac(physeq_fall, weighted=FALSE, 
                             normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_fall  <- melt(as.matrix(dis_unifrac_fall))

dis_unifrac_spri <- UniFrac(physeq_spri, weighted=FALSE, 
                            normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_spri  <- melt(as.matrix(dis_unifrac_spri))

dis_unifrac_summ <- UniFrac(physeq_summ, weighted=FALSE, 
                            normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_summ  <- melt(as.matrix(dis_unifrac_summ))

dis_unifrac_wint <- UniFrac(physeq_wint, weighted=FALSE, 
                            normalized=TRUE, parallel=FALSE, fast=TRUE)
df_unifrac_wint  <- melt(as.matrix(dis_unifrac_wint))

#find the average for each sample
df_unifrac_fall <- df_unifrac_fall[,-1]
df_unifrac_fall_sum <- df_unifrac_fall %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_spri <- df_unifrac_spri[,-1]
df_unifrac_spri_sum <- df_unifrac_spri %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_summ <- df_unifrac_summ[,-1]
df_unifrac_summ_sum <- df_unifrac_summ %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

df_unifrac_wint <- df_unifrac_wint[,-1]
df_unifrac_wint_sum <- df_unifrac_wint %>% group_by(Var2) %>% 
  summarize_at(c('value'), funs(mean))

#summarize beta-div. 
mean(df_unifrac_fall_sum$value)
sd(df_unifrac_fall_sum$value)

mean(df_unifrac_spri_sum$value)
sd(df_unifrac_spri_sum$value)

mean(df_unifrac_summ_sum$value)
sd(df_unifrac_summ_sum$value)

mean(df_unifrac_wint_sum$value)
sd(df_unifrac_wint_sum$value)

# B
beta_diversity_calc_season <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Season, data = sampledf)))
}

beta_dispersion_calc_season <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Season)
  print(return(permutest(beta)))
}

beta_diversity_calc_season(physeq_npn)
beta_dispersion_calc_season(physeq_npn)

# C
season_pairwise <- function(physeq) {
  physeq_F <- subset_samples(physeq, Season == 'Fall')
  physeq_W <- subset_samples(physeq, Season == 'Winter')
  physeq_Sp <- subset_samples(physeq, Season == 'Spring')
  physeq_Su <- subset_samples(physeq, Season == 'Summer')
  return(list(physeq_FW <- merge_phyloseq(physeq_F, physeq_W),
              physeq_FSp <- merge_phyloseq(physeq_F, physeq_Sp),
              physeq_FSu <- merge_phyloseq(physeq_F, physeq_Su),
              physeq_WSp <- merge_phyloseq(physeq_W, physeq_Sp),
              physeq_WSu <- merge_phyloseq(physeq_W, physeq_Su),
              physeq_SpSu <- merge_phyloseq(physeq_Sp, physeq_Su)))
}

season_list <- season_pairwise(physeq_npn)

for(i in 1:length(season_list)) {
  print(beta_diversity_calc_season(season_list[[i]]))
}



# Table S11 ---------------------------------------------------------------
#create function to pull OTU table for ANCOM analysis, and put it into usable format
# for function
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

otu_ancom <- otu_ancom_make(physeq_npn)

# add metadata into correct format for ANCOM function
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#blank list to collect results from ANCOM function in
ancom_results_season <- list()

for(i in 1:length(season_list)) {
  otu_ancom <- otu_ancom_make(season_list[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Season",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  w_values$otu.names <- gsub('X', '', w_values$otu.names)
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_season[[i]] <- ancom_sign_taxa
}


# Table S12 ---------------------------------------------------------------

#random forest classification of pig
#using OOB error rate and 1000 decision trees
random_foresting_season <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Season)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=2000)
  return(Forest)
}

m1 <- random_foresting_season(physeq_npn)
plot(m1)

m1


# Table S13 ---------------------------------------------------------------
#A) alpha-diversity interaction among years and seasons

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
model = lm(Observation ~ Season:Year,
           data = kw_values[[1]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year,
           data = kw_values[[2]])
Anova(model,
      type = "II")

model = lm(Observation ~ Season:Year,
           data = kw_values[[3]])
Anova(model,
      type = "II")


#B) beta-diversity interaction effect among years and seasons 
beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Season * Year, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)


# Table S14 ---------------------------------------------------------------

##alpha-diversity among microbiome types 
# A) Mean and standard deviation of alpha diversity metrics
# B) Kruskal-Wallis test among alpha diversity metrics

# A
#organize data
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(sample_type) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

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

#list to hold alpha div metrics 
kw_values <- prepare_samples_kw(rich)

# KW
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ sample_type, data = kw_values[[i]]))
}


# Table S15 ---------------------------------------------------------------

##beta-diversity among microbiome types (internal, external)
# A) summary stats
# B) PERMANOVA results for comparing beta-diversity and beta-dispersion for 999 permuations. 

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


# Table S16 ---------------------------------------------------------------

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

# Table S17 ---------------------------------------------------------------
#random forest classification of microbiome type (external, internal)
#using OOB error rate and 2000 decision trees
random_foresting_inex <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$sample_type)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=2000)
  return(Forest)
}

m1 <- random_foresting_inex(physeq_npn)
plot(m1)

m1


# Table S18 ---------------------------------------------------------------

#interaction among microbiome type (external, internal) and date collected for
# A. alpha-diversity
# B. beta-diversity

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
#PERMANOVA
beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ date_collected * sample_type, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)

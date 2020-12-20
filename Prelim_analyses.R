#preliminary analyses
#Kaszubinski 2021

#set directory up for success
rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Lucas_grant_project/R_scripts/")

#packages
library(car)
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

#import sequencing information
seq_info <- read.csv("seq_production_stats.csv")

sum(seq_info$PF_Reads)
mean(seq_info$PF_Reads)
sd(seq_info$PF_Reads)

mean(seq_info$R1_Ave_Q)
sd(seq_info$R1_Ave_Q)

mean(seq_info$R2_Ave_Q)
sd(seq_info$R2_Ave_Q)

#import three files: otu table, taxonomy table, and phylo tree
otu <- read.csv("table.csv")
tax <- read.csv("tax_format.csv")
tree <- read_tree('tree.nwk')

#load in metadata
metadata=(read.csv("aquatic_bones_lucas_grant_metadata.csv",header=TRUE))
metadata$pig <- as.factor(metadata$pig)
metadata$Year <- as.factor(metadata$Year)
metadata$date_collected <- factor(metadata$date_collected, levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                                                      '4/23/2019', '7/26/2019', '10/26/2019',
                                                                      '1/24/2020'))
#format it into phyloseq format
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

#put otu table in phyloseq format
rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
#merge tree and otu tables
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

#put taxonomy table in phyloseq format
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


#determine what to rarefy to
#show rarefaction curves
tiff("SFIG1.TIF", width = 2000, height = 1500, res=300)
rarecurve(t(otu_table(physeq_trim)), step=1000, cex=0.3)
dev.off()

#rarefy
physeq_7000 <- rarefy_even_depth(physeq_trim, sample.size = 7000)
physeq_7000
#2475, 65 samples

#removes positive and negative control samples
physeq_in <- subset_samples(physeq_7000, sample_type == 'internal_microbiome')
physeq_ex <- subset_samples(physeq_7000, sample_type == 'external_microbiome')

physeq_npn <- merge_phyloseq(physeq_in, physeq_ex)

###Analyses
#Q1: Question 1: what variation are we capturing in this study, 
#and what variation is useful for forensic practitioners should consider? 


# PIG ---------------------------------------------------------------------


#test and see if the pigs differ from each other
#alpha, beta, ANCOM, and RF

##alpha-diversity
# A) Mean and standard deviation of alpha diversity metrics across pigs.
# B) Kruskal-Wallis test among alpha diversity metrics across pigs.
# C) Pair-wise posthoc Nemenyi test

# A
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
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

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ pig, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$pig, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

##beta-diversity
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 

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

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="pig", shape = 'sample_type') +
  scale_color_manual(values = c('#2F6E2B', '#E37D3E', '#886F98', '#E4BE18', '#A41400')) +
  scale_shape_manual(values=c(17,8)) + geom_point(size = 5) + 
  labs(shape = "Sample Type")

ordplot

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
ancom_results_pig <- list()

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
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_pig[[i]] <- ancom_sign_taxa
}


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







# Year --------------------------------------------------------------------


#test and see if years differ from each other
#alpha, beta, ANCOM, and RF

##alpha-diversity
# A) Mean and standard deviation of alpha diversity metrics across years.
# B) Kruskal-Wallis test among alpha diversity metrics across years.
# C) Pair-wise posthoc Nemenyi test

# A
erich <- estimate_richness(physeq_7000, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(Year) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

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

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Year, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Year, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

##beta-diversity
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 

# A
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


#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="Year", shape = 'sample_type') +
  scale_color_manual(values = c('#2F6E2B', '#E37D3E')) +
  scale_shape_manual(values=c(17,8)) + geom_point(size = 5) + 
  labs(shape = "Sample Type")

ordplot

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
ancom_results_year <- list()

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
  tax <- data.frame(tax_table(physeq_npn))
  tax <- tax %>% select(Kingdom, Phylum, Order, Class, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_year[[1]] <- ancom_sign_taxa



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


# Season ------------------------------------------------------------------

#test and see if the pigs differ from each other
#alpha, beta, ANCOM, and RF

##alpha-diversity
# A) Mean and standard deviation of alpha diversity metrics across pigs.
# B) Kruskal-Wallis test among alpha diversity metrics across pigs.
# C) Pair-wise posthoc Nemenyi test

# A
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
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

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Season, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Season, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

##beta-diversity
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 

# A
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

# B
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

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="Season", shape = 'sample_type') +
  scale_color_manual(values = c('#2F6E2B', '#E37D3E', '#886F98', '#E4BE18')) +
  scale_shape_manual(values=c(17,8)) + geom_point(size = 5) + 
  labs(shape = "Sample Type")

ordplot

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








# Year and Season Interaction ---------------------------------------------

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




beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Season * Year, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)










# In vs. Ex Microbiome ----------------------------------------------------


#Q2: Question 2: Was our hypothesis founded, did the internal microbiome 
#resemble the external microbiome over time?

##alpha-diversity
# A) Mean and standard deviation of alpha diversity metrics across pigs.
# B) Kruskal-Wallis test among alpha diversity metrics across pigs.
# C) Pair-wise posthoc Nemenyi test

# A
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(sample_type) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

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

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ sample_type, data = kw_values[[i]]))
}

##beta-diversity
# A) PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 
# B) Pair-wise PERMANOVA results for comparing beta-diversity and beta-dispersion among pigs for 999 permuations. 

# A
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


beta_diversity_calc_inex_int1 <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ sample_type * Season, data = sampledf)))
}

beta_diversity_calc_inex_int1(physeq_npn)


beta_diversity_calc_inex_int2 <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ sample_type * date_collected, data = sampledf)))
}

beta_diversity_calc_inex_int2(physeq_npn)


#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
cols = c('#32C428', '#28C45B', '#28C486', '#28B2C4', '#2862C4', "#2B28C4", "#6428C4")
ordplot=plot_ordination(physeq_npn, ord, shape="sample_type", color = 'date_collected') +
  scale_shape_manual(values=c(8, 17)) + geom_point(size = 5) +
  scale_color_manual(values = cols) +
  labs(shape = "Sample Type", color = 'Date Collected') 

ordplot


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
  #change based on MOD/COD
  new_df <- mtotus %>% group_by(sample_type) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}

#start with MoD
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

write_phyloseq(core_phy, 'TAXONOMY')

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

#show the numbers of otus present at certain times. If there is an increase,
# then the core microbiome is growing over time. if not, no trend
#with int and ext, should be decreasing over time. 

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

core_df <- find_if_core_increases(core_phy)
ggplot(core_df, aes(x=date, y=total_norm)) + geom_point() 
core_df$date_num <- as.numeric(core_df$date)
m1 <- lm(date_num ~ poly(total_norm), data = core_df)
summary(m1)
m2 <- lm(date_num ~ total_norm, data = core_df)
summary(m2)
####best model
#look at slope, increasing at different rates? 
m3 <- lm(date_num ~ poly(total), data = core_df)
summary(m3)
m3
m4 <- lm(date_num ~ total, data = core_df)
summary(m4)

ex_df <- find_if_core_increases(ex_phy)
ggplot(ex_df, aes(x=date, y=total_norm)) + geom_point() 
ex_df$date_num <- as.numeric(ex_df$date)
m1 <- lm(date_num ~ poly(total_norm), data = ex_df)
summary(m1)
m2 <- lm(date_num ~ total_norm, data = ex_df)
summary(m2)
m3 <- lm(date_num ~ poly(total), data = ex_df)
summary(m3)
m3
m4 <- lm(date_num ~ total, data = ex_df)
summary(m4)

in_df <- find_if_core_increases(int_phy)
ggplot(in_df, aes(x=date, y=total_norm)) + geom_point()  
in_df$date_num <- as.numeric(in_df$date)
m1 <- lm(date_num ~ poly(total_norm), data = in_df)
summary(m1)
m2 <- lm(date_num ~ total_norm, data = in_df)
summary(m2)
m3 <- lm(date_num ~ poly(total), data = in_df)
summary(m3)
m3
m4 <- lm(date_num ~ total, data = in_df)
summary(m4)


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
ancom_results_inex <- list()

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

#random forest classification of pig
#using OOB error rate and 1000 decision trees
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

forest_predictors_inex <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.100 <- imp.sort[1:500, ]
  imp.100$predictors <- gsub('X', '', imp.100$predictors)
  tax <- data.frame(tax_table(physeq_npn))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.100 <- merge(imp.100, tax)
  return(imp.100)
}

for_pred_inex <- forest_predictors_inex(m1)
inex_pred_phy <- pop_taxa(physeq_npn, for_pred_inex$predictors)
inex_pred_phy

find_rf_preds_time <- function(physeq) {
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
  
  plot_data_df <- data.frame(date,total)
  plot_data_df$date <- factor(plot_data_df$date, 
                              levels = c('8/20/2018', '10/25/2018', '1/22/2019',
                                         '4/23/2019', '7/26/2019', '10/26/2019',
                                         '1/24/2020'))
  
  return(plot_data_df)
}

rf_pred_df <- find_rf_preds_time(inex_pred_phy)
ggplot(rf_pred_df, aes(x=date, y=total)) + geom_point() 
rf_pred_df$date_num <- as.numeric(rf_pred_df$date)
m1 <- lm(date_num ~ poly(total), data = rf_pred_df)
summary(m1)
m1
#essentially no difference over time 
m2 <- lm(date_num ~ total, data = rf_pred_df)
summary(m2)
m2

# Sample type and Time interaction ----------------------------------------

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




beta_diversity_calc_int <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ date_collected * sample_type, data = sampledf)))
}

beta_diversity_calc_int(physeq_npn)



# pmsi modeling -----------------------------------------------------------

#let's looks at alpha and beta-diversity
# A
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
erich_sums <- merge(erich, metadata)
erich_sums %>% group_by(date_collected) %>% summarise_at(c('Observed', 'Chao1', "Shannon"), funs(mean, sd))

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

kw_values <- prepare_samples_kw(rich)

# KW and post-hoc Nemenyi 
for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ date_collected, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$date_collected, 
                                      dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

#try logistic regression with alpha div values
kw_values[[1]]$date_num <- as.numeric(kw_values[[1]]$date_collected)
m1 <- lm(kw_values[[1]]$date_num ~ poly(kw_values[[1]]$Observation))
summary(m1)

kw_values[[2]]$date_num <- as.numeric(kw_values[[2]]$date_collected)
m2 <- lm(kw_values[[2]]$date_num ~ poly(kw_values[[2]]$Observation))
summary(m2)

kw_values[[3]]$date_num <- as.numeric(kw_values[[3]]$date_collected)
m3 <- lm(kw_values[[3]]$date_num ~ poly(kw_values[[3]]$Observation))
summary(m3)


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

beta_diversity_calc_time(physeq_npn)
beta_dispersion_calc_time(physeq_npn)

# B
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

date_list <- date_pairwise(physeq_npn)

for(i in 1:length(date_list)) {
  print(beta_diversity_calc_time(date_list[[i]]))
}




# Random forest regression of ADH using out-of-bag (OOB) classification and 500 decision trees.
# Variation explained, mean square residuals, and mean squre error are reported. 

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
sqrt(m1$mse[which.min(m1$mse)])

#testing test set error rate. it ended up being slightly higher than
#the OOB error. 500 trees seems legit 
valid_split <- initial_split(otu, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "date_collected")]
y_test <- otu_valid$date_collected

rf_oob_comp <- randomForest(
  formula = date_collected~ .,
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


#pretty accurate! 

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
sqrt(m1$mse[which.min(m1$mse)])

pred <- m1$predicted
pred <- as.data.frame(pred)
colnames(pred) <- 'Predicted'

#using the mean square error to determine the true value
pred[pred$Predicted < 7.851, "True"] <- '7'
pred[pred$Predicted < 6.851, "True"] <- '6'
pred[pred$Predicted < 5.851, "True"] <- '5'
pred[pred$Predicted < 4.851, "True"] <- '4'
pred[pred$Predicted < 3.851, "True"] <- '3'
pred[pred$Predicted < 2.851, "True"] <- '2'
pred[pred$Predicted < 1.851, "True"] <- '1'

pred$True <- as.numeric(pred$True)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") 
}

d <- ggplotRegression(lm(Predicted ~ True, data = pred)) +
  geom_pointrange(aes(ymin = Predicted-0.851, ymax = Predicted+0.851), 
                  position=position_jitter(width=0.5), alpha = 0.5) 

lm_eqn <- function(df){
  m <- lm(Predicted ~ True, df);
  eq <- substitute(italic(Predicted) == a + b %.% italic(True)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

d1 <- d + geom_text(x = 2, y = 6, label = lm_eqn(pred), parse = TRUE)
d1


#things to do: take out pig 1, try just one year, try with one 
#community type

random_forest_date <- function(physeq) {
  otu <- as.data.frame(t(otu_table(physeq)))
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
  
  print(m1)
  print(sqrt(m1$mse[which.min(m1$mse)]))
}

#try it without pig 1
physeq_nopig1 <- subset_samples(physeq_npn, pig != '1')
random_forest_date(physeq_nopig1)

#different community types

#internal microbiome
random_forest_date(physeq_in)

#external microbiome
random_forest_date(physeq_ex)

#year 1
physeq_year1 <- subset_samples(physeq_npn, Year == '1')
random_forest_date(physeq_year1)

#year 2
physeq_year2 <- subset_samples(physeq_npn, Year == '2')
random_forest_date(physeq_year2)



# ADH ---------------------------------------------------------------------

add <- read.csv('meantemp_ADH_Lucas.csv')
add_sort <- add %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d")) %>%
  arrange(Date)

add_model1 <- lm(X ~ CumulativeTemp, data = add)
summary(add_model1)

add_model2 <- lm(X ~ poly(CumulativeTemp), data = add)
summary(add_model2)

a <- ggplot(add, aes(y=X, x=CumulativeTemp)) +
  geom_point(color = '#3E8E18', size=1) + stat_smooth(method = "lm", formula = y ~ x + I(x^3), size = 1) +
  ylab('Day Since Placement') + xlab('Accumulated Degree Days (ADD)')
a

#random forest
otu <- as.data.frame(t(otu_table(physeq_npn)))
otu$SampleID <- rownames(otu)
meta_sa <- metadata %>% select(SampleID, date_collected)

#ADD instead of date
meta_sa_sort <- meta_sa %>%
  mutate(date_collected = as.Date(date_collected, "%m/%d/%Y")) %>%
  arrange(date_collected)
meta_sa_sort$Date <- meta_sa_sort$date_collected

meta_sa_fin <- merge(add_sort, meta_sa_sort, by = 'Date' )
meta_sa_fin <- meta_sa_fin[,4:5]

otu <- merge(meta_sa_fin, otu, by = 'SampleID')
otu <- otu[,-1]
names(otu) <- make.names(names(otu))

m1 <- randomForest(
  formula = CumulativeTemp ~ .,
  data    = otu,
  ntree= 500
)

m1
sqrt(m1$mse[which.min(m1$mse)])

            
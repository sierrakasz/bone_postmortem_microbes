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
#removes positive and negative control samples
physeq_in <- subset_samples(physeq_7000, sample_type == 'internal_microbiome')
physeq_ex <- subset_samples(physeq_7000, sample_type == 'external_microbiome')
physeq_npn <- merge_phyloseq(physeq_in, physeq_ex)

#calculate beta div with unifrac dissimilarity matrix
ord_in = ordinate(physeq_in, method="PCoA", distance="unifrac")
ord_ex = ordinate(physeq_ex, method="PCoA", distance="unifrac")
#create plot
cols <- c('#C88D39', '#5EA849', '#4972A8', '#8149A8', '#A84949')

ordplot_ex=plot_ordination(physeq_ex, ord_ex, color="pig") +
  scale_color_manual(values = cols) +
  labs(color = 'Pig Replicate') + geom_point(shape = 17, size = 5,  alpha = 0.7)

ordplot_ex

ordplot_in=plot_ordination(physeq_in, ord_in, color="pig") +
  scale_color_manual(values = cols) +
  labs(color = 'Pig Replicate') + geom_point(shape = 16, size = 5, alpha = 0.7)

ordplot_in

theme_set(theme_classic(base_size = 12))
tiff("SFIG2.TIF", width = 3500, height = 2000, res=300)
ggarrange(ordplot_ex, ordplot_in, 
          labels = c("A.", "B."),
          nrow = 2, ncol = 2)
dev.off()


# Figure S3 ---------------------------------------------------------------

#first panel; alpha-diversity amoung years

#calculate diversity metrics
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)

#organize data
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

#make the plot
p <- ggplot(rich, aes(x=Year, y=Observation, fill=Year)) +
  geom_boxplot(outlier.size = 3)
a <- p + facet_wrap(~Index, scales="free") + 
  scale_fill_manual(values = c("#3984C8", "#C83953")) +
  scale_x_discrete(labels=c('1' = 'Year 1', '2' = 'Year 2')) +
  xlab("") + ylab('Counts') + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom")

a

#second panel; beta-diversity amoung years

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="Year", shape = 'sample_type') +
  scale_color_manual(values = c("#3984C8", "#C83953")) +
  scale_shape_manual(values=c(17,16)) + geom_point(alpha = 0.5, size = 5) + 
  labs(shape = "Sample Type")

ordplot

theme_set(theme_classic(base_size = 16))
tiff("SFIG3.TIF", width = 2500, height = 3500, res=300)
ggarrange(a, ordplot, 
          labels = c("A.", "B."),
          nrow = 2, ncol = 1)
dev.off()


# Figure S4  ----------------------------------------------------------------
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

ancom_taxa_year <- filter(w_values, detected_0.9 == TRUE)

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

phy_ancom_year <- pop_taxa(physeq_npn, ancom_taxa_year$otu.names)

#to make things easier, going to plot phyla of ANCOM taxa
GPrPhylum=tax_glom(phy_ancom_year, "Phylum")
#make abundances out of 100
PhylumLevel_year = transform_sample_counts(GPrPhylum, function(x) x / sum(x))

#get abundances into format that can be graphed
df <- psmelt(PhylumLevel_year) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Year"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Trtdata <-Trtdata[order(Trtdata$Year, Trtdata$mean),]
Trtdata$Phylum <- factor(Trtdata$Phylum, 
                         levels = c('Planctomycetes', 'Elusimicrobia', 'Firmicutes',
                                    'Unassigned', 'Kiritimatiellaeota', 'WS2', 'FBP',
                                    'Proteobacteria','Bacteroidetes', 'Acidobacteria', 
                                    'Patescibacteria'))


theme_set(theme_classic(base_size = 16))
tiff("SFIG4.TIF", width = 2500, height = 2000, res=300)
ggplot(Trtdata, aes(x = Year, y = Phylum)) + 
  geom_point(aes(color = Year, size = mean), alpha = 0.5) +
  scale_color_manual(values = c("#3984C8", "#C83953")) +
  scale_size(range = c(1, 12)) + 
  labs(size = "Relative Abundance") +
  ylab('Phyla')
dev.off()



# Figure S5 ---------------------------------------------------------------

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

ancom_taxa_inex <- filter(w_values, detected_0.9 == TRUE)

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

phy_ancom_inex <- pop_taxa(physeq_npn, ancom_taxa_inex$otu.names)

#to make things easier, going to plot phyla of ANCOM taxa
GPrPhylum=tax_glom(phy_ancom_inex, "Phylum")
#make abundances out of 100
PhylumLevel_inex = transform_sample_counts(GPrPhylum, function(x) x / sum(x))

#get abundances into format that can be graphed
df <- psmelt(PhylumLevel_inex) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "sample_type", "date_collected"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

#going to choose top 9 taxa, there are too many to show in a graph
top_9 <- c('Lentisphaerae', 'Elusimicrobia', 'Chlamydiae', 'Proteobacteria',
           'Latescibacteria', 'Firmicutes', 'Chloroflexi', 'Planctomycetes', 
           'Verrucomicrobia')

Trtdata_filter <- subset(Trtdata, Phylum %in% top_9)

theme_set(theme_classic(base_size = 16))
tiff("SFIG5.TIF", width = 2500, height = 3000, res=300)
ggplot(Trtdata_filter, aes(x = sample_type, y = Phylum)) + 
  geom_point(aes(color = sample_type, size = mean), alpha = 0.5) +
  scale_color_manual(values = c("#3984C8", "#C83953")) +
  scale_size(range = c(1, 12)) + facet_wrap(~ date_collected) +
  labs(size = "Relative Abundance", color = 'Microbial Communities') +
  ylab('Phyla') + xlab('') + theme(axis.ticks.x=element_blank(),
                                                        axis.text.x=element_blank())
dev.off()


# Figure S6 ---------------------------------------------------------------

#first panel; alpha-diversity amoung years

#calculate diversity metrics
erich <- estimate_richness(physeq_npn, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata)

#organize data
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata)

#make the plot
p <- ggplot(rich, aes(x=date_collected, y=Observation, fill=date_collected)) +
  geom_boxplot(outlier.size = 3)
a <- p + facet_grid(Index~sample_type, scales = 'free') +
  scale_fill_manual(values = c("#E38F0F", "#CE0900", "#2F94AE", "#2FAE49",
                               "#E3C90F", "#AE2F4C", "#562FAE")) +
  scale_x_discrete(labels=c('1' = 'Year 1', '2' = 'Year 2')) +
  xlab("") + ylab('Alpha-diversity') + labs(fill = 'Date Collected') + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "right", axis.ticks.x=element_blank(),
  axis.text.x=element_blank())

a

#second panel; beta-diversity amoung years

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="date_collected", shape = 'sample_type') +
  scale_color_manual(values =  c("#E38F0F", "#CE0900", "#2F94AE", "#2FAE49",
                                 "#E3C90F", "#AE2F4C", "#562FAE")) +
  scale_shape_manual(values=c(17,16)) + geom_point(alpha = 0.5, size = 5) + 
  labs(color = 'Date Collected', shape = "Sample Type")

ordplot

theme_set(theme_classic(base_size = 16))
tiff("SFIG6.TIF", width = 2500, height = 3500, res=300)
ggarrange(a, ordplot, 
          labels = c("A.", "B."),
          nrow = 2, ncol = 1)
dev.off()


# Figure S7 ---------------------------------------------------------------

x <- transform(dietswap, "compositional")
otu <- abundances(x)
metadata <- meta(x)

rda.result <- vegan::rda(t(otu) ~ factor(metadata$nationality),
                         na.action = na.fail, scale = TRUE)

plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = metadata$time)
points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = metadata$time)
pl <- ordihull(rda.result, metadata$nationality, scaling = 3, label = TRUE)


data_pca <- transform(physeq_npn, "compositional")
otu_pca <- abundances(physeq_npn)
meta_pca <- meta(physeq_npn)

rda.res <- vegan::rda(t(otu_pca) ~ factor(meta_pca$date_collected) + factor(meta_pca$sample_type))
screeplot(rda.res)

R2 <- RsquareAdj(rda.res)$r.squared
R2 

plot(rda.res, scaling = 1)
spe.sc <- scores(otu_pca, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

ggord(rda.res, meta_pca$date_collected) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

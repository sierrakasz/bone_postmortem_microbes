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

# Figure 2 --------------------------------------------------------------
#first panel; alpha-diversity amoung microbial types (interal, external)

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
p <- ggplot(rich, aes(x=sample_type, y=Observation, fill = sample_type)) +
  geom_violin(color = 'white') + geom_jitter(shape=16, position=position_jitter(0.2), color = 'black')
a <- p + facet_wrap(~Index, scales = 'free') +
  scale_fill_manual(values = c("#3958BD", "#4EC8E1"), labels = c("External Microbiome", 
                                                                 "Internal Microbiome")) +
  scale_x_discrete(labels=c('internal_microbiome' = 'Internal', 'external_microbiome' = 'External')) + 
  xlab("") + ylab('Alpha-diversity') + labs(fill = 'Microbiome Type') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "left", axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

a

#second panel; beta-diversity amoung microbiome types

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn, method="PCoA", distance="unifrac")
#create plot
ordplot=plot_ordination(physeq_npn, ord, color="sample_type") +
  scale_color_manual(values =  c("#4EC8E1", "#3958BD")) + 
  geom_point(size = 5) + 
  labs(color = 'Microbiome Type') + 
  stat_ellipse(alpha = 1, size=2, aes(color= sample_type)) + 
  theme(legend.position = "none")

ordplot

#third panel; ancom results amoung microbiome types

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
Trtdata <- ddply(df, c("Phylum", "sample_type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

#going to choose top 12 taxa, there are too many to show in a graph
Trtdata_filter <- subset(Trtdata, mean > 1)

c <- ggplot(Trtdata_filter, aes(x = sample_type, y = Phylum)) + 
  geom_point(aes(color = sample_type, size = mean), alpha = 0.5) +
  scale_color_manual(values = c("#4EC8E1", "#3958BD"), labels = c('External Microbiome',
                                                                  'Internal Microbiome')) +
  scale_size(range = c(1, 12)) +
  labs(size = "Relative Abundance", color = 'Microbial Communities') +
  ylab('Phyla') + xlab('') + theme(axis.ticks.x=element_blank(),
                                   axis.text.x=element_blank())  + 
  theme(legend.position = "right")
c

#fourth panel; alpha-diversity metrics eukaryotes
#calculate diversity metrics
erich <- estimate_richness(physeq_npn_euk, measures = c("Observed", 'Chao1', "Shannon"))
erich <- add_rownames(erich, "SampleID")
#summary stats
erich_sums <- merge(erich, metadata_euk)

#organize data
erich <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon"), na.rm = TRUE)
rich = merge(erich, metadata_euk)

#make the plot
p <- ggplot(rich, aes(x=sample_type, y=Observation, fill = sample_type)) +
  geom_violin(color = 'white') + geom_jitter(shape=16, position=position_jitter(0.2), color = 'black')
d <- p + facet_wrap(~Index, scales = 'free') +
  scale_fill_manual(values = c("#00AA54", "#F1C138"), labels = c("External Communities", 
                                                               "Internal Communities")) +
  scale_x_discrete(labels=c('internal_microbiome' = 'Internal', 'external_microbiome' = 'External')) + 
  xlab("") + ylab('Alpha-diversity') + labs(fill = 'Mycobiome Type') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "left", axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

d

#fifth panel; beta-diversity amoung eukaryotes

#calculate beta div with unifrac dissimilarity matrix
ord = ordinate(physeq_npn_euk, method="PCoA", distance="unifrac")
#create plot
ordplot_euk=plot_ordination(physeq_npn_euk, ord, color="sample_type") +
  scale_color_manual(values =  c("#00AA54", "#F1C138")) + 
  geom_point(size = 5) + 
  labs(color = 'Microbiome Type') + 
  stat_ellipse(alpha = 1, size=2, aes(color= sample_type)) + 
  theme(legend.position = "none")

ordplot_euk

theme_set(theme_classic(base_size = 16))
tiff("FIG2.TIF", width = 6500, height = 3000, res=300)
ggarrange(a, ordplot, c, d, ordplot_euk, 
          labels = c("A.", "B.", "C.", "D.", "E."),
          nrow = 2, ncol = 3)
dev.off()


# Figure 3 ----------------------------------------------------------------
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

#partition number of internal OTUs over time
in_df <- find_if_core_increases(int_phy)

#change data type for linear/quadratic models
in_df$date_num <- as.numeric(in_df$date)

a <- ggplot(in_df, aes(y=total, x=date_num)) +
  geom_point(color = '#4EC8E1', size=4) + geom_smooth(method = "lm", 
                                                      formula = y ~ x + I(x^2), 
                                                      size = 1, color = 'Black') +
  ylab('Number of Internal Taxa') + xlab('Collection Date') 
a

#partition number of core OTUs over time
core_df <- find_if_core_increases(core_phy)

#change data type for linear/quadratic models
core_df$date_num <- as.numeric(core_df$date)

b <- ggplot(core_df, aes(y=total, x=date_num)) +
  geom_point(color = '#D172C7', size=4) + geom_smooth(method = "lm", 
                                                      formula = y ~ x + I(x^2), 
                                                      size = 1, color = 'Black') +
  ylab('Number of Core Taxa') + xlab('Collection Date') 
b

#partition number of external OTUs over time
ex_df <- find_if_core_increases(ex_phy)

#change data type for linear/quadratic models
ex_df$date_num <- as.numeric(ex_df$date)

c <- ggplot(ex_df, aes(y=total, x=date_num)) +
  geom_point(color = "#3958BD", size=4) + geom_smooth(method = "lm", 
                                                      formula = y ~ x + I(x^2), 
                                                      size = 1, color = 'Black') +
  ylab('Number of External Taxa') + xlab('Collection Date') 
c

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

#pull out just taxa names
impTaxa <- int_otus$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
int_phy <- pop_taxa(physeq_npn_euk, impTaxa)
int_phy

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

#partition number of internal OTUs over time
in_df <- find_if_core_increases(int_phy)

#change data type for linear/quadratic models
in_df$date_num <- as.numeric(in_df$date)

d <- ggplot(in_df, aes(y=total, x=date_num)) +
  geom_point(color = "#F1C138", size=4) + geom_smooth(method = "lm", 
                                                      formula = y ~ x + I(x^2), 
                                                      size = 1, color = 'Black') +
  ylab('Number of Internal Taxa') + xlab('Collection Date') 
d

theme_set(theme_classic(base_size = 14))
tiff("FIG3.TIF", width = 4000, height = 2000, res=300)
ggarrange(a,b,c,d, 
          labels = c("A.", "B.", 'C.', 'D.'),
          nrow = 2, ncol = 3)
dev.off()


# Figure 4 ----------------------------------------------------------------

pig1 <- subset_samples(physeq_npn, pig == '1')
pig2 <- subset_samples(physeq_npn, pig == '2')
pig3 <- subset_samples(physeq_npn, pig == '3')
pig4 <- subset_samples(physeq_npn, pig == '4')
pig5 <- subset_samples(physeq_npn, pig == '5')

#changes these for all the pigs 1-5
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

FEAST_avg$Percent_Source <- FEAST_avg$Percent_Source * 100

p1 <- ggline(FEAST_avg, x='date', y='Percent_Source', color = 'SourceSink', 
            add = c("mean_se", "dotplot"),
            palette = c("#00AFBB", "#E7B800")) + xlab('Collection Date') + 
            ylab('Percent Source Contribution (%)') +
            labs(col = 'Source') 

p1

##mycobiome
pig1 <- subset_samples(physeq_npn_euk, pig == '1')
pig2 <- subset_samples(physeq_npn_euk, pig == '2')
pig3 <- subset_samples(physeq_npn_euk, pig == '3')
pig4 <- subset_samples(physeq_npn_euk, pig == '4')
pig5 <- subset_samples(physeq_npn_euk, pig == '5')

#changes these for all the pigs 1-5
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

FEAST_avg$Percent_Source <- FEAST_avg$Percent_Source * 100

p2 <- ggline(FEAST_avg, x='date', y='Percent_Source', color = 'SourceSink',
            add = c("mean_se", "dotplot"),
            palette = c("#78B474", "#86759B")) + xlab('Collection Date') + 
            ylab('Percent Source Contribution (%)') +
            labs(col = 'Source')  

p2

theme_set(theme_classic(base_size = 18))
tiff("FIG4.TIF", width = 4000, height = 2000, res=300)
ggarrange(p1,p2, 
          labels = c("A.", "B."),
          ncol = 2)
dev.off()




# Figure 5 ----------------------------------------------------------------

feast_output <- read.csv('timeseries_FEAST.csv')

feast_op_micro <- feast_output %>% filter(Community == 'Microbiome')

feast_op_micro$Compare <- factor(feast_op_micro$Compare, levels = c('Comp1', 'Other1',
                                                                    'Comp2', 'Other2',
                                                                    'Comp3', 'Other3',
                                                                    'Comp4', 'Other4',
                                                                    'Comp5', 'Other5',
                                                                    'Comp6', 'Other6'
                                                                    
))

p <- ggplot(feast_op_micro, aes(x=Compare, y=Contribution, fill = Compare)) + 
  geom_boxplot() + scale_fill_manual(values = c ('#2EABA8', '#87B4B3',
                                                 '#254FB6', '#8794B4',
                                                 '#5A25B6', '#927DB6',
                                                 '#A412B8', '#AF7DB6',
                                                 '#B81292', '#B972A9',
                                                 '#BC2357', "#B9728A"
  )) + theme(legend.position = "none", 
             axis.text.x = element_text(angle = 90)) +
  xlab('Collection Timepoint Comparisons') + ylab('Proportion of Contribution') +
  scale_x_discrete(labels=c("Comp1" = "1 to 2", "Other1" = "1 to 3-7",
                            "Comp2" = "2 to 3", "Other2" = "2 to 1 & 4-7",
                            "Comp3" = "3 to 4", "Other3" = "3 to 1-2 & 5-7",
                            "Comp4" = "4 to 5", "Other4" = "4 to 1-3 & 6-7",
                            "Comp5" = "5 to 6", "Other5" = "5 to 1-4 & 7",
                            "Comp6" = "6 to 7", "Other6" = "6 to 1-5"))
p


feast_op_myco <- feast_output %>% filter(Community == 'Mycobiome')

feast_op_myco$Compare <- factor(feast_op_myco$Compare, levels = c('Comp1', 'Other1',
                                                                  'Comp2', 'Other2',
                                                                  'Comp3', 'Other3',
                                                                  'Comp4', 'Other4',
                                                                  'Comp5', 'Other5',
                                                                  'Comp6', 'Other6'
                                                                  
))

q <- ggplot(feast_op_myco, aes(x=Compare, y=Contribution, fill = Compare)) + 
  geom_boxplot() + scale_fill_manual(values = c ('#BE5821', '#D19A7C',
                                                 '#D5B01B', '#D1C07C',
                                                 '#8CD51B', '#B4D581',
                                                 '#16CF23', '#81D587',
                                                 '#16CF82', '#7BCBAA',
                                                 '#1ECCC8', "#7BCBC9"
  )) + theme(legend.position = "none", 
             axis.text.x = element_text(angle = 90)) +
  xlab('Collection Timepoint Comparisons') + ylab('Proportion of Contribution') +
  scale_x_discrete(labels=c("Comp1" = "1 to 2", "Other1" = "1 to 3-7",
                            "Comp2" = "2 to 3", "Other2" = "2 to 1 & 4-7",
                            "Comp3" = "3 to 4", "Other3" = "3 to 1-2 & 5-7",
                            "Comp4" = "4 to 5", "Other4" = "4 to 1-3 & 6-7",
                            "Comp5" = "5 to 6", "Other5" = "5 to 1-4 & 7",
                            "Comp6" = "6 to 7", "Other6" = "6 to 1-5"))
q

theme_set(theme_classic(base_size = 16))
tiff("FIG5.TIF", width = 3000, height = 4000, res=300)
ggarrange(p,q, 
          labels = c("A.", "B."),
          nrow = 2, ncol = 2)
dev.off()


# Figure 6 --------------------------------------------------------------



#build random forest model
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

#pull out the predicted values from the RF regression
pred <- m1$predicted
pred <- as.data.frame(pred)
colnames(pred) <- 'Predicted'

#using the mean square error to determine the true value
pred[pred$Predicted < 7.800, "True"] <- '7'
pred[pred$Predicted < 6.800, "True"] <- '6'
pred[pred$Predicted < 5.800, "True"] <- '5'
pred[pred$Predicted < 4.800, "True"] <- '4'
pred[pred$Predicted < 3.800, "True"] <- '3'
pred[pred$Predicted < 2.800, "True"] <- '2'
pred[pred$Predicted < 1.800, "True"] <- '1'

pred$True <- as.numeric(pred$True)

#plotting predicted vs. true
ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1]), 
         color = ) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") 
}

#plot it
d <- ggplotRegression(lm(Predicted ~ True, data = pred)) +
  geom_pointrange(aes(ymin = Predicted-0.800, ymax = Predicted+0.800), 
                  position=position_jitter(width=0.5), alpha = 0.5) + xlab ('True Collection Date') +
  ylab('Predicted Collecion Date')
d
#gives the regression line
lm_eqn <- function(df){
  m <- lm(Predicted ~ True, df);
  eq <- substitute(italic(Predicted) == a + b %.% italic(True)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

d1 <- d + geom_text(x = 2, y = 6, label = lm_eqn(pred), parse = TRUE)

theme_set(theme_classic(base_size = 18))
tiff("FIG6.TIF", width = 2500, height = 2500, res=300)
d1
dev.off()




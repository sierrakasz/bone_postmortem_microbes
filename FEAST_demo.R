library(FEAST)
library(RCurl)
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
lapply(Packages, library, character.only = TRUE)


meta <-Load_metadata(metadata_path = "metadata_example_multi.txt")
otus <- Load_CountMatrix(CountMatrix_path = "otu_example_multi.txt")

FEAST_output <- FEAST(C = otus, metadata = meta, 
                      different_sources_flag = 1,
                      dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                      outfile="demo")

PlotSourceContribution(SinkNames = rownames(FEAST_output)[c(5:8)],
                       SourceNames = colnames(FEAST_output), dir_path = "C:/Users/sierr/Documents/Lucas_grant_project/R_scripts",
                       mixing_proportions = FEAST_output, Plot_title = "Test_",Same_sources_flag = 0, N = 4)

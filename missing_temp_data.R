set.seed(1234)

library(mice)
library(Hmisc)
library(DMwR)

#load data
temp_data <- read.csv("LucasHOBOdataMerged.csv")

#check to see those NAs
md.pattern(temp_data)

#use k nearest neighbor to generate temp data
knnOutput <- knnImputation(temp_data[, !names(temp_data) %in% "medv"])  # perform knn imputation.
anyNA(knnOutput)

mean(knnOutput$Temperature)
sd(knnOutput$Temperature)


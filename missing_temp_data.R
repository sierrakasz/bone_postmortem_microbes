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


library(plyr)
library(stringr)
#calculate ADH 
knnOutput_sort <- knnOutput %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
  arrange(Date)

Trtdata <- ddply(knnOutput_sort, c("Date"), summarise,
                 MeanDailyTemp=mean(Temperature)
)

Trtdata$CumulativeTemp<-cumsum(Trtdata$MeanDailyTemp) #Use the cumsum function from base R to add together temperatures
Trtdata

mean(Trtdata$CumulativeTemp)
sd(Trtdata$CumulativeTemp)

write.csv(Trtdata, file = 'meantemp_ADH_Lucas.csv')

levels = c('8/20/2018', '10/25/2018', '1/22/2019',
           '4/23/2019', '7/26/2019', '10/26/2019',
           '1/24/2020')

Trtdata_collection <- Trtdata %>% filter(Date == as.Date('2018-08-20'))

b%>% filter(Date = '2018-10-25') %>% 
  filter(Date = '2019-01-22') %>% filter(Date = '2019-04-23') %>% filter(Date = '2019-07-26') %>% 
  filter(Date = '2019-10-26') %>% filter(Date = '2020-01-24')

class(Trtdata$Date)

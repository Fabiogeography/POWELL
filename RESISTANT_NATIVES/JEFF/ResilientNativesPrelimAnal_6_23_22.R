#Explore graphs of Native abundance for Resilient Natives project
#What is widespread enough to be "common"? What would our cutoff be?
#Look on an ecoregion-by-ecoregion basis

library (dplyr)
library (tidyr)
library(ggplot2)
#library(sf)
library(sp)

options(scipen = 999)

setwd ("~/Documents/PristinePlotsProject")

#Read in datafile
Data <- read.csv ("spcisFullEco_withTRY_24May2022.csv")
Divdata <- read.csv ("FullDatabase_diversity.csv")
#PlotData <-  read.csv("SPCIS_plots.csv")
#Spdata <- read.csv ("SPCIS_plant_taxa.csv")
#both = left_join(Data, PlotData, by=c("Plot", "Year"))
BothDiv = left_join(Data, Divdata, by = c("Plot", "Year"))

glimpse (BothDiv)

#Desert ecoregion only
DesertNat <- filter (BothDiv, Ecoregion == "NORTH AMERICAN DESERTS", NativeStatus == "N")
n_distinct (DesertNat$SpCode) #3865 spp
n_distinct (DesertNat$Plot) #27,488 plots

#Count observations
Num <- count(DesertNat, SpCode)
Num <- Num [order(-Num$n), ]

#Visualize distribution
ggplot (Num, aes (x= n)) +
  geom_histogram(binwidth = 100) 

#How many in bins?
#3865 spp total

n50 <- Num %>% filter (n >50) #597
n100 <- Num %>% filter (n >100) #306
n200 <- Num %>% filter (n >200) #181
n300 <- Num %>% filter (n >300) #113
n400 <- Num %>% filter (n > 400) #90

#Visualize relationship between native cover ~ plot invasion
counts <- Num %>% filter (n >400)
RandomSpp <- DesertNat %>% filter (SpCode %in% sample (counts$SpCode, 12))
RandomSpp %>% filter (RelCov_I >2) %>%
  ggplot (aes(RelCov_I, RelCov)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")

#Eastern Forests ecoregion only
EForNat <- filter (BothDiv, Ecoregion == "EASTERN TEMPERATE FORESTS", NativeStatus == "N")
n_distinct (EForNat$SpCode) #5577 spp
n_distinct (EForNat$Plot) #30591 plots

#Count observations
Num <- count(EForNat, SpCode)
Num <- Num [order(-Num$n), ]

#Visualize distribution
ggplot (Num, aes (x= n)) +
  geom_histogram(binwidth = 100) 

#How many in bins?
#5577 spp total

n50 <- Num %>% filter (n >50) #1927
n100 <- Num %>% filter (n >100) #1403
n200 <- Num %>% filter (n >200) #904
n300 <- Num %>% filter (n >300) #674
n400 <- Num %>% filter (n > 400)#525

#Visualize relationship between native cover ~ plot invasion
counts <- Num %>% filter (n >400)
RandomSpp <- EForNat %>% filter (SpCode %in% sample (counts$SpCode, 12))
RandomSpp %>% filter (RelCov_I >2) %>%
  ggplot (aes(RelCov_I, RelCov)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")

#Great Plains ecoregion only
GPNat <- filter (BothDiv, Ecoregion == "GREAT PLAINS", NativeStatus == "N")
n_distinct (GPNat$SpCode) #2908 spp
n_distinct (GPNat$Plot) #5599 plots

#Count observations
Num <- count(GPNat, SpCode)
Num <- Num [order(-Num$n), ]

#Visualize distribution
ggplot (Num, aes (x= n)) +
  geom_histogram(binwidth = 100) 

#How many in bins?
#2908 spp total

n50 <- Num %>% filter (n >50) #323
n100 <- Num %>% filter (n >100) #104
n200 <- Num %>% filter (n >200) #61
n300 <- Num %>% filter (n >300) #38
n400 <- Num %>% filter (n > 400) #24

#Visualize relationship between native cover ~ plot invasion
counts <- Num %>% filter (n >400)
RandomSpp <- GPNat %>% filter (SpCode %in% sample (counts$SpCode, 12))
RandomSpp %>% filter (RelCov_I >2) %>%
  ggplot (aes(RelCov_I, RelCov)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")

#Make table of number of common natives in each ecoregion
EcoRegList <- unique (Natives$Ecoregion)
Natives <- filter (BothDiv, NativeStatus == "N")
Table <- c()
for (i in 1:length(unique(Natives$Ecoregion))){
  NextEcoregion <- filter (Natives, Ecoregion == EcoRegList[i])
  SppNum <- n_distinct (NextEcoregion$SpCode)
  Plots <- n_distinct (NextEcoregion$Plot)
  Num <- count(NextEcoregion, SpCode)
  n100 <- nrow(Num %>% filter (n >100))
  n200 <- nrow(Num %>% filter (n >200))
  n400 <- nrow(Num %>% filter (n > 400))
  Table <- rbind (Table, c(EcoRegList[i], Plots, SppNum, n100, n200, n400))
  print (paste (EcoRegList[i], Plots, SppNum, n100, n200, n400))
}
colnames (Table) <- c("Ecoregion", "PlotNum", "Spp in each ecoregion", "Spp w >100 obs", "Spp w >200 obs", "Spp w >400 obs")
write.table (Table, "TablePerEcoregion.csv")
EForNat <- filter (BothDiv, Ecoregion == "EASTERN TEMPERATE FORESTS", NativeStatus == "N")
n_distinct (EForNat$SpCode) #5577 spp
n_distinct (EForNat$Plot) #30591 plots

#Make list of all species that reach at least 400 observations across the whole dataset
Natives <- filter (BothDiv, NativeStatus == "N")
Num <- count(EForNat, SpCode)
n400 <- nrow(Num %>% filter (n > 400)) #451 species
counts <- Num %>% filter (n >400)
RandomSpp <- Natives %>% filter (SpCode %in% sample (counts$SpCode, 12))
RandomSpp %>% filter (RelCov_I >2) %>%
  ggplot (aes(RelCov_I, RelCov)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")
write.csv (counts, "SpeciesList_7_4_22.csv")

#Make table of number of species that are widespread and reach at least >10% RelCover
EcoRegList <- unique (Natives$Ecoregion)
Natives <- filter (BothDiv, NativeStatus == "N")
#Get names of natives with at least 10% RelCov
Common <- filter (EForNat, RelCov >10)
NatCommon <- filter (EForNat, SpCode %in% Common$SpCode)
Table <- c()
for (i in 1:length(unique(NatCommon$Ecoregion))){
  NextEcoregion <- filter (NatCommon, Ecoregion == EcoRegList[i])
  SppNum <- n_distinct (NextEcoregion$SpCode)
  Plots <- n_distinct (NextEcoregion$Plot)
  Num <- count(NextEcoregion, SpCode)
  n25 <- nrow(Num %>% filter (n >25))
  n100 <- nrow(Num %>% filter (n >100))
  n400 <- nrow(Num %>% filter (n > 400))
  Table <- rbind (Table, c(EcoRegList[i], Plots, SppNum, n25, n100, n400))
  print (paste (EcoRegList[i], Plots, SppNum, n25, n100, n400))
}

#Eastern Forests ecoregion only - abundance > 10%
#Get names of natives with at least 10% RelCov
EForNat <- filter (BothDiv, Ecoregion == "EASTERN TEMPERATE FORESTS", NativeStatus == "N")
Common <- filter (EForNat, RelCov>10)
n_distinct (Common$SpCode) #5577 spp
#NatCommon <- filter (EForNat, SpCode %in% Common$SpCode)
EForNat <- filter (BothDiv, Ecoregion == "EASTERN TEMPERATE FORESTS", NativeStatus == "N", SpCode %in% Common$SpCode)
n_distinct (EForNat$SpCode) #2170 spp
n_distinct (EForNat$Plot) #30585 plots

#Count observations
Num <- count(EForNat, SpCode)
Num <- Num [order(-Num$n), ]

#Visualize relationship between native cover ~ plot invasion for abundant species
counts <- Num %>% filter (n >400)
AbundEF <- filter (EForNat, SpCode %in% counts$SpCode)
n_distinct (AbundEF$SpCode) #5577 spp
RandomSpp <- EForNat %>% filter (SpCode %in% sample (counts$SpCode, 12))
RandomSpp %>% filter (RelCov_I >2) %>%
  ggplot (aes(RelCov_I, RelCov)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")


###############################################################################
##### EXPLORE SPECIES ABUNDANCES ########################
##### Written by: Regan Early ##########################################
##### Written on: 6th June 2022 #########################################
##### Modified on: 4th July 2022 ##########################################
########################################################################

library(dplyr)
library(ggplot2)
library(purrr) ## map_dfr
library(tidyverse) ## group_split

wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS/"
wd.out <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/RESILIENT_NATIVES/"

##### Read in species data and join environment data for each plot #####
dat <- read.csv(paste0(wd.dat, "FULLDatabase_05272022.csv")) ## 1874220 rows
dat$SpCode <- as.factor(dat$SpCode) ## easier to glance at data this way
dat <- dat[!is.na(dat$SpCode),]
length(unique(dat$SpCode)) ## 14162 species

## plots file contains the invader metrics for each plot. 
## This includes plots with > 10% relative cover of unknown species - might want to remove later. 
plots <- read.csv(paste0(wd.dat, "FULLDatabase_05272022_plotsenv4Jul2022.csv"))
plots <- plots[,!names(plots) %in% c("Dataset", "Original.Plot","Site","Original.Site","FuzzedCoord","Zone","PlotArea.m2","SamplingMethod","Resampled","Long","Lat")] ## Remove names that are also in dat

dat <- merge(plots, dat, by.x=c("Plot", "Year"), by.y=c("Plot", "Year")) ## merge based on plot and year. Plots that don't have coordinates are not retained. 1770529 rows

##### Add invader metrics for each plot ######
dat.I <- read.csv(paste0(wd.dat, "FullDatabase_diversity.csv"))
dat.I <- dat.I[,!names(dat.I) %in% c("Dataset","Site","PlotArea.m2","Long","Lat")] ## Remove names that are also in dat. except Plot and Year, the join columns

table(dat$Plot %in% dat.I$Plot) ## 143280 Plot IDs that are not in FullDatabase_diversity.csv, approx. 10% of all plots (1627249 that are)

dat <- left_join(dat, dat.I, by = c("Plot", "Year")) ## Function used by Eve, rather than merge. Weird that gained 56 rows. 

##### Identify commonest natives in the remaining plots #####
frq <- table(dat$SpCode)

length(cmn <- names(frq[frq > 400])) ## 914 species found in > 400 plots of the surviving dataset
length(cmn <- cmn[cmn %in% dat$SpCode[dat$NativeStatus=="N"]]) ## 801 of those species are native

## Calculate the mean and max coverage of the co-ocurring non-native species
length(dat[is.na(dat$NA_L1NAME),]) ## 73 points fall outside any ecoregion
dat <- dat[!is.na(dat$NA_L1NAME),] 

dat.smry <- dat[dat$SpCode %in% cmn,] %>% 
  group_by(SpCode, NA_L1NAME) %>%
  summarise(
    num_plots = n(),
    mean_cov = mean(PctCov_100),
    max_cov = max(PctCov_100)
  )

## Pick out the native species with > 400 records in **any given ecoregion**, not overall.
dim(dat.smry.400 <- dat.smry[dat.smry$num_plots>400,]) ## 753 species * ecoregion combinations where the species is present in > 400 plots
length(cmn.nms <- unique(dat.smry.400$SpCode)) ## 675 species are identified as common in at least one ecoregion

dat.cmn <- dat[dat$SpCode %in% cmn.nms,] ## Note that some of these native species are recorded as NI in some of the plots. 

## Exclude plots whose RelCov_I was <2% or that are uninvaded, to avoid the "trace" cover artifact.
## Don't do this now as it removes presences that are needed for the presence absence data frame
# dat <- dat[dat$RelCov_I>=2  & !is.na(dat$RelCov_I),] ## Lose over 2/3rds of plots. This seems high - checking with Eve.

## Save the file that contains the cover of the common natives in each plot
write.csv(dat.cmn, paste0(wd.out, "FULLDatabase_05272022_commonCov_env_IAS.csv"), row.names=F)

##### Build dataframes that include presences and absences as well as abundances for each common species #####
## For each species pull out the plots it is found and some plots in the level 4 ecoregion where it is absent 
## The number of absences in each L4 ecoregion matches the number of presences, i.e. prevalence in each L4 ecoregion is 0.5
## Alternative is to select absences from within a geographic or environmental buffer, but compare_env_ecoregions.R Found that L4 ecoregions explains the vast majority of variation in aridity, npp, human modification and min temp.
##  Level IV ecoregions are intended for large geographic extents (i.e. states, multiple counties, or river basin: https://bit.ly/3yDxR9s
## There are 967 ecoregions in the conterminous U.S.

set.seed(0) ## So that the absence plots sampled will be the same each time. Remove if you want to make multiple different samples. 

# sp <- "ACNE2"
for (sp in cmn.nms[543:length(cmn.nms)]) {

  ## Make the presence dataframe
  pres <- dat.cmn[dat.cmn$SpCode==sp,]
  pres$occ <- 1
  
  ## Get number of plots in each L4 ecoregion where the species is found
  l4.count <- as.data.frame(table(dat.cmn[dat.cmn$SpCode==sp,"US_L4NAME"])) 
  colnames(l4.count) <- c("US_L4NAME", "num_pres_plots")
  
  ## Select the same number of absences from each L4 ecoregion. This keeps the prevalence at 50% within each ecoregion.
  abs <- dat[dat$SpCode!=sp,] ## Plots where the species is absent
  abs <- left_join(abs, l4.count, by="US_L4NAME") ## Join the number of plots in each L4 ecoregion where the sepcies is found
  abs <- abs[!is.na(abs$num_pres_plots),]
  abs$occ <- 0
  
  abs <- abs %>% 
    group_split(US_L4NAME) %>% ## Returns a list of tibbles. Each tibble contains the rows of .tbl for the associated group and all the columns, including the grouping variables.
    map_dfr(~{ ## Apply a function to each element
      sample_n(.x, num_pres_plots[1], replace = F) ## sample absence plots from each L4 ecoregion in proportion to the presence plots in that L4 ecoregion.
      ## Sampling without replacement should be fine as there should be more plots where a given species is absent than present
    })
      
  ## Checking the numbers are correct
  print(paste(sp, sum(l4.count$num_pres_plots) == nrow(abs))) ## Ensure the presence plots in total matches the number absence plots selected
  
  abs$num_pres_plots <- NULL ## Remove the column indicating how many absences to select
  
  assign(paste0(sp, "_occ"), rbind(pres, abs))
  
  write.csv(get(paste0(sp, "_occ")), paste0(wd.out, sp, "_occ.csv"), row.names=F)

}

## Failed at SAVE4. Probably need to set replace=T

##### Explore with plotting #####
#Eastern Forests ecoregion only
EForNat <- filter (dat.cmn.I, NA_L1NAME == "NORTH AMERICAN DESERTS")
RandomSpp <- EForNat %>% filter (SpCode %in% sample (cmn.nms, 12))
RandomSpp %>%
  ggplot (aes(RelCov_I, PctCov_100)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + facet_wrap(~SpCode, scales = "free_y")




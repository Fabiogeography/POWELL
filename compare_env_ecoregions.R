#################################################################################
##### Visualise env vars  between and within L4 ecoregions #####
##### Written by Regan Early
##### Written on 17th June 2022
##### Modified on 17th June 2022
#################################################################################

# library(sp)
# library(rgdal)
# library(raster)
# library(maptools)
# library(rgeos)
# library(maps)

# wd.env <- "E:/NON_PROJECT/NCEAS2/ENV_DATA/ENV_LAYERS"
wd.dat <- "D:/NON_PROJECT/WORKSHOPS/POWELL/DATA/"

dat <-read.csv(paste0(wd.dat, "SPCIS_plots_env_17thJune2022.csv"))

boxplot(dat$ar ~ dat$US_L4NAME) ## 747 ecoregions
boxplot(dat$npp ~ dat$US_L4NAME)
boxplot(dat$hm ~ dat$US_L4NAME)
boxplot(dat$tmin ~ dat$US_L4NAME)

summary(lm.ar <- lm(dat$ar ~ dat$US_L4NAME)) ## L4 explains almost all the variation (96%)
summary(lm.npp <- lm(dat$npp ~ dat$US_L4NAME)) ## L4 explains slightly less variation (89%)
summary(lm.hm <- lm(dat$hm ~ dat$US_L4NAME)) ## L4 explains less of the variation (73%)
summary(lm.tmin <- lm(dat$tmin ~ dat$US_L4NAME)) ## )

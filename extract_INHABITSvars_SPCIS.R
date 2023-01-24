#################################################################################
##### Extract INHABITS environmental data at plot locations #####
##### Written by Regan Early
##### INHABITS info: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0263056
#################################################################################

.libPaths("C:/SOFTWARE/R/R-4.1.2/library")
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(maps)
library(ff)
library(prioritizr) ## fast extraction from raster
library(dplyr)

wd.env <- "F:/NON_PROJECT/WORKSHOPS/POWELL/ENVIRONMENT/INHABIT/VARIABLES/"
wd.spcis <- "F:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS/"

crs.W84 <- sp::CRS(SRS_string="EPSG:4326") ## Mercator projection string that gives planar coordinates. https://epsg.io/3857

##### Prepare INHABIT environmental data #####
## Get file names and import
rastlist <- list.files(path = wd.env, pattern='.tif', all.files=TRUE, full.names=FALSE)
setwd(wd.env)
allrasters <- lapply(rastlist, raster)
names(allrasters) <- tools::file_path_sans_ext(basename(rastlist))

## Harmonise resolution and extent. Which rasters have different crs to the others?
for (i in allrasters) {
  print(paste0(proj4string(i), "   ", names(i)))
}
## gHM, PPT_ETo_Jun_Aug_1981_2018, TMAXmean_Jun_Aug_1981_2018

ref.rast <- allrasters$Americas_N_LCM_Cat100 ## Raster with reference crs and extent

# allrasters$gHM <- projectRaster(allrasters$gHM, ref.rast)
# if(res(allrasters$gHM) == res(ref.rast)){"match"}; extent(allrasters$gHM); extent(ref.rast)
# writeRaster(allrasters$gHM, paste0(wd.env,"gHM_proj.tif"))
allrasters$gHM <- raster(paste0(wd.env,"gHM_proj.tif"))

# allrasters$TMAXmean_Jun_Aug_1981_2018 <- projectRaster(allrasters$TMAXmean_Jun_Aug_1981_2018, ref.rast)
# writeRaster(allrasters$TMAXmean_Jun_Aug_1981_2018, paste0(wd.env,"TMAXmean_Jun_Aug_1981_2018_proj.tif"))
allrasters$TMAXmean_Jun_Aug_1981_2018 <- raster(paste0(wd.env,"TMAXmean_Jun_Aug_1981_2018_proj.tif"))

# allrasters$PPT_ETo_Jun_Aug_1981_2018 <- projectRaster(allrasters$PPT_ETo_Jun_Aug_1981_2018, ref.rast) 
# if(res(allrasters$PPT_ETo_Jun_Aug_1981_2018) == res(ref.rast)){"match"}; extent(allrasters$PPT_ETo_Jun_Aug_1981_2018); extent(ref.rast)
# writeRaster(allrasters$PPT_ETo_Jun_Aug_1981_2018, paste0(wd.env,"PPT_ETo_Jun_Aug_1981_2018_proj.tif"))
allrasters$PPT_ETo_Jun_Aug_1981_2018 <- raster(paste0(wd.env,"PPT_ETo_Jun_Aug_1981_2018_proj.tif"))

env <- stack(allrasters)

## Ecoregions - not from IHNABIT dataset
ecoR <- readOGR("F:/NON_PROJECT/NCEAS2/ENV_DATA/ENV_LAYERS", "ecoR_l48")
ecoR <- spTransform(ecoR, crs(ref.rast))

##### Extract environmental data #####

## Plots with diversity metrics
div <- read.csv(paste0(wd.spcis, "FULLDatabase_10272022_diversity_PctCov_10272022.csv"))
div <- div[!is.na(div$Long),]
coordinates(div) <- ~ Long + Lat
crs(div) <- crs.W84
div <- spTransform(div,crs(ref.rast))

div.env <- base::as.data.frame(rbind(prioritizr::fast_extract(env,div))) ## INHABITS data
colnames(div.env) <- rastlist 

div.ecoR <- as.data.frame(over(div, ecoR)) ## Ecoregion data 

div.env <- as.data.frame(cbind(div, div.env, div.ecoR)) ## Combine

write.csv(div.env, paste0(wd.spcis,"FULLDatabase_10272022_diversity_PctCov_10272022_INHABIT.csv"), row.names=F) 

## Plots with diversity metrics, with pct cover constrained to 100"
div100 <- read.csv(paste0(wd.spcis, "FULLDatabase_10272022_diversity_PctCov_100_10272022.csv"))
div100 <- div100[!is.na(div100$Long),]
coordinates(div100) <- ~ Long + Lat
crs(div100) <- crs.W84
div100 <- spTransform(div100,crs(ref.rast))

div100.env <- base::as.data.frame(rbind(prioritizr::fast_extract(env,div100))) ## INHABITS data
colnames(div100.env) <- rastlist 

div100.ecoR <- as.data.frame(over(div100, ecoR)) ## Ecoregion data 

div100.env <- as.data.frame(cbind(div100, div100.env, div100.ecoR)) ## Combine

write.csv(div100.env, paste0(wd.spcis,"FULLDatabase_10272022_diversity_PctCov_100_10272022_INHABIT.csv"), row.names=F) 

## Data for each individual native and non-native species
plot <- read.csv(paste0(wd.spcis, "FULLDatabase_10272022.csv"))
plot <- div100[!is.na(plot$Long),]
coordinates(plot) <- ~ Long + Lat
crs(plot) <- crs.W84
plot <- spTransform(plot,crs(ref.rast))

plot.env <- base::as.data.frame(rbind(prioritizr::fast_extract(env,plot))) ## INHABITS data
colnames(plot.env) <- rastlist 

plot.ecoR <- as.data.frame(over(plot, ecoR)) ## Ecoregion data 

plot.env <- as.data.frame(cbind(plot, plot.env, plot.ecoR)) ## Combine

write.csv(plot.env, paste0(wd.spcis,"FULLDatabase_10272022_INHABIT.csv"), row.names=F) 


##### NEON domains - requested by Ines #####
# dat <- read.csv(paste0(out.wd, "/env_plots18Mar21_nlcd.csv"))
dom <- read.csv("F:/NON_PROJECT/NCEAS2/LOCATION_DATA/NEON_Field_Site_Metadata_20210226_0.csv")
dat$field_site_id <- substr(dat$plot, 1, 4)

dat <- merge(dat, dom[,c("field_site_id", "field_domain_id")],
             by.x="field_site_id", by.y="field_site_id",
             all.x=T)

dat$field_site_id <- NULL


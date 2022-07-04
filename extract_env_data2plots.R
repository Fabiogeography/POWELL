#################################################################################
##### Extract environmental data at plot locations #####
##### Written by Regan Early
##### Written on 20th May 2020
##### Modified on 4th July 2022
#################################################################################

library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(maps)
library(dplyr)

wd.env <- "E:/NON_PROJECT/NCEAS2/ENV_DATA/ENV_LAYERS"
wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS/"

proj.wgs <-
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

crs.merc <- CRS("+init=epsg:3857") ## Mercator projection string that gives planar coordinates. https://epsg.io/3857

##### Geographic plot locations #####
dat <- read.csv(paste0(wd.dat, "FULLDatabase_05272022.csv"))
dat <- dat[,c("Dataset","Plot","Year","Original.Plot","Site","Original.Site", "Long","Lat","FuzzedCoord","Zone","PlotArea.m2","SamplingMethod","Resampled")] ## The field names from Lais' original plots file
dat <- unique(dat) ## 87810 unique plots (87808 in the SPCIS file for publication)

# datx <- read.csv(paste0(wd.dat, "SPCIS_plots_env_17thJune2022.CSV")) ## from Eve, uploaded to Drive 17th June 2022
dat <- dat[!is.na(dat$Long), ] ##  remove plots with no coordinates. 84308 remaining.
coordinates(dat) <- ~ Long + Lat ## convert to shapefile
proj4string(dat) <- proj.wgs ## project

##### Load in environmental data #####
## Rasters
ras <- c("npp", "ar", "hm", "swd", "pt", "pwd") ## see google drive doc for definitions
for (r in ras) {
  assign(paste0(r), raster(paste0(wd.env, "/", r, "_l48.tif")))
}

# nlcd <- raster(paste0(wd.env, "/nlcd_l48.tif"))
# nlcd <- projectRaster(nlcd, to=npp) ## nlcd is in a different projection. Project once (very lengthy) and save
# writeRaster(nlcd, paste0(wd.env, "/nlcd_l48_wgs84.tif"))
nlcd <- raster(paste0(wd.env, "/nlcd_l48_wgs84.tif")) ## Projecting the raster didn't work in R. Use version projected in Arcmap instead

## Raster stacks
sta <- c("clim_wc", "prism800m") ## worldclim and prism data. Will be superceded by Terra data.
for (s in sta) {
  assign(paste0(s), stack(paste0(wd.env, "/", s, "_l48.grd")))
}

## Shapefiles
shp <- c("ecoR", "hard")
for (s in shp) {
  assign(paste0(s), readOGR(wd.env, paste0(s, "_l48")))
}

##### Extract environmental data #####
##### Rasters #####  
ras <- c("npp", "ar", "hm", "swd", "pt", "pwd", "nlcd")
for(r in ras) {
  a <- as.data.frame(extract(get(r), dat))
  colnames(a) <- r
  
  ### Find nearest env data value to points that don't intersect env data
  # Not running as failing, lengthy and not important right now
  # ## Identify points with no env data attached
  # datX <- dat[is.na(a[,1]) == T, ]
  # datX$ID <- row.names(datX)
  # datX <- spTransform(datX, crs.merc)
  # 
  # ## Draw a buffer around the points with no environmental data (distance in m)
  # # 1km buffer adds environmental data to ~200 cells but leaves ~900 points still without data, but these lie outside the study area
  # datX.buf <- buffer(datX, width = 1000, dissolve = F)
  # datX.buf$ID <- row.names(datX)
  # 
  # ## Get the numbers of the raster cells that fall in the buffers
  # if(r == "nlcd") {
  #   s.merc <- raster(paste0(wd.env, "/nlcd_l48_merc.tif"))
  # }  else { 
  #   s.merc <- projectRaster(get(r), crs=crs.merc)
  # }
  # aX.cell <- cellFromPolygon(s.merc, datX.buf)
  # names(aX.cell) <- row.names(datX)
  # 
  # ## Extract the values of the raster cells that fall in the buffers
  # aX.val <- extract(s.merc, datX.buf)
  # names(aX.val) <- row.names(datX)
  # 
  # for (i in names(aX.cell)) {
  #   ## Get the XYs of the raster cells within the buffer and convert to a data frame
  #   if (!is.null(aX.cell[[i]])) {
  #     aX.cellXY <- as.data.frame(xyFromCell(s.merc, aX.cell[[i]]))
  #     df <- na.omit(cbind(aX.val[[i]], aX.cellXY))
  #     
  #     if(nrow(df)>0) {
  #       ## Make the data frame spatial
  #       coordinates(df) <- ~ x + y ## projection not important here
  #       proj4string(df) <- crs.merc
  #       
  #       ## Calculate the distance from each point to the raster cells that fall in the buffer
  #       df$dist <-
  #         sapply(1:nrow(df), function(x) {
  #           gDistance(datX[i, ], df[x, ]) ## Returns minimum distance.
  #         }) 
  #       vals <- as.data.frame(df[which.min(df$dist), ]) ## Find the row that represents the raster cell that is closest to the point
  #       vals[, c("x", "y", "dist")] <- NULL
  #     }
  #   }
  #   ## Make a dummy vals variable if the previous code could not make it
  #   if(exists("vals")==F) { vals <- rep(NA, ncol(a)) }
  #   
  #   a[i, ] <- vals ## Add in the environmental data from the nearest raster cell to the data frame
  #   print(dim(vals))
  #   rm(vals)
  # }
  # 
  dat <- spCbind(dat, a)
}

##### Stacks #####
for (s in sta) {
  a <- as.data.frame(extract(get(s), dat))
  
  ### Find nearest env data value to points that don't intersect env data
  # Not running as failing, lengthy and not important right now
  # ## Identify points with no env data attached
  # datX <- dat[is.na(a[,1]) == T, ]
  # datX$ID <- row.names(datX)
  # datX <- spTransform(datX, crs.merc)
  # 
  # ## Draw a buffer around the points with no environmental data (distance in m)
  # # 1km buffer adds environmental data to ~200 cells but leaves ~900 points still without data, but these lie outside the study area
  # datX.buf <- buffer(datX, width = 1000, dissolve = F)
  # datX.buf$ID <- row.names(datX)
  # 
  # ## Get the numbers of the raster cells that fall in the buffers
  # s.merc <- projectRaster(get(s), crs=crs.merc)
  # aX.cell <- cellFromPolygon(s.merc, datX.buf)
  # names(aX.cell) <- row.names(datX)
  # 
  # ## Extract the values of the raster cells that fall in the buffers
  # aX.val <- extract(s.merc, datX.buf)
  # names(aX.val) <- row.names(datX)
  # 
  # # i <- names(aX.cell)[1]
  # for (i in names(aX.cell)) {
  #   ## Get the XYs of the raster cells within the buffer and convert to a data frame
  #   if (!is.null(aX.cell[[i]])) {
  #     aX.cellXY <- as.data.frame(xyFromCell(s.merc, aX.cell[[i]]))
  #     df <- na.omit(cbind(aX.val[[i]], aX.cellXY))
  #     
  #     if(nrow(df)>0) {
  #       ## Make the data frame spatial
  #       coordinates(df) <- ~ x + y 
  #       proj4string(df) <- crs.merc
  #       
  #       ## Calculate the distance from each point to the raster cells that fall in the buffer
  #       df$dist <-
  #         sapply(1:nrow(df), function(x) {
  #           gDistance(datX[i, ], df[x, ]) ## Returns minimum distance.
  #         }) 
  #       vals <- as.data.frame(df[which.min(df$dist), ]) ## Find the row that represents the raster cell that is closest to the point
  #       vals[, c("x", "y", "dist")] <- NULL
  #     }
  #   }
  #   ## Make a dummy vals variable if the previous code could not make it
  #   if(exists("vals")==F) { vals <- rep(NA, ncol(a)) }
  #   
  #   a[i, ] <- vals ## Add in the environmental data from the nearest raster cell to the data frame
  #   rm(vals)
  # }
  # 
  dat <- spCbind(dat, a)
}

##### Shapefiles #####
for (s in shp) {
  a <- as.data.frame(over(dat, get(s)))
  
  ### Find nearest env data value to points that don't intersect env data
  # Not running as failing, length and not important right now
  # ## Identify points with no env data attached
  # datX <- dat[is.na(a[,1]) == T, ] ## 682 points with no ecoR polygon
  # datX$ID <- row.names(datX)
  # 
  # ##  Project points into a planar coordinate system (here UTM zone 32)
  # datX <- spTransform(datX, crs.merc)
  # 
  # ## Draw a dissolved buffer around the points with no environmental data (distance in m)
  # # 1km buffer adds environmental data to ~200 cells but leaves ~900 points still without data, but these lie outside the study area
  # datX.buf <- buffer(datX, width = 1000, dissolve = T) 
  # 
  # ##  Project buffer into a planar coordinate system (here UTM zone 32)
  # datX.buf <- spTransform(datX.buf, crs.merc)
  # 
  # ## Clip out the polygons that intersect with teh buffers around points without data
  # s.merc <- spTransform(get(s), crs.merc)
  # aX.poly <- intersect(s.merc, datX.buf)
  # 
  # for (i in datX$ID) {
  #   g <- gDistance(datX[i,], aX.poly, byid=T)
  #   vals <- as.data.frame(aX.poly[which.min(g),])
  #   a[i, ] <- vals ## Add the environmental data from the nearest polygon in to the data frame
  #   rm(vals)
  # }
  # 
  dat <- spCbind(dat, a)
}

##### Disturbance layer - can't be reprojected so project env plots instead #####
### Lengthy and vdist info not useful, so skip
# vdist.wd <- "E:\\GIS_DATA\\USA\\DISTURBANCE\\US_VegDIST2014_10312018\\Grid"
# vdist <- raster(paste0(vdist.wd,"\\us_vdist2014"), RAT=F)
# 
# dat.p <- spTransform(dat, proj4string(vdist))
# a <- as.data.frame(extract(vdist, dat.p))
# colnames(a) <- "vdist"
# ### 453 cells with no values, but these seem to be outside the study region, and > 1000m from the nearest vdist cell, so wouldn't be assigned an env value even if buffer process used. 
# dat.safe <- dat
# dat <- spCbind(dat, a)

##### NEON domains #####
dom <- read.csv("E:/NON_PROJECT/NCEAS2/LOCATION_DATA/NEON_Field_Site_Metadata_20210226_0.csv")
# dat.safe <- dat ## in case the merge fails!
dat$field_site_id <- substr(dat$Plot, 1, 4)

dat <- merge(dat, dom[,c("field_site_id", "field_domain_id")],
              by.x="field_site_id", by.y="field_site_id",
              all.x=T)

dat$field_site_id <- NULL

##### write data #####
write.csv(dat, paste0(wd.dat, "FULLDatabase_05272022_plotsenv4Jul2022.csv"), row.names = F)

writeOGR(dat, dsn="E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS", layer="FULLDatabase_05272022_plotsenv4Jul2022", driver="ESRI Shapefile") ## Note that NAs are converted to 0.

##### Plot geographic locations #####
data(stateMapEnv)

map('state')
points(dat$Long, dat$Lat, pch=16, cex=0.5) # , col=dat$colrs


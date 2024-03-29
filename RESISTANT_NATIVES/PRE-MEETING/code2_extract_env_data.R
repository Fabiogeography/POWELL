#################################################################################
##### Extract environmental data at plot locations #####
##### Written by Regan Early
##### Written on 20th May 2020
##### Modified on 7th June 2022
#################################################################################

library(sp)
library(rgdal)
library(raster)
library(maptools)
library(rgeos)
library(maps)

out.wd <- "E:/NON_PROJECT/NCEAS2/ENV_DATA/ENV_LAYERS"
wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"

proj.wgs <-
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

crs.merc <- CRS("+init=epsg:3857") ## Mercator projection string that gives planar coordinates. https://epsg.io/3857

##### Geographic plot locations #####
dat <- read.csv(paste0(wd.dat, "IAS_common_rarest_24May2022.csv"))
# dat <- read.csv("E:\\NON_PROJECT\\NCEAS2\\LOCATION_DATA\\AllPlotLocations_18March2021.csv")
dat <- dat[!is.na(dat$Long), ] ##  in case there are plots with no long or lat data
coordinates(dat) <- ~ Long + Lat
proj4string(dat) <- proj.wgs

##### Load in environmental data #####
## Rasters
ras <- c("npp", "ar", "hm", "swd", "pt", "pwd")
for (r in ras) {
  assign(paste0(r), raster(paste0(out.wd, "/", r, "_l48.tif")))
}

# nlcd <- raster(paste0(out.wd, "/nlcd_l48.tif"))
# nlcd <- projectRaster(nlcd, to=npp) ## nlcd is in a different projection. Project once (very lengthy) and save
# writeRaster(nlcd, paste0(out.wd, "/nlcd_l48_wgs84.tif"))
nlcd <- raster(paste0(out.wd, "/nlcd_l48_wgs84.tif")) ## Projecting the raster didn't work in R. Use version projected in Arcmap instead

## Raster stacks
sta <- c("clim_wc", "prism800m")
for (s in sta) {
  assign(paste0(s), stack(paste0(out.wd, "/", s, "_l48.grd")))
}

## Shapefiles
shp <- c("ecoR", "hard")
for (s in shp) {
  assign(paste0(s), readOGR(out.wd, paste0(s, "_l48")))
}

##### Extract environmental data #####
##### Rasters #####  
ras <- c("npp", "ar", "hm", "swd", "pt", "pwd", "nlcd")
for(r in ras) {
  a <- a.safe <- as.data.frame(extract(get(r), dat))
  colnames(a) <- r
  
  # Not running as failing ad not important right now
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
  #   s.merc <- raster(paste0(out.wd, "/nlcd_l48_merc.tif"))
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
  
  # Not running as failing ad not important right now
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
  
  # Not running as failing ad not important right now
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

##### NEON domains - requested by Ines #####
# dat <- read.csv(paste0(out.wd, "/env_plots18Mar21_nlcd.csv"))
dom <- read.csv("E:/NON_PROJECT/NCEAS2/LOCATION_DATA/NEON_Field_Site_Metadata_20210226_0.csv")
dat.safe <- dat
dat$field_site_id <- substr(dat$Plot, 1, 4)

dat <- merge(dat, dom[,c("field_site_id", "field_domain_id")],
              by.x="field_site_id", by.y="field_site_id",
              all.x=T)

dat$field_site_id <- NULL

##### write data #####
write.csv(dat, paste0(wd.dat, "IAS_common_rarest_env_24May2022.csv"), row.names = F)

writeOGR(dat, dsn="E:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS", layer="IAS_common_rarest_env_24May2022", driver="ESRI Shapefile") ## Note that NAs are converted to 0.

##### Plot geographic locations #####
data(stateMapEnv)
# dat.plot <- merge(dat, colrs)
dat[!is.na(dat$SpCode.NR), "colrs"] <- "hotpink"
dat[!is.na(dat$SpCode.NC), "colrs"] <- "purple"
dat[!is.na(dat$SpCode.I), "colrs"] <- "darkgreen"

map('state')
points(dat$Long, dat$Lat, pch=16, cex=0.5, col=dat$colrs)

colrs <- data.frame(source=c("Rarest", "Commonest", "Introduced"), colrs=c("hotpink","purple","darkgreen"), stringsAsFactors=F)
legend(-126, 35, legend=colrs$source,
       col=colrs$colrs, pch=16, cex=0.7, bty="n")

map('state')
dat.NR <- dat[!is.na(dat$SpCode.NR),]
points(dat.NR$Long, dat.NR$Lat, pch=16, cex=0.5, col="hotpink")

map('state')
dat.NC <- dat[!is.na(dat$SpCode.NC),]
points(dat.NC$Long, dat.NC$Lat, pch=16, cex=0.5, col="purple")

map('state')
dat.I <- dat[!is.na(dat$SpCode.I),]
points(dat.I$Long, dat.I$Lat, pch=16, cex=0.5, col="darkgreen")



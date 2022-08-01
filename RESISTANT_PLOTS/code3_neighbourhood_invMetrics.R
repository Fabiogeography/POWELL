###############################################################################
##### CALCULATE METRICS OF INVADER COVER, RICHNESS, DIVERSITY IN THE AREA SURROUNDING EACH PLOT ########################
##### Written by: Regan Early ##########################################
##### Written on: 26th July 2022 #########################################
##### Modified on:  ##########################################
########################################################################

### start replacing raster and sp functions with terra - see jonas lembrechts twitter thread on 26th july 2022

.libPaths("D:/SOFTWARE/R-4.1.1/library")
# library(spatstat) ## nndist
library(nabor) # knn
library(magrittr)
library(rlang)
library(dplyr)
library(raster) ## buffer
library(rgdal)
library(spatialEco) ## point.in.poly
library(ncf) ## correlog
library(gstat) ## variogram

wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS/"
wd.out <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/RESISTANT_PLOTS/"

proj.wgs <-
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

proj.albers <- 
  "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs " ## Albers Equal Area Projection for North America. https://epsg.io/102008

##### Read in plot data, join environment data for each plot, project #####

## plots file contains the invader metrics for each plot. 
## This includes plots with > 10% relative cover of unknown species - might want to remove later. 
plots <- read.csv(paste0(wd.dat, "FULLDatabase_05272022_plotsenv4Jul2022.csv")) ## Made by E:\NON_PROJECT\WORKSHOPS\POWELL\R\ENV_VARS\extract_env_data2plots.R

### Load invader metrics for each plot 
dat.I <- read.csv(paste0(wd.dat, "FullDatabase_diversity_July4.csv"))
dat.I <- dat.I[,!names(dat.I) %in% c("Dataset","Site","PlotArea.m2","Long","Lat")] ## Remove names that are also in plots except Plot and Year, the join columns
table(plots$Plot %in% dat.I$Plot) ## Checking that all the plots are in the invader metrics file

### Join plot data invasion metrics
dat <- left_join(plots, dat.I, by = c("Plot", "Year")) 

### Reduce table so each plot only represented once
## Therefore, When summarising invasion in the surrounding landscape, if a plot was sampled in multiple years, I used the mean values across all sample years.
dat <- dat %>% group_by(Plot) %>%
  summarise(
    dataset = first(Dataset),
    num_yrs = n_distinct(Year),
    num_reps = n(),
    PlotArea.m2 = mean(PlotArea.m2),
    TotalPctCover_I = mean(TotalPctCover_I),
    Richness_I = mean(Richness_I),
    RelCov_I = mean(RelCov_I),
    DiversityShannon_I = mean(DiversityShannon_I),
    DiversityInvSimpson_I = mean(DiversityInvSimpson_I),
    EvennessEvar_I = mean(EvennessEvar_I),
    EvennessPielou_I = mean(EvennessPielou_I),
    Original.Long = first(Original.Long),
    Original.Lat = first(Original.Lat)
  )

## Using this code within summarise confirms that each plot has a uniquw long and lat, so using 'first' is ok.
# num_x = n_distinct(x),
# num_y = n_distinct(y)
## Comparing summary of orignal and modified dat confirms teh long and lat are maintained correctly by first()

### Make df spatial + project to equal area projection
coordinates(dat) <- ~ Original.Long + Original.Lat
proj4string(dat) <- proj.wgs
# plot(dat)
dat <- spTransform(dat, CRS = proj.albers) ## Spatial units are now in metres
# plot(dat)

### Ensure coordinates are in data as well, which can be helpful
dat$x <- coordinates(dat)[,"Original.Long"]
dat$y <- coordinates(dat)[,"Original.Lat"] 

##### Testing whether each plot has distinct coordinates#####
### Email with Lais and Dan B 26th July 2022. 
##  NPS contains different plots but with a single set of coordinates for them all.
## I also know similar situation happens for some plots of VNHP and WVNHP.
## The coordinates are assigned at the Site level
test <- nabor::knn(data=as.matrix(coordinates(dat)), k=2) 
table(test$nn.dists[,2]==0) ## 14361 (of 80782) plots are in the exact same location as another plot with a different ID (less than when Long and Lat were used - was 23120)
test2 <- dat[test$nn.dists[,2]==0,] ## which plots have duplicated coordinates?
table(test2$dataset) ## Shows which datasets contain plots with the same coordinates. Doesn't include FIA plots now

# CVS IL_CTAP     NPS    NWCA    VNHP   WVNHP 
# 696     572    2000   10893     155      45 

## All of the plots in NWCA have non-unique coordinates. Are they locationally important?
plot(test2$x[test2$dataset=="NWCA"], test2$y[test2$dataset=="NWCA"], pch=16, cex=0.5, col="red")
points(test2$x[test2$dataset=="NPS"], test2$y[test2$dataset=="NPS"], pch=16, cex=0.5, col="green")
points(test2$x[test2$dataset=="CVS"], test2$y[test2$dataset=="CVS"], pch=16, cex=0.5, col="blue")
points(test2$x[test2$dataset=="IL_CTAP"], test2$y[test2$dataset=="IL_CTAP"], pch=16, cex=0.5, col="pink")
points(test2$x[test2$dataset=="VNHP"], test2$y[test2$dataset=="VNHP"], pch=16, cex=0.5, col="purple")
points(test2$x[test2$dataset=="WVNHP"], test2$y[test2$dataset=="WVNHP"], pch=16, cex=0.5, col="orange")
points(dat$Original.Long[!(dat$Plot %in% test2$Plot)], dat$Original.Lat[!(dat$Plot %in% test2$Plot)], pch=16, cex=0.5)
## Yes, NWCA points occur in locations that other datasets don't cover. No other dataset is important. 
## Emailed Lais and Dan to see if could resolve

### For the time being, remove all plots with non-unique coordinates.
dat <- dat[!(dat$Plot %in% test2$Plot),]

##### Nearest neighbour calculations (m) #####
## Separate the invaded and uninvaded plots
plots.I <- dat[dat$RelCov_I>0, "Plot"]
plots.U <- dat[dat$RelCov_I==0, "Plot"]

## Save to try to figure out why some plots seem to share coordinates
writeOGR(plots.I, dsn="E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/RESISTANT_PLOTS", layer="plotsI", driver="ESRI Shapefile")
writeOGR(plots.U, dsn="E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/RESISTANT_PLOTS", layer="plotsU", driver="ESRI Shapefile")

## Distance from each invaded plots to nearest invaded plot
I <- as.matrix(coordinates(plots.I))
rownames(I) <- plots.I$Plot
dist.I <- nabor::knn(data=I, k=2) ## Doesn't contain any 0 distances
in500m <- dist.I$nn.dists[,2]
length(in500m <- in500m[in500m<=5000]) ## 15739 invaded plots have an invaded neighbour within 500m, 21689 within 2km, 31180 within 5km
length(dist.I$nn.dists[,2]) ## 35898 invaded plots


## Distance from each UNinvaded plots to nearest invaded plot
U <- as.matrix(coordinates(plots.U))
rownames(U) <- plots.U$Plot
dist.U <- nabor::knn(data=I, query=U, k=1) ## this way round finds the distance FROM each uninvaded plot. Doesn't contain any 0 distances
length(in500m <- dist.U$nn.dists[dist.U$nn.dists<=5000]) ## 7391 uninvaded plots have an invaded neighbour within 500m, 16286 within 2km, 23293 within 5km
length(dist.U$nn.dists) ## 30523 uninvaded plots

## Plot distance to nearest invader.
par(mfrow=c(1,2))
hist(dist.I$nn.dists[,2], #xlim=c(0,100), 
     xlim=c(0,200000),  breaks=seq(0,550000, 10000), # 10km bins
     main="Distance from each\ninvaded plot", xlab="Distance to nearest\ninvaded neighbour (m)") ## Note x-axis is truncated as there are a very few plots 550000m from the nearest invaded plot
hist(dist.U$nn.dists,
     xlim=c(0,200000), breaks=seq(0,750000, 10000), # 10km bins
     main="Distance from each\nuninvaded plot", xlab="Distance to nearest\ninvaded neighbour (m)")

##### OBTAIN INVADER METRICS FOR NEAREST NEIGHBOUR #####
## Add the ID of the nearest neighbour to invaded plots
dist.I.df <- as.data.frame(dist.I)
dist.I.df <- dist.I.df[,c("nn.idx.2","nn.dists.2")]
colnames(dist.I.df) <- c("nn.idx","nn.dist")

dist.I.df$Plot.I <- plots.I$Plot[dist.I.df$nn.idx] ## Find Plot id of the neighbour plot

plots.I.nn <- cbind(plots.I, dist.I.df)
plots.I.nn <- merge(plots.I.nn, dat, by.x="Plot", by.y="Plot")

## Add the ID of the nearest invaded neighbour to uninvaded plots. 
## Need to merge based on the Plot id of the invaded neighbour plot
dist.U.df <- as.data.frame(dist.U)
dist.U.df$Plot.I <- plots.I$Plot[dist.U.df$nn.idx] ## Find Plot id of the neighbour plot
plots.U.nn <- cbind(dist.U.df, plots.U) ## Bind distance data to the uninvaded plots data
plots.U.nn <- merge(plots.U.nn, dat, by.x="Plot.I", by.y="Plot") ## Merge with data on the invasive cover in the nearby plots

par(mfrow=c(2,2))
hist(plots.I.nn$RelCov_I, xlim=c(0,100), 
     breaks=seq(0,100, 10), 
     main="Invaded plots", xlab="Invader RelCov in nearest\ninvaded neighbour")

hist(plots.U.nn$RelCov_I, xlim=c(0,100), 
     breaks=seq(0,100, 10), 
     main="Uninvaded plots", xlab="Invader RelCov in nearest\ninvaded neighbour")

hist(plots.I.nn$Richness_I, xlim=c(0,50), 
     breaks=seq(0, 50, 5), 
     main="Invaded plots", xlab="Number of invaders in nearest\ninvaded neighbour")

hist(plots.U.nn$Richness_I, xlim=c(0,50), 
     breaks=seq(0, 50, 5), 
     main="Uninvaded plots", xlab="Number of invaders in nearest\ninvaded neighbour")

### Repeat spline correlograms from code2 and 2b now that plots with duplicate coordinates are removed
# If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram.
# If observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.

set.seed(0)
dat.trim <- dat[sample(nrow(dat), 10000),] ## Otherwise very lengthy and uses up too much memory to run

### Across all points
## RelCov_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$RelCov_I, resamp=25, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="RelCov_I") 

## TotalPctCover_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$TotalPctCover_I, resamp=5, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="TotalPctCover_I") 

## Richness_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$Richness_I, resamp=5, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="Richness_I") 

## DiversityShannon_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$DiversityShannon_I, resamp=5, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="DiversityShannon_I") 

### Across points within 10k radius of each other
## RelCov_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$RelCov_I, resamp=25, quiet=FALSE, latlon=F, xmax=1000)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="RelCov_I") 

## TotalPctCover_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$TotalPctCover_I, resamp=5, quiet=FALSE, latlon=F, xmax=1000)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="TotalPctCover_I") 

## Richness_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$Richness_I, resamp=5, quiet=FALSE, latlon=F, xmax=1000)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="Richness_I") 

## DiversityShannon_I
cgm.sp <- spline.correlog(dat.trim$Original.Long, dat.trim$Original.Lat,  dat.trim$DiversityShannon_I, resamp=5, quiet=FALSE, latlon=F, xmax=1000)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="DiversityShannon_I") 

##### Sample (empircal) variograms #####
## gamma is half the average squared distance between plots within a given distance bin
## If have outliers use robust estimator is designed to downweight observations that are unusually large or small compared to neighboring observations.
## Divides data up into 15 distance bins (default?) and calculates the variogram value gamma for the pairs of plots within each distance bin
par(mfrow=c(2,2))

vario1 <- variogram(RelCov_I ~ 1, data = dat, cutoff=2000) ## Lengthy unless data restricted or cutoff is implemented. For quicker calculations use dat[1:1000,]
plot(vario1$dist, vario1$gamma, xlab="Distance, m", ylab="gamma", main="RelCov_I") ## gamma (i.e. difference in RelCov_I between plots) peaks at 400m then drops. 

vario2 <- variogram(TotalPctCover_I ~ 1, data = dat, cutoff=2000) 
plot(vario2$dist, vario2$gamma, xlab="Distance, m", ylab="gamma", main="TotalPctCover_I") 

vario3 <- variogram(Richness_I ~ 1, data = dat, cutoff=2000) 
plot(vario3$dist, vario3$gamma, xlab="Distance, m", ylab="gamma", main="Richness_I") 

vario4 <- variogram(DiversityShannon_I ~ 1, data = dat[!is.na(dat$DiversityShannon_I),], cutoff=2000) ## Continues to increase with distance
plot(vario4$dist, vario4$gamma, xlab="Distance, m", ylab="gamma", main="DiversityShannon_I") 

##### Calculate metrics within 500m radius of focal plot #####

## Make a shapefile containing a buffer for each focal plot
radius <- 500
d.buf <- buffer(dat, width = radius, dissolve = F)
d.buf@data <- d.buf@data[,"Plot"] ## Just retain the plot id
colnames(d.buf@data) <- "plot.focal"

## Join the points that are found within the radius to the focal plot. 
pp <- point.in.poly(dat, d.buf)@data ## Just the data frame, not spatial data
## 199767 (was 588099 before dups removed) rows so on average 3 plots inside every buffer (one of which is the focal plot)

# for(i in 1:5) { ## 1:nrow(dat)
inv.metric.fn <- function(i) {
  plt <- as.character(dat@data[i,"Plot"]) ## focal plot
  
  ## Summarise some data for plots within 500m of the focal plots where invaders are present
  inv.dat <- pp[pp$plot.focal==plt & !is.na(pp$plot.focal) & pp$RelCov_I>0,] %>% 
    group_by(plot.focal) %>%
    filter(Plot != plt) %>% ## remove rows where the invaded plot is the same as the focal plot
    summarise(
      num_inv = n(), ## number of invaded plots
      PlotArea.m2_sum = sum(PlotArea.m2), ## total area of the invaded plots
      RelCov_I_sum = sum(RelCov_I), ## Could multiply by area to get aerial coverage of invaders
      TotalPctCover_I_sum = sum(TotalPctCover_I),
      Richness_I_mean = mean(Richness_I),
      DiversityShannon_I_mean = mean(DiversityShannon_I, na.rm=T), ## Diversity summaries will only be calcualted for plots with non-NA values of the metric
      DiversityInvSimpson_I_mean = mean(DiversityInvSimpson_I, na.rm=T),
      EvennessEvar_I_mean = mean(EvennessEvar_I, na.rm=T),
      EvennessPielou_I_mean = mean(EvennessPielou_I, na.rm=T)
    )
  # print(glimpse(inv.dat))
  if(nrow(inv.dat)==0) {inv.dat <- c(plt, rep(NA, ncol(inv.dat)-1))}
  inv.dat
}

inv.summary <- lapply(c(1:nrow(dat)), inv.metric.fn)
inv.summary <- as.data.frame(do.call(rbind, inv.summary))

write.csv(inv.summary, file=paste0(wd.out, "FULLDatabase_invasion_summary_,"radius,".csv"), row.names=F)


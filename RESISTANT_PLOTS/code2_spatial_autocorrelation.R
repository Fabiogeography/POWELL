###############################################################################
##### INVESTIGATE SPATIAL CORRELATION IN INVADER METRICS ###########
##### Written by: Regan Early #################################################
##### Written on: 25th July 2022 ##############################################
##### Modified on:  ##########################################
###############################################################################

.libPaths("C:/SOFTWARE/R/R-4.1.2/library")
library(dplyr) ## left_join
library(sp) ## coordinates
library(gstat) ## variogram
library(ncf) ## correlog
library(nlme) ## Variogram
library(MASS) ## glmmPQL

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

### Add invader metrics for each plot 
dat.I <- read.csv(paste0(wd.dat, "FullDatabase_diversity_July4.csv"))
dat.I <- dat.I[,!names(dat.I) %in% c("Dataset","Site","PlotArea.m2","Long","Lat")] ## Remove names that are also in plots except Plot and Year, the join columns

table(plots$Plot %in% dat.I$Plot) ## Checking that all the plots are in the invader metrics file

dat <- left_join(plots, dat.I, by = c("Plot", "Year")) ## Function used by Eve, rather than merge. Weird that gained 56 rows. 

### Make df spatial + project to equal area projection
coordinates(dat) <- ~ Long + Lat
proj4string(dat) <- proj.wgs
# plot(dat)
dat <- spTransform(dat, CRS = proj.albers) ## Spatial units are now in metres
# plot(dat)

## Optionally restrict the data to plots that contain invaders
dat <- dat[dat$RelCov_I>0,]

##### Sample (empircal) variograms #####
## gamma is half the average squared distance between plots within a given distance bin
## If have outliers use robust estimator is designed to downweight observations that are unusually large or small compared to neighboring observations.
## Divides data up into 15 distance bins (default?) and calculates the variogram value gamma for the pairs of plots within each distance bin
par(mfrow=c(2,2))

vario <- variogram(RelCov_I ~ 1, data = dat[1:1000,]) ## Lengthy. For quicker calculations use dat[1:1000,]
plot(vario$dist, vario$gamma, xlab="Distance, m", ylab="gamma", main="RelCov_I") ## gamma (i.e. difference in RelCov_I between plots) peaks at 400m then drops. 

vario <- variogram(TotalPctCover_I ~ 1, data = dat[1:1000,]) 
plot(vario$dist, vario$gamma, xlab="Distance, m", ylab="gamma", main="TotalPctCover_I") 

vario <- variogram(Richness_I ~ 1, data = dat[1:1000,]) 
plot(vario$dist, vario$gamma, xlab="Distance, m", ylab="gamma", main="Richness_I") 

vario <- variogram(DiversityShannon_I ~ 1, data = dat[!is.na(dat$DiversityShannon_I),][1:1000,]) ## Continues to increase with distance
plot(vario$dist, vario$gamma, xlab="Distance, m", ylab="gamma", main="DiversityShannon_I") 

##### Correlograms #####
## Useful webpage comparing functions: https://www.r-bloggers.com/2013/05/spatial-correlograms-in-r-a-mini-overview/#:~:text=Spatial%20correlograms%20are%20great%20to,or%20Geary's%20c)%20against%20distance.
## Checks randomness in a data set. If random, autocorrelations should be near zero for any and all time-lag separations. If non-random, then one or more of the autocorrelations will be significantly non-zero.
## The spatial (cross-)correlogram and Mantel (cross-)correlogram estimates the spatial dependence at discrete distance classes.
## The region-wide similarity forms the reference line (the zero-line); the x-intercept is thus the distance at which object are no more similar than that expected by-chance-alone across the region.
## If the data are univariate, the spatial dependence is measured by Moran's I.
## If it is multivariate, it is measured by the centred Mantel statistic. (Use correlog.nc if the non-centered multivariate correlogram is desired).
## Sig. positive score mean dataset is more spatially clustered than would be expected if underlying spatial processes were random.
## Sig. negative score means dataset is more spatially dispersed than would be expected if underlying spatial processes were random.

set.seed(0)
dat.trim <- dat[sample(nrow(dat), 1000),] ## Otherwise very lengthy

### RelCov_I
cgm <- correlog(dat.trim$Long, dat.trim$Lat,  dat.trim$RelCov_I, increment=100, resamp=5, quiet=FALSE, latlon=F)
summary(cgm$p) ## Number of significant correlations increases with resamp. Don't use a value that is too low.

plot(cgm$correlation ~ cgm$mean.of.class, type="b", pch=21,
     xlab="Distance (km)", ylab="Moran's I", main="RelCov_I", ## Distances are in km if latlon=T and gcdist is used to calculate distances
     bg=ifelse(cgm$p>0.05,"white","red"), ## fill colour is red if correlation is significant
     xlim=c(0,100000) ## Zone in on neighbour distances < 500km
)
legend("topright", c("NS","Significant"), pch=c(1,16), col=c("white","red"))
abline(h=0, lty=2)

## If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram.
## If observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.
cgm.sp <- spline.correlog(dat.trim$Long, dat.trim$Lat,  dat.trim$RelCov_I, resamp=25, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="RelCov_I") 
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="RelCov_I 500km", xlim=c(0,500000)) 

## Doesn't seem to be spatial autocorrelation at all, even when play around with increments

### TotalPctCover_I
cgm <- correlog(dat.trim$Long, dat.trim$Lat, dat.trim$TotalPctCover_I, increment=100, resamp=5, quiet=FALSE, latlon=F)
summary(cgm$p) ## Number of significant correlations increases with resamp. Don't use a value that is too low.

plot(cgm$correlation ~ cgm$mean.of.class, type="b", pch=21,
     xlab="Distance (m)", ylab="Moran's I", main="TotalPctCover_I", ## Distances are in km if latlon=T and gcdist is used to calculate distances
     bg=ifelse(cgm$p>0.05,"white","red"), ## fill colour is red if correlation is significant
     xlim=c(0,100000) ## Zone in on neighbour distances < 500km
)
legend("topright", c("NS","Significant"), pch=c(1,16), col=c("white","red"))
abline(h=0, lty=2)

## If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram.
## If observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.
cgm.sp <- spline.correlog(dat.trim$Long, dat.trim$Lat,  dat.trim$TotalPctCover_I, resamp=25, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="TotalPctCover_I") 
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="TotalPctCover_I 500km", xlim=c(0,500000)) 

### Richness_I
cgm <- correlog(dat.trim$Long, dat.trim$Lat, dat.trim$Richness_I, increment=100, resamp=5, quiet=FALSE, latlon=F)
summary(cgm$p) ## Number of significant correlations increases with resamp. Don't use a value that is too low.

plot(cgm$correlation ~ cgm$mean.of.class, type="b", pch=21,
     xlab="Distance (m)", ylab="Moran's I", main="Richness_I", ## Distances are in km if latlon=T and gcdist is used to calculate distances
     bg=ifelse(cgm$p>0.05,"white","red"), ## fill colour is red if correlation is significant
     xlim=c(0,300000) ## Zone in on neighbour distances < 500km
)
legend("topright", c("NS","Significant"), pch=c(1,16), col=c("white","red"))
abline(h=0, lty=2)

## If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram.
## If observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.
cgm.sp <- spline.correlog(dat.trim$Long, dat.trim$Lat,  dat.trim$Richness_I, resamp=25, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="Richness_I") 
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="Richness_I 500km", xlim=c(0,500000)) 

### DiversityShannon_I
dat.trim <- dat.trim[!is.na(dat.trim$DiversityShannon_I),] ## 522 records
cgm <- correlog(dat.trim$Long, dat.trim$Lat, dat.trim$DiversityShannon_I, increment=100, resamp=5, quiet=FALSE, latlon=F)
summary(cgm$p) ## Number of significant correlations increases with resamp. Don't use a value that is too low.

plot(cgm$correlation ~ cgm$mean.of.class, type="b", pch=21,
     xlab="Distance (m)", ylab="Moran's I", main="DiversityShannon_I", ## Distances are in km if latlon=T and gcdist is used to calculate distances
     bg=ifelse(cgm$p>0.05,"white","red"), ## fill colour is red if correlation is significant
     xlim=c(0,300000) ## Zone in on neighbour distances < 500km
)
legend("topright", c("NS","Significant"), pch=c(1,16), col=c("white","red"))
abline(h=0, lty=2)

## If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram.
## If observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.
cgm.sp <- spline.correlog(dat.trim$Long, dat.trim$Lat,  dat.trim$DiversityShannon_I, resamp=25, quiet=FALSE, latlon=F)
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="DiversityShannon_I") 
plot(cgm.sp, xlab="Distance (m)", ylab="Moran's I", main="DiversityShannon_I 500km", xlim=c(0,500000)) 

##### Spatial glms using penalized quasi-likelihood in MASS package #####
## See Variogram.lme for details on Variogram method used

dat.model <- dat[!is.na(dat$ar) & !is.na(dat$hm),] ## remove 0s. Can't use na.action="na.exclude" in model as then Variogram will fail
dat.model <- dat.model[!duplicated(dat.model$Long), ] ## Rewmove the few duplicate plots as otherwise can't incorporate a correlation structure in models
dat.model <- dat.model[sample(nrow(dat.model), 10000),] ## Variogram won't calculate if distance matrix too large

# m <- lme(RelCov_I ~ ar, data=dat.model, random = ~1|Dataset)

formula <- as.formula(RelCov_I ~ poly(ar,2) + poly(hm, 2))
m <- MASS::glmmPQL(formula, random=~1|Dataset,
                  data=dat.model,
                  family="gaussian")

## Plot the spatial autocorrelation in the residuals from the model
## Only the autocorrelation within each Dataset is calculated, so the x-axis is the max distance covered by the largest dataset
## By default plots the Pearson standardized residuals (raw residuals divided by the corresponding standard errors) 
## By default plots the 20 variogram values for 5% quantiles of the data
## if robust = T uses an approach thought to be robust against outliers.
plot(nlme::Variogram(m, data=dat.model, robust=T), main = "No correlation structure included in model") ## See Variogram.lme for details on Variogram method used


### exp correlation structure
m.exp <- MASS::glmmPQL(formula, random=~1|Dataset,
                   data=dat.model,
                   correlation=corExp(form=~Long+Lat), 
                   family="gaussian") ## Very lengthy

plot(nlme::Variogram(m.exp, data=dat.model, robust=T), main = "Exp correlation structure included in model") ## See Variogram.lme for details on Variogram method used




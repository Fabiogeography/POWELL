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
# library(SpatialEpi)
library(rlang)
library(dplyr)
# library(RANN)
library(raster) ## buffer
library(rgdal)
library(spatialEco) ## point.in.poly

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
    Long = first(Long),
    Lat = first(Lat)
  )

## Using this code within summarise confirms that there is only 1 long and lat per plot. 
# num_x = n_distinct(x),
# num_y = n_distinct(y)
## Comparing summary of orignal and modified dat confirms teh long and lat are maintained correctly by first()

### Make df spatial + project to equal area projection
coordinates(dat) <- ~ Long + Lat
proj4string(dat) <- proj.wgs
# plot(dat)
dat <- spTransform(dat, CRS = proj.albers) ## Spatial units are now in metres
# plot(dat)

### Ensure coordinates are in data as well, which can be helpful
dat$x <- coordinates(dat)[,"Long"]
dat$y <- coordinates(dat)[,"Lat"] 

##### Testing whether each plot has distinct coordinates#####
test <- knn(data=as.matrix(coordinates(dat)), k=2) 
table(test$nn.dists[,2]==0) ## 23120 (of 80775) plots are in the exact same location as another plot with a different ID
test2 <- dat[test$nn.dists[,2]==0,]
table(test2$dataset)
test3 <- test2[test2$dataset=="NWCA",]
### See email to Lais and Dan B 26th July 2022. Going to continue with existinc data for now.

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
dist.I <- knn(data=I, k=2)
# dist.I <- dist.I2nn.dists[,2]

## Distance from each UNinvaded plots to nearest invaded plot
U <- as.matrix(coordinates(plots.U))
rownames(U) <- plots.U$Plot
dist.U <- knn(data=I, query=U, k=1) ## this way round finds the distance FROM each uninvaded plot
# dist.U <- dist.U$nn.dists

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

##### Calculate metrics within 500m radius of focal plot #####

## Make a shapefile containing a buffer for each focal plot
d.buf <- buffer(dat, width = 500, dissolve = F)
d.buf@data <- d.buf@data[,"Plot"] ## Just retain the plot id
colnames(d.buf@data) <- "plot.focal"

## Join the points that are found within the radius to the focal plot. 
pp <- point.in.poly(dat, d.buf)@data ## Just the data frame, not spatial data
## 588099 rows so on average 7.2 plots inside every buffer (one of which is the focal plot)

# for(i in 1:5) { ## 1:nrow(dat)
inv.metric.fn <- function(i) {
  plt <- as.character(dat@data[i,"Plot"])
  
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

write.csv(inv.summary, file=paste0(wd.out, "FULLDatabase_invasion_summary_500.csv"))

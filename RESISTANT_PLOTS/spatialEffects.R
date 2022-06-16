# library(dplyr)
library(spatstat) ## nndist
library(nabor) # knn
library(magrittr, lib="D:/SOFTWARE/R-4.1.1/library")
library(SpatialEpi)
library(rlang, lib="D:/SOFTWARE/R-4.1.1/library")
library(dplyr, lib="D:/SOFTWARE/R-4.1.1/library")

wd.dat <- "D:/NON_PROJECT/WORKSHOPS/POWELL/DATA/"
wd.out <- "D:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"

dat <- read.csv(paste0(wd.dat, "spcisFullEco_withTRY_24May2022.csv"))
# dat <- dat[dat$Ecoregion=="NORTH AMERICAN DESERTS",]

## Extract the invaded and uninvaded plots only
plots.I <- na.omit(unique(dat[dat$NativeStatus=="I",c("Plot", "Long", "Lat")]))
plots.U <- na.omit(unique(dat[!(dat$Plot %in% plots.I$Plot), c("Plot", "Long", "Lat")]))

## Convert to kilometre-based grid coordinates
plots.I[,c("x","y")] <- latlong2grid(plots.I[,c("Long", "Lat")])
plots.U[,c("x","y")] <- latlong2grid(plots.U[,c("Long", "Lat")])

### Calculate nearest neighbour distances within the uninvaded and invaded plots
# dist.I <- nndist(na.omit(as.matrix(plots.I[,c("x", "y")])))
# dist.U <- nndist(na.omit(as.matrix(plots.U[,c("x", "y")])))

## Invaded plots
dist.I2 <- knn(data=na.omit(as.matrix(plots.I[,c("x", "y")])), k=2)
# dist.I2 <- dist.I2$nn.dists[,2]

##Uninvaded plots
A <- plots.U[,c("x", "y")]
B <- plots.I[,c("x", "y")]
dist.U2 <- knn(data=B, query=A, k=1) ## this way round finds the distance FROM each uninvaded plot
# dist.U2 <- dist.U2$nn.dists

## Plot distance to nearest invader.
par(mfrow=c(1,2))
hist(dist.I2$nn.dists[,2], xlim=c(0,100), 
     breaks=seq(0,600, 1), 
     main="Distance from each\ninvaded plot", xlab="Distance to nearest\ninvaded neighbour (km)")
hist(dist.U2$nn.dists,
     xlim=c(0,100), breaks=seq(0,1100, 1),
     main="Distance from each\nuninvaded plot", xlab="Distance to nearest\ninvaded neighbour (km)")

### Find invader abundance of nearest invaded plot.
## Invasion cover of each plot
RelCov.I <- dat[dat$NativeStatus=="I",] %>%
  group_by(Plot) %>%
  summarise(
    nn.inv.relcov = sum(RelCov),
    num_I = n()
  )

## Add the ID of the nearest neighbour to invaded plots
dist.I2.df <- as.data.frame(dist.I2)
dist.I2.df <- dist.I2.df[,c("nn.idx.2","nn.dists.2")]
colnames(dist.I2.df) <- c("nn.idx","nn.dist")
dist.I2.df$Plot.I <- plots.I$Plot[dist.I2.df$nn.idx] ## Find Plot id of the neighbour plot
plots.I <- cbind(plots.I, dist.I2.df)
plots.I <- merge(plots.I, RelCov.I, by.x="Plot.I", by.y="Plot")

## Add the ID of the nearest neighbour to uninvaded plots. 
## Need to merge based on the Plot id of the neighbour plot
dist.U2.df <- as.data.frame(dist.U2)
dist.U2.df$Plot.I <- plots.I$Plot[dist.U2.df$nn.idx] ## Find Plot id of the neighbour plot
plots.U <- cbind(plots.U, dist.U2.df) ## Bind distance data to the uninvaded plots data
plots.U <- merge(plots.U, RelCov.I, by.x="Plot.I", by.y="Plot") ## Merge with data on the invasive cover in the nearby plots

par(mfrow=c(1,2))
hist(plots.I$nn.inv.relcov, xlim=c(0,150), 
     breaks=seq(0,200, 10), 
     main="Invaded plots", xlab="Invader RelCov in nearest\ninvaded neighbour")

hist(plots.U$nn.inv.relcov, xlim=c(0,150), 
     breaks=seq(0,200, 10), 
     main="Uninvaded plots", xlab="Invader RelCov in nearest\ninvaded neighbour")

hist(plots.I$num_I, xlim=c(0,60), 
     breaks=seq(0, 70, 5), 
     main="Invaded plots", xlab="Number of invaders in nearest\ninvaded neighbour")

hist(plots.U$num_I, xlim=c(0,60), 
     breaks=seq(0, 70, 5), 
     main="Uninvaded plots", xlab="Number of invaders in nearest\ninvaded neighbour")




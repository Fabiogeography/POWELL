###############################################################################
##### EXPLORE SPECIES ABUNDANCES ########################
##### Written by: Regan Early ##########################################
##### Written on: 6th June 2022 #########################################
##### Modified on: 4th July 2022 ##########################################
########################################################################

library(dplyr)

wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/SPCIS/"
wd.out <- "E:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"

dat <- read.csv(paste0(wd.dat, "FULLDatabase_05272022.csv"))

length(unique(dat$SpCode)) ## 14163 species

##### Identify commonest natives #####
frq <- table(dat$SpCode)

cmn <- names(frq[frq > 400]) ## 941 species
cmn <- cmn[cmn %in% dat$SpCode[dat$NativeStatus=="N"]] ## 827 of those species are native


##### OLD #####
frq <- table(dat$SpCode)
frq <- frq[frq<1000] ## 367 species with an inordinately high number of plots - remove

par(mfrow=c(2,2))
hist(dat$PctCov[dat$PctCov<100], xlab="Pct cover < 100", main="")
hist(frq, xlab="Num. plots < 1000", main="")
boxplot(dat$PctCov ~ dat$NativeStatus)
boxplot(log(dat$PctCov) ~ dat$NativeStatus)


##### Identify rare and rarest natives #####
### Rare species - need to be sufficiently widespread but locally rare 
rare <- names(frq[frq > 150 & frq < 250]) ## 643 species are found in 150-250 plots
rare <- rare[rare %in% dat$SpCode[dat$NativeStatus=="N"]] ## 583 of those species species are native

## Calculate the mean and max coverage of the species in 150-250 plots
rare <- dat[dat$SpCode %in% rare,] %>%
  group_by(SpCode) %>%
  summarise(
    num_plots = n(),
    mean_cov = mean(PctCov),
    max_cov = max(PctCov)
  )
    
par(mfrow=c(3,2))
hist(rare$mean_cov, xlab="Mean PctCov", main="Species in 150-250 plots")
hist(rare$max_cov, xlab="Max PctCov", main="Species in 150-250 plots")

## Reduce the rare species list to those with <5% mean PctCov and <20% max PctCov
rare <- rare[rare$mean_cov<5 & rare$max_cov<20,] ## 232 species
hist(rare$mean_cov, xlab="Mean PctCov", main="Species in 150-250 plots\n<5% mean PctCov, <20% max PctCov")
hist(rare$max_cov, xlab="Max PctCov", main="Species in 150-250 plots\n<5% mean PctCov, <20% max PctCov")

## Reduce the rare species list to those with <1% mean PctCov and <5% max PctCov
rarest <- rare[rare$mean_cov<1 & rare$max_cov<5,] ## 59 species
hist(rarest$mean_cov, xlab="Mean PctCov", main="Species in 150-250 plots\n<1% mean PctCov, <5% max PctCov")
hist(rarest$max_cov, xlab="Max PctCov", main="Species in 150-250 plots\n<1% mean PctCov, <5% max PctCov")

jpeg(paste0(wd.out, "PctCov_rarest.jpg"))
par(mfrow=c(1,2))
hist(rarest$mean_cov, xlab="Mean PctCov", main="Species in 150-250 plots\n<1% mean PctCov, <5% max PctCov")
hist(rarest$max_cov, xlab="Max PctCov", main="Species in 150-250 plots\n<1% mean PctCov, <5% max PctCov")
dev.off()

## Save the data for the rare and rarest species
dat.rare <- dat[dat$SpCode %in% rare$SpCode,]
# write.csv(dat.rare, paste0(wd.out, "dat_rare.csv"), row.names=F)

dat.rarest <- dat[dat$SpCode %in% rarest$SpCode,]
# write.csv(dat.rarest, paste0(wd.out, "dat_rarest.csv"), row.names=F)

##### ID the non-natives that co-occur with the rare and rarest natives #####
## Non-natives that co-occur at least once with the rarest natives
plots.rarest <- dat$Plot[dat$Plot %in% dat.rarest$Plot] ## Plots that contain the rarest species
IAS.rarest <- unique(na.omit(dat[dat$Plot %in% plots.rarest & dat$NativeStatus=="I",c("Plot","SpCode")]))
length(unique(IAS.rarest$SpCode)) ## 583 co-occurring non-natives

## Non-natives that co-occur in at least 50 plots with the rarest natives
# test <- dat[dat$SpCode %in% rarest$SpCode ... dat$SpCode %in% IAS.rarest$SpCode]
IAS.ntv.rarest <- unique(merge(IAS.rarest, dat.rarest[,c("Plot","SpCode")], by="Plot")) ## 19469 plots that contain one of the rarest and one or more non-native species
colnames(IAS.ntv.rarest)[2:3] <- c("IAS","ntv")

IAS.ntv.rarest <- IAS.ntv.rarest %>%
  group_by(IAS, ntv) %>%
  summarise(
    num_plots = n(),
  )

hist(IAS.ntv.rarest$num_plots)
table(IAS.ntv.rarest$num_plots)
sum(table(IAS.ntv.rarest$num_plots[IAS.ntv.rarest$num_plots > 50])) ## 63 IAS:ntv species combinations where both appear in at least 50 plots together
     
IAS.ntv.rarest <- IAS.ntv.rarest[IAS.ntv.rarest$num_plots > 50,]

## Calculate the mean and max coverage of the co-ocurring non-native species
IAS.rarest <- dat[dat$SpCode %in% IAS.ntv.rarest$IAS,] %>% ## 23 co-occurring non-native species
  group_by(SpCode) %>%
  summarise(
    num_plots = n(),
    mean_cov = mean(PctCov),
    max_cov = max(PctCov)
  )

hist(IAS.rarest$num_plots, xlab="Num plots", main="Non-natives co-occurring with rarest natives")
IAS.names <- unique(IAS.rarest$SpCode[IAS.rarest$num_plots>1000]) ## 9 co-occuring non-natives in > 1000 plots

par(mfrow=c(2,2))
hist(IAS.rarest$mean_cov, xlab="Mean PctCov", main="Non-natives co-occurring with rarest\nnatives")
hist(IAS.rarest$max_cov, xlab="Max PctCov", main="Non-natives co-occurring with rarest\nnatives")

hist(IAS.rarest$mean_cov[IAS.rarest$num_plots>1000], xlab="Mean PctCov", main="Non-natives co-occurring with rarest\nnatives & in >1000 plots")
hist(IAS.rarest$max_cov[IAS.rarest$num_plots>1000], xlab="Max PctCov", main="Non-natives co-occurring with rarest\nnatives & in >1000 plots")

jpeg(paste0(wd.out, "PctCov_non-ntv.jpg"))
par(mfrow=c(1,2))
hist(IAS.rarest$mean_cov[IAS.rarest$num_plots>1000], xlab="Mean PctCov", main="Non-natives co-occurring with rarest\nnatives & in >1000 plots")
hist(IAS.rarest$max_cov[IAS.rarest$num_plots>1000], xlab="Max PctCov", main="Non-natives co-occurring with rarest\nnatives & in >1000 plots")
dev.off()

## Isolate the plots for the 9 IAS that co-occur with the rarest natives
dat.IAS.rarest <- dat[dat$SpCode %in% IAS.names,]

##### Find the commonest species that co-occur in the plots selected for the rarest species #####
## Natives that co-occur at least once with the rarest natives
common.rarest <- unique(na.omit(dat[dat$Plot %in% plots.rarest & dat$NativeStatus=="N",c("Plot","SpCode","PctCov")]))
common.rarest <- common.rarest[!(common.rarest$SpCode %in% rarest$SpCode),] ## remove the rarest species
length(unique(common.rarest$SpCode)) ## 7496 co-occurring natives
hist(common.rarest$PctCov)

### Select the 59 natives with the highest mean abundance in the plots containing the rarest species
## Calculate the mean and max coverage of the co-occurring native species
common.rarest <- common.rarest %>%
  group_by(SpCode) %>%
  summarise(
    num_plots = n(),
    mean_cov = mean(PctCov),
    max_cov = max(PctCov)
  ) %>%
arrange(desc(mean_cov), .by_group = TRUE) ## Sort by the mean abundance
        
jpeg(paste0(wd.out, "PctCov_commonest.jpeg"))
par(mfrow=c(1,2))
hist(common.rarest$mean_cov, xlab="Mean PctCov", main="Natives co-occurring with rarest\nnatives")
hist(common.rarest$max_cov, xlab="Max PctCov", main="Natives co-occurring with rarest\nnatives")
dev.off()

common.names <- common.rarest$SpCode[1:59] ## Select the 59 species with the highest mean abundance

##### Create data #####
sp.vars <- c("SpCode","PctCov","RelCov", "Duration","Growth.Habit", "Woodiness") ## Variables we want for each species
plot.vars <- c("Long", "Lat", "Dataset", "Site", "Ecoregion","PlotArea.m2","SamplingMethod") ## Variables we want for each plot

## Isolate the plots for the 9 IAS that co-occur with the rarest natives
# dat.IAS.rarest <- dat[dat$SpCode %in% IAS.names,] ## already done above
dat.IAS.ntv.rarest <- merge(dat.IAS.rarest[,c("Plot",sp.vars, plot.vars)], dat.rarest[,c("Plot",sp.vars, plot.vars)], by=c("Plot", plot.vars), all=T) ## 42193 plots. Removing all=T shows that 5207 plots that contain one of the rarest and one or more co-occurring non-native species

colnames(dat.IAS.ntv.rarest) <- gsub(".x", ".I", colnames(dat.IAS.ntv.rarest))
colnames(dat.IAS.ntv.rarest) <- gsub(".y", ".NR", colnames(dat.IAS.ntv.rarest))

## Isolate the plots for the 59 common natives that co-occur with the rarest natives
sp.vars.I <- paste0(c("SpCode","PctCov","RelCov", "Duration","Growth.Habit", "Woodiness"),".I") ## Variables we want for each species
sp.vars.NR <- paste0(c("SpCode","PctCov","RelCov", "Duration","Growth.Habit", "Woodiness"),".NR") ## Variables we want for each species
sp.vars.INR <- c(sp.vars.I, sp.vars.NR)

dat.common.rarest <- dat[dat$SpCode %in% common.names,] ## Not very many rows - the very locally common species are not very widespread!
dat.IAS.common.ntv.rarest <- merge(dat.IAS.ntv.rarest[,c("Plot",sp.vars.INR, plot.vars)], dat.common.rarest[,c("Plot",sp.vars, plot.vars)], by=c("Plot", plot.vars), all=T) ## 55781 plots. Removing all=T shows that ... plots that contain one of the rarest and one or more co-occurring non-native species
colnames(dat.IAS.common.ntv.rarest)[21:26] <- paste0(colnames(dat.IAS.common.ntv.rarest[21:26]), ".NC")

## Inspect
par(mfrow=c(3,1))
hist(dat.IAS.common.ntv.rarest$PctCov.NR)
hist(dat.IAS.common.ntv.rarest$PctCov.NC)
hist(dat.IAS.common.ntv.rarest$PctCov.I)

## Save
write.csv(dat.IAS.common.ntv.rarest, paste0(wd.out, "IAS_common_rarest_24May2022.csv"), row.names=F)



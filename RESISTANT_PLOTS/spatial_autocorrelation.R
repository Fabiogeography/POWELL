##### INVESTIGATE SPATIAL PATTERN IN INVASION #####
##### Try fitting spatial parameter in glmmmtmb #####

library(glmmTMB, lib="D:/SOFTWARE/R-4.1.1/library")
library(dplyr, lib="D:/SOFTWARE/R-4.1.1/library")
library(sp, lib="D:/SOFTWARE/R-4.1.1/library")

wd.dat <- "D:/NON_PROJECT/WORKSHOPS/POWELL/DATA/"
wd.out <- "D:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"


# Load data

# plot variogram for residuals. Once you decide the neighbourhood distance create a grouping variable to put plots within that distance into the model
library(gstat)
vario <- variogram(resid ~ 1, data = div_sub) ## at what distance do residuals look different
plot(vario$dist, vario$gamma)
# ok!! so there's definitely a plateau but then wobble again
plot(vario$dist, vario$gamma, xlim = c(0,500000), ylim=c(0.2,0.25))
# plateau at ~150 km distance...
# try correlogram - same thing but more options
library(ncf)
#Correl_GLM <- correlog(div_sub$Long, div_sub$Lat,
#                      div_sub$resid, 
#                     increment = 100,
#                    resamp=10,quiet=FALSE, latlon=TRUE)
#plot(Correl_GLM$correlation~Correl_GLM$mean.of.class, type="b", pch=21,
#    xlab="Distance", ylab="Moran's I",
#   bg=ifelse(Correl_GLM$p>0.05,0,1),
#  ylim=c(-1,1))
#legend("topright",c("Significant","NS"),pch=c(16,1))
#abline(h=0)
Spline_Correl_GLM <- spline.correlog(div_sub$Long, div_sub$Lat,
                                     div_sub$resid, resamp=10, latlon=TRUE)

plot(Spline_Correl_GLM) # very beginning and wonky at the end, but flat in the middle not bad
plot(Spline_Correl_GLM, xlim=c(0,200)) # agrees with other stuff - autcorrelation for plots within 50-100km

## Calculate autocovariate(s)
# Make a matrix of the coordinates
div_sub <- as.data.frame(div_sub)
coords<-as.matrix(cbind(div_sub$Long,div_sub$Lat))
# calculate using a couple different weighting schemes we can compare with AIC
library(spdep)
div_sub$AC_one <- autocov_dist(div_sub$Richness_N_trans, coords, nbs = 100,
                               type = "one", zero.policy=TRUE,longlat=TRUE)
div_sub$AC_inverse <- autocov_dist(div_sub$Richness_N_trans, coords,
                                   nbs = 100, type = "inverse",
                                   zero.policy=TRUE,longlat=TRUE)
div_sub$AC_inversesq <- autocov_dist(div_sub$Richness_N_trans, coords,
                                     nbs = 100, type = "inverse.squared",
                                     zero.policy=TRUE,longlat=TRUE)
summary(div_sub[,c("AC_one","AC_inverse","AC_inversesq")])

# add cov to models
mod2 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_one +
                  (1|Group), div_sub)
summary(mod2)
mod3 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_inverse +
                  (1|Group), div_sub)
summary(mod3)
mod4 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_inversesq +
                  (1|Group), div_sub)
summary(mod4)
# none have significant effects... could be neighborhood size?
AIC(mod1, mod2, mod3,mod4)
# basically the same
# try changing neighborhood size?
div_sub$AC_one_v2 <- autocov_dist(div_sub$Richness_N_trans, coords, nbs = 250,
                                  type = "one", zero.policy=TRUE,longlat=TRUE)
div_sub$AC_inverse_v2 <- autocov_dist(div_sub$Richness_N_trans, coords,
                                      nbs = 250, type = "inverse",
                                      zero.policy=TRUE,longlat=TRUE)
div_sub$AC_inversesq_v2 <- autocov_dist(div_sub$Richness_N_trans, coords,
                                        nbs = 250, type = "inverse.squared",
                                        zero.policy=TRUE,longlat=TRUE)

mod4 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_one +
                  (1|Group), div_sub)

mod5 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_one_v2 +
                  (1|Group), div_sub)
summary(mod5)
mod6 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_inverse_v2 +
                  (1|Group), div_sub)
summary(mod6)
mod7 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                  RelCov_I_trans*hm_z + 
                  RelCov_I_trans*Richness_I_zz + Dataset + AC_inversesq_v2 +
                  (1|Group), div_sub)
summary(mod7)

AIC(mod1,mod4,mod5,mod6,mod7)
# mod 5 reduces AIC by quite a lot
summary(mod5)
# test for autocorrelation!
div_sub <- div_sub %>% filter(!is.na(npp), !is.na(hm))
div_sub_sp <- st_multipoint(matrix(cbind(div_sub$Long, div_sub$Lat),,2))
box <- st_polygon(list(rbind(c(min(div_sub$Long) - 1,min(div_sub$Lat) - 1),
                             c(max(div_sub$Long) + 1,min(div_sub$Lat) - 1),
                             c(max(div_sub$Long) + 1,max(div_sub$Lat) + 1),
                             c(min(div_sub$Long) - 1,max(div_sub$Lat) + 1),
                             c(min(div_sub$Long) - 1,min(div_sub$Lat) - 1))))
# compute Voronoi polygons:
pts = st_as_sf(div_sub, coords = c("Long", "Lat"))
pols = st_collection_extract(st_voronoi(div_sub_sp, st_sfc(box)))
# match them to points:
pts$pols = pols[unlist(st_intersects(pts, pols))]
# create neighbor matrix
w <- poly2nb(pts$pols, row.names=pts$UniqueID)
summary(w)
# generate the weights matrix based on neighbor assignments
ww <-  nb2listw(w, style='B')
glimpse(pts)
pts$residuals <- residuals(mod5)
moran(pts$residuals, ww, n=length(ww$neighbours), S0=Szero(ww))
# much lower but still residual autocorrelation
set.seed(1234)
moran.mc(pts$residuals, ww, nsim=99) # significant residual autoc

# try changing neb distance again? now we know what parameter works
AC_values <- expand.grid(Neighborhood=seq(150,350, by=50),w="one")

for(i in 1:nrow(AC_values)){
  div_sub$AC <- autocov_dist(div_sub$Richness_N_trans, coords,
                             nbs=AC_values[i,"Neighborhood"],
                             type=AC_values[i,"w"], zero.policy=TRUE,
                             longlat=TRUE)
  AC_values[i,"AIC"] <- AIC(glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                                      RelCov_I_trans*hm_z + 
                                      RelCov_I_trans*Richness_I_zz + Dataset + AC, div_sub))
}
# Add the AIC from the model with no autocovariate
AC_values <- rbind(AC_values, data.frame(Neighborhood="0",w="None",AIC=AIC(mod1)))
AC_values$Delta <- AC_values$AIC-min(AC_values$AIC)
head(AC_values[order(AC_values$AIC),])
# AIC decreases with Neighborhood distance?
library(ggplot2)
AC_values %>%
  mutate(NB = as.numeric(Neighborhood)) %>%
  ggplot(aes(Neighborhood, AIC)) +
  geom_point()
# check moran's with highest nb
mod_350 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                     RelCov_I_trans*hm_z + 
                     RelCov_I_trans*Richness_I_zz + Dataset + AC, div_sub)
summary(mod_350)
pts$residuals <- residuals(mod_350)
moran(pts$residuals, ww, n=length(ww$neighbours), S0=Szero(ww))
# same val as other neighborhood distance - doesn't do much
set.seed(1234)
moran.mc(pts$residuals, ww, nsim=99) # significant residual autoc
# try increasing nb one more time and see?
AC_values <- expand.grid(Neighborhood=seq(350,500, by=50),w="one")

for(i in 1:nrow(AC_values)){
  div_sub$AC <- autocov_dist(div_sub$Richness_N_trans, coords,
                             nbs=AC_values[i,"Neighborhood"],
                             type=AC_values[i,"w"], zero.policy=TRUE,
                             longlat=TRUE)
  AC_values[i,"AIC"] <- AIC(glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                                      RelCov_I_trans*hm_z + 
                                      RelCov_I_trans*Richness_I_zz + Dataset + AC, div_sub))
}
# Add the AIC from the model with no autocovariate
AC_values <- rbind(AC_values, data.frame(Neighborhood="0",w="None",AIC=AIC(mod1)))
AC_values$Delta <- AC_values$AIC-min(AC_values$AIC)
head(AC_values[order(AC_values$AIC),])
# keeps dropping with neighborhood size...at what point does increasing neighborhood size get silly?
# p sure 500 km is way too big
mod_500 <- glmmTMB(Richness_N_trans ~ RelCov_I_trans*npp_zz + 
                     RelCov_I_trans*hm_z + 
                     RelCov_I_trans*Richness_I_zz + Dataset + AC, div_sub)
summary(mod_500)
pts$residuals <- residuals(mod_500)
moran(pts$residuals, ww, n=length(ww$neighbours), S0=Szero(ww))
# barely lower - def not the way to go
set.seed(1234)
moran.mc(pts$residuals, ww, nsim=99) # significant residual autoc
# no go!

### Now try spatial glmms
rm(list=ls())
library(MASS)
library(ggplot2)
library(dplyr)
library(nlme)
div_env <- read.csv("rawdatamodel_Jan28.csv")
glimpse(div_env)
# Clean up data
div_sub <- div_env %>% 
  mutate(Group = paste(NA_L1CODE, nlcd_veg_name1, sep="_"),
         loc = paste(Long, Lat, sep=",")) %>% 
  filter(!is.na(npp), !is.na(hm), !is.na(NA_L1NAME), 
         !is.na(nlcd_veg_name1))
div_sub <- div_sub[!duplicated(div_sub$loc),]
plot_counts <- div_sub %>%
  mutate(Group = paste(NA_L1CODE, nlcd_veg_name1, sep="_")) %>%
  group_by(Group) %>% summarise(nplot = n()) %>% filter(nplot > 20)
div_sub <- div_sub %>% filter(Group %in% plot_counts$Group)
neon <- div_sub %>% filter(Dataset=="NEON")
# fit PQL
formula <- Richness_N_trans ~ RelCov_I_trans * npp_zz + RelCov_I_trans *hm_z + RelCov_I_trans * Richness_I_zz + Dataset
model.base<- glmmPQL(formula, random=~1|Group,
                     data=div_sub,
                     family="gaussian")
summary(model.base)
plot(Variogram(model.base), main = "No Correlation Structure") # a little weird/wobbly
# exp correlation structure
PQL_model_exp <- glmmPQL(formula, random=~1|Group,
                         data=div_sub,
                         correlation=corExp(form=~Long+Lat), family="gaussian")
summary(PQL_model_exp)
plot(Variogram(PQL_model_exp)) # weirddddd
#Gaussian correlation structure
PQL_model_gaus <- glmmPQL(formula, random=~1|Group,
                          data=neon,
                          correlation=corGaus(form=~Long+Lat), family="gaussian")
summary(PQL_model_gaus)
plot(Variogram(PQL_model_gaus)) # also weirddddd

#spherical correlation structure
PQL_model_spher <- glmmPQL(formula, random=~1|Group,
                           data=neon,
                           correlation=corSpher(form=~Long+Lat), family="gaussian")
summary(PQL_model_spher)
plot(Variogram(PQL_model_spher)) # same shape as others
# compare
AIC(model.base, PQL_model_exp,PQL_model_gaus, PQL_model_spher)
# no AIC val?
# check moran?
neon_sp <- st_multipoint(matrix(cbind(neon$Long, neon$Lat),,2))
box <- st_polygon(list(rbind(c(min(neon$Long) - 1,min(neon$Lat) - 1),
                             c(max(neon$Long) + 1,min(neon$Lat) - 1),
                             c(max(neon$Long) + 1,max(neon$Lat) + 1),
                             c(min(neon$Long) - 1,max(neon$Lat) + 1),
                             c(min(neon$Long) - 1,min(neon$Lat) - 1))))
# compute Voronoi polygons:
pts = st_as_sf(neon, coords = c("Long", "Lat"))
pols = st_collection_extract(st_voronoi(neon_sp, st_sfc(box)))
# match them to points:
pts$pols = pols[unlist(st_intersects(pts, pols))]
# create neighbor matrix
w <- poly2nb(pts$pols, row.names=pts$UniqueID)
summary(w)
# generate the weights matrix based on neighbor assignments
ww <-  nb2listw(w, style='B')
glimpse(pts)
pts$residuals <- residuals(PQL_model_exp)
moran(pts$residuals, ww, n=length(ww$neighbours), S0=Szero(ww))
# reduced quite a bit!
set.seed(1234)
moran.mc(pts$residuals, ww, nsim=99) # NO LONGER SIGNIFICANT!!!
# still sig
pts$residuals_gaus <- residuals(PQL_model_gaus)
moran(pts$residuals_gaus, ww, n=length(ww$neighbours), S0=Szero(ww))
# lower than the other
set.seed(1234)
moran.mc(pts$residuals_gaus, ww, nsim=99) # still sig
# still sig
pts$residuals_spher <- residuals(PQL_model_spher)
moran(pts$residuals_spher, ww, n=length(ww$neighbours), S0=Szero(ww))
# higher than gaus
set.seed(1234)
moran.mc(pts$residuals_spher, ww, nsim=99) # still sig

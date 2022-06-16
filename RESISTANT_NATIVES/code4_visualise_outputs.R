###############################################################################
##### VISUALISE OUTPUTS FROM HMSC ########################
##### Written by: Regan Early ##########################################
##### Written on: 12th June 2022 #########################################
##### Modified on:  ##########################################
########################################################################

# library(dplyr, lib="D:/SOFTWARE/R-4.1.1/library")
# library(Hmsc, lib="D:/SOFTWARE/R-4.1.1/library")
# library(corrplot, lib="D:/SOFTWARE/R-4.1.1/library")
# set.seed(1)
library(rlang, lib="D:/SOFTWARE/R-4.1.1/library")
library(ggplot2, lib="D:/SOFTWARE/R-4.1.1/library")

# wd.dat <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/"
wd.out <- "D:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"

dat <- read.csv(paste0(wd.out,"IAS_vs_NC_pairwise_hm_npp_quad_SAFE.csv"))

### How many species pairs have a supported association? 0.95 level - increase?
d.con <- na.omit(dat[dat$support.con>=0.95, c("I", "NC", "assoc.con")])
dim(d.con) ## 4 of 211 species pairs
table(d.con$I) ## ROMU responsible for 3 and TAOF for 1

jpeg(paste0(wd.out, "assoc_constant_rarest_NC.jpg"))
boxplot(d.con$assoc.con, main="Species associations") ## Add jiggered points if possible
dev.off()
## The seven supported associations are positive

### For how many species pairs does human modification alter the association?
d.cd <- na.omit(dat[dat$support.cd>=0.95, c("I", "NC", "assoc.cd")])
dim(d.cd) ## 12 of 211 species pairs
table(d.cd$I) ## ROMU responsible for 3 and TAOF for 4, LOJA for 3

jpeg(paste0(wd.out, "assoc_hm-dep_rarest_NC.jpg"))
boxplot(d.cd$assoc.cd, main="Effect of HM on\nspecies associations") ## Add jiggered points if possible 
dev.off()
## Increasing human modification always makes the assocation more negative

### Represent difference between thm, e.g. how many negative interactions turned positive?
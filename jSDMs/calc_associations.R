###############################################################################
##### CALCULATE SPECEIS ASSOCATIONS ########################
##### Uses file obtained using find_correlations.R ######
##### Written by: Regan Early ##########################################
##### Written on: 3rd Jan 2020 #########################################
##### Modified on: 7th Jan 2020  ##########################################
########################################################################

library(MRFcov)
library(gower) ## for function cv_MRF_diag_rep
library(raster)
library(dplyr) ## unique function
library(tidyr) ## spread function
library(ggplot2) ## plots

wd <- "E:\\NON_PROJECT\\NCEAS2\\NEON\\PLANTS\\OLD\\HERB_FORMATING_DEC2019\\"
dat.all <- read.csv(paste0(wd, "plantdat1m_taxonomy_DBarnettErrorsRm.csv"), as.is=T)

### Just the two species of interest (competitors?)
dat <- dat.all[dat.all$Accepted.Symbol %in% c("POPR","CAINH2"),] # ,"BRIN2","GABO2","EUES"
dat <- dat[,c("Accepted.Symbol", "decimalLatitude", "decimalLongitude","plotID", "subplotID")] ## plot and subplotID needed to identify unique values below.

### Some random absence data - the same number of rows as in dat.mrf (but the precise number isn't really important)
dat.abs <- dat.all[sample(nrow(dat.all), 724),c("decimalLatitude", "decimalLongitude")]
# test <- unique(dat.abs)  # Only half of the records are in unique locations

##### Climate #####
clim.wd <- "E:\\GIS_DATA\\CLIMATE\\WORLDCLIM\\1950-2000\\"
clim <- raster(paste0(clim.wd,"bio1.bil"))
for(i in 2:19) {
  assign(paste0("bio",i), raster(paste0(clim.wd,"bio",i,".bil")))
  clim <- addLayer(clim, get(paste0("bio",i)))
}

dat.clim <- raster::extract(clim, dat[,c("decimalLongitude","decimalLatitude")])
dat <- cbind(dat, dat.clim)

abs.clim <- raster::extract(clim, dat.abs[,c("decimalLongitude","decimalLatitude")])
abs.clim <- abs.clim[,c(paste0("bio",c(5,6,13, 14)))]
add <- data.frame(matrix(0, nrow(abs.clim), 2, dimnames=list(c(), c("CAINH2", "POPR"))))
abs.clim <- cbind(abs.clim, add)
abs.clim <- abs.clim[,c(5,6,1,2,3,4)]

##### Format for MRFcov #####
dat$quadrat <- paste0(dat$plotID, ".", dat$subplotID)
dat <- dat[,c("Accepted.Symbol", paste0("bio",c(5,6,13, 14)), "quadrat")] ## reducing the columns

## Convert from long to wide form. http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
dat$Accepted.Symbol <- factor(dat$Accepted.Symbol)
dat <- unique(dat)

dat$occ <- 1
dat.mrf <- spread(dat, Accepted.Symbol, occ)
dat.mrf <- dat.mrf[,c(6,7,1,2,3,4)] ## Removes quadrat column as well
dat.mrf$CAINH2[is.na(dat.mrf$CAINH2)==T] <- 0
dat.mrf$POPR[is.na(dat.mrf$POPR)==T] <- 0

dat.mrf <- rbind(dat.mrf, abs.clim)

##### Run MRFcov #####
# ### The species associations only
# MRF_sponly <- MRFcov(data = dat.mrf[, c(1:2)], n_nodes=2, family='binomial') ## Fails!! Perhaps too few species?
# 
# ### The species associations and the environmental variables
MRF_env <- MRFcov(data=dat.mrf, n_nodes=2, family='binomial') ## Works
# 
# MRF_env$key_coefs ## Poa pratensis has a negative association with Carex inops ssp. heliophila, which becomes stronger with precip in the driest month (bio14). This could be interpreted as C. inops outcompeting P. pratensis in more mesic environments. Max temp is the most important variable for both species (positive effect). Min temperature has a negative effect on both species. Precip in wettest month "benefits" P. pratensis but negatively affects C. inops. Precip in driest month "benefits" C. inops, but negatively affects P. pratensis.

### Calculate model classification accuracy
mod_fits <- cv_MRF_diag_rep(data=dat.mrf, n_nodes=2,
                            n_cores=1, family='binomial', plot=F, 
                            compare_null=F, ## Fails when ask to compare against null
                            n_folds=10)
mean(mod_fits$mean_sensitivity)
mean(mod_fits$mean_specificity)
# ##### Use MRF functions to plot results. Not very exciting #####
# dev.off() ## Plots fail unless this is done.
# plotMRF_hm(MRF_mod=MRF_env,
#            node_names = c('Poa pratensis', 'Carex inops ssp. heliphila')) ##, plot_observed_vals=T, data=dat.mrf)
# 
# net <- igraph::graph.adjacency(MRF_env$graph, weighted = T, mode = "undirected")
# igraph::plot.igraph(net, layout = igraph::layout.circle,
#                     edge.width = abs(igraph::E(net)$weight),
#                     edge.color = ifelse(igraph::E(net)$weight < 0, 
#                                         'blue',
#                                         'red'))
# 
#### Bootstrap models with environmental covariates to estimate confidence intervals of effect sizes, and plot #####
MRF_env_boot <- bootstrap_MRF(data=dat.mrf, n_nodes=2, family='binomial', n_bootstraps = 10, n_cores = 1)
# MRF_env_boot
# MRF_env_boot$direct_coef_means
# MRF_env_boot$direct_coef_lower90
# MRF_env_boot$direct_coef_upper90

# ### Make plot dataframe for CAINH2
# plot.df <- data.frame(matrix(vector(), nrow=ncol(MRF_env_boot$direct_coef_means)-2, ncol=4))
# colnames(plot.df) <- c("var","meanfx","lCI90", "uCI90")
# plot.df$var <- names(MRF_env_boot$direct_coef_means["CAINH2",3:ncol(MRF_env_boot$direct_coef_means)])
# plot.df$meanfx <- MRF_env_boot$direct_coef_means["CAINH2",3:ncol(MRF_env_boot$direct_coef_means)]
# plot.df$lCI90 <- MRF_env_boot$direct_coef_lower90["CAINH2",3:ncol(MRF_env_boot$direct_coef_lower90)]
# plot.df$uCI90 <- MRF_env_boot$direct_coef_upper90["CAINH2",3:ncol(MRF_env_boot$direct_coef_upper90)]
# 
# # plot.df <- plot.df[plot.df$meanfx!=0,]
# plot.df <- plot.df[plot.df$var %in% c("POPR","bio14_POPR","bio14","bio5","bio6","bio13"),] ## For the grant proposal retain only the major effects
# 
# ggplot(plot.df, aes(x=var, y=meanfx)) +
#   geom_errorbar(aes(ymin=lCI90, ymax=uCI90), width=.1) +
#   geom_line() +
#   geom_point() +
#   ylab(" Standardised coefficient") +
#   xlab("Variable") +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle = 45, hjust = 1))

### Try to plot CAINH" and POPR together
nr <- 2*(ncol(MRF_env_boot$direct_coef_means)-1) ## Make enough rows for all variables the intercept
plot.df <- data.frame(matrix(vector(), nrow=nr, ncol=5))
colnames(plot.df) <- c("var","species","meanfx","lCI90", "uCI90")
plot.df$species[1:ncol(MRF_env_boot$direct_coef_means)-1] <- "CAINH2"
plot.df$species[ncol(MRF_env_boot$direct_coef_means):nrow(plot.df)] <- "POPR"

plot.df$var <- names(MRF_env_boot$direct_coef_means["CAINH2",2:ncol(MRF_env_boot$direct_coef_means)])

plot.df$meanfx[plot.df$species=="CAINH2"] <- MRF_env_boot$direct_coef_means["CAINH2",2:ncol(MRF_env_boot$direct_coef_means)]
plot.df$lCI90[plot.df$species=="CAINH2"] <- MRF_env_boot$direct_coef_lower90["CAINH2",2:ncol(MRF_env_boot$direct_coef_lower90)]
plot.df$uCI90[plot.df$species=="CAINH2"] <- MRF_env_boot$direct_coef_upper90["CAINH2",2:ncol(MRF_env_boot$direct_coef_upper90)]

plot.df$meanfx[plot.df$species=="POPR"] <- MRF_env_boot$direct_coef_means["POPR",2:ncol(MRF_env_boot$direct_coef_means)]
plot.df$lCI90[plot.df$species=="POPR"] <- MRF_env_boot$direct_coef_lower90["POPR",2:ncol(MRF_env_boot$direct_coef_lower90)]
plot.df$uCI90[plot.df$species=="POPR"] <- MRF_env_boot$direct_coef_upper90["POPR",2:ncol(MRF_env_boot$direct_coef_upper90)]

plot.df <- plot.df[plot.df$var %in% c("POPR","bio14_POPR","bio14","bio5","bio6"),] ## For the grant proposal retain only the major effects, and one of the BI directions

### Rename the interaction rows
plot.df$species[plot.df$var=="POPR"] <- "Both"
plot.df$species[plot.df$var=="bio14_POPR"] <- "Both" ## Rename the interaction rows
plot.df$var[plot.df$var=="POPR"] <- "Biotic Interaction"
plot.df <- plot.df[plot.df$meanfx!=0,] ## Get rid of POPR effect on itself

### Rename the variables
plot.df$var[plot.df$var=="bio14"] <- "Precip dry month"
plot.df$var[plot.df$var=="bio6"] <- "Min temp"
plot.df$var[plot.df$var=="bio5"] <- "Max temp"
plot.df$var[plot.df$var=="bio14_POPR"] <- "Precip dry month * BI"
plot.df$species[plot.df$species=="POPR"] <- "Poa"
plot.df$species[plot.df$species=="CAINH2"] <- "Carex"

### Put the levels in the right order for plotting
plot.df$species <- factor(plot.df$species, levels = c("Poa", "Carex", "Both"))

plot.df$o <- c("d","a","b","c","e","a","b","c")
plot.df <- plot.df[order(plot.df$o),]
plot.df$var <- factor(plot.df$var, levels=c("Max temp", "Min temp", "Precip dry month", "Biotic Interaction", "Precip dry month * BI"))

# jpeg("E:\\NON_PROJECT\\FUNDING\\PROJECTS_LED\\NERC_JAN2020\\FIGS\\MRFcov_fig.jpg", width=200, height=480)
p <- ggplot(plot.df, aes(x=var, y=meanfx, colour=species)) +
  geom_errorbar(aes(ymin=lCI90, ymax=uCI90), width=.1) +
  geom_line() +
  geom_point(size=3) +
  ylab("Standardised coefficient") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        text = element_text(size=20)) +
  theme(legend.justification=c(1,0), legend.position=c(0.4,0)) + # Position legend in bottom right and draw a box around it
  scale_color_manual(values = c("light grey", "dark grey", "black"), labels = levels(plot.df$species)) +
  geom_hline(yintercept=0, linetype="dashed")

p + coord_flip() 

dev.off()

## Save presence-absence dtaa associated with figure:
save(MRF_env_boot , file="E:\\NON_PROJECT\\FUNDING\\PROJECTS_LED\\NERC_JAN2020\\FIGS\\MRFcov_fig_data.R")

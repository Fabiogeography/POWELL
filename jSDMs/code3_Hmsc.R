###############################################################################
##### RUN HMSC ########################
##### Uses a modified version of computeAssociations in the file compuseAssociations.r #####
##### Written by: Regan Early ##########################################
##### Written on: 6th June 2022 #########################################
##### Modified on:  ##########################################
########################################################################

library(dplyr, lib="D:/SOFTWARE/R-4.1.1/library")
library(Hmsc, lib="D:/SOFTWARE/R-4.1.1/library")
library(corrplot, lib="D:/SOFTWARE/R-4.1.1/library")
set.seed(1)

wd.dat <- "D:/NON_PROJECT/WORKSHOPS/POWELL/DATA/"
wd.out <- "D:/NON_PROJECT/WORKSHOPS/POWELL/RESULTS/"

dat <- read.csv(paste0(wd.out, "IAS_common_rarest_env_24May2022.csv"))

### Convergence rules
nChains = 2
test.run = T
if (test.run){
  #with this option, convergence runs fast but results are not reliable - use if doing a test run
  thin = 1
  samples = 3
  transient = 5
  verbose = 0
} else {
  thin = 10
  samples = 50
  transient = 5 #*thin
  verbose = 0
}

##### Make a list of species combinations to run through later #####
combos.r <- na.omit(unique(dat[,c("SpCode.I","SpCode.NC")]))
dim(combos.r) ## 211 combinations of invader and rarest species

combos.c <- na.omit(unique(dat[,c("SpCode.I","SpCode.NC")]))
dim(combos.c) ## 85 combinations of invader and commonest species

##### Iterate over all species combinations #####
out <- as.data.frame(matrix(NA,nrow(combos.c),8))
colnames(out) <- c("I", "NC", "assoc.con", "support.con", "assoc.cd", "support.cd", "r2.I", "r2.NC")

for(i in 1:nrow(combos.c)) { ## 
  
  ##### Make a model for the environmental effects on species' abundances #####
  ### Make the data frame for the two focal species only
  I <- combos.c[i,"SpCode.I"]
  NC <- combos.c[i,"SpCode.NC"]
  
  d <- dat[dat$SpCode.I==I | dat$SpCode.NC==NC,]
  
  ## Set the cover to 0 where the focal species is absent
  d$RelCov.I[d$SpCode.I != I] <- 0
  d$RelCov.NC[d$SpCode.NC != NC | is.na(d$SpCode.NC)] <- 0
  
  ## Set the species names to NA where the focal species are absent
  d$SpCode.I[d$SpCode.I != I] <- "X"
  d$SpCode.NC[d$SpCode.NC != NC | is.na(d$SpCode.NC)] <- "X"
  
  ## Remove rows where there are NAs in the selected env vars
  d <- d[!is.na(d$npp),]
  d <- d[!is.na(d$hm),]
  
  ## Remove rows where neither of the focal species are present
  d <- d[d$SpCode.I==I | d$SpCode.NC==NC,]
  
  ### Make the model 
  Y <- as.matrix(d[,c("RelCov.I","RelCov.NC","RelCov.NC")])
  # npp = d$npp ## environmental predictor 1
  # hm = d$hm ## environmental predictor 2
  XData = data.frame(npp=d$npp,hm=d$hm)
  
  ## Inspect data to see what kind of link function might work
  par(mfrow=c(2,2))
  plot(XData$npp,Y[,1], xlab="npp", main=I)
  plot(XData$npp,Y[,2], xlab="npp", main=NC)
  plot(XData$hm,Y[,1], xlab="hm", main=I)
  plot(XData$hm,Y[,2], xlab="hm", main=NC)
  
  ##### Make and test the model object #####
  ## Including species-species associations
  studyDesign <- as.data.frame(matrix(NA,nrow(Y),3))
  studyDesign[,1] <- as.factor(d$Dataset)
  studyDesign[,2] <- as.factor(d$SamplingMethod )
  studyDesign[,3] <- as.factor(1:nrow(d))
  colnames(studyDesign) <- c("Dataset","SamplingMethod", "species")
  
  rL.Dataset <- HmscRandomLevel(units = studyDesign$Dataset)
  rL.species <- HmscRandomLevel(units = studyDesign$species)
  # rL.SamplingMethod <- HmscRandomLevel(units = studyDesign$SamplingMethod)
  
  ## Including covariate-dependent species-species associations. 
  ## The linear effect of human modification
  rL.species.hm = HmscRandomLevel(xData=data.frame(x1=rep(1,length(XData$hm)),x2=XData$hm)) ## see "working" example here: See https://github.com/hmsc-r/HMSC/issues/31

  ## To avoid fitting excessive latent factors, we constrain the minimum and maximum number of latent factors to 5 and 10, respectively.
  rL.species.hm$nfMin = 5
  rL.species.hm$nfMax = 10

  ## Set the model formula
  # XFormula <- ~x1+x2 ## Simplest - additive effect
  XFormula <- ~ poly(npp, degree = 2, raw = TRUE) + poly(hm, degree = 2, raw = TRUE) ## implement polynomial effects
  
  m <- Hmsc(Y = Y, XData = XData, XFormula = XFormula,
           studyDesign = studyDesign, ranLevels = list("species"=rL.species.hm, "Dataset"=rL.Dataset), YScale=T) ## the order of the random levels changes the order of the output of the computeAssociations and gelmanDiag functions (and others)
  
  m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, verbose = verbose, alignPost=FALSE) ## In order to implement covariate-dependent species associations alignPost=FALSE must be set. See https://github.com/hmsc-r/HMSC/issues/31. 
  
  ##### Check model assumptions are met #####
  preds <- computePredictedValues(m)

  ## I species
  preds.mean.I <- apply(preds[,1,], FUN=mean, MARGIN=1)
  nres.I <- scale(Y[,1]-preds.mean.I) ## residuals for the species combinations

  ## NC species
  preds.mean.NC <- apply(preds[,2,], FUN=mean, MARGIN=1)
  nres.NC <- scale(Y[,2]-preds.mean.NC) ## residuals for the species combinations

  # ##### Check MCMC convergence diagnostics. #####
  # mpost = convertToCodaObject(m) ## converts the posteriors into a named list of the posteriors of the various elements of the model
  # plot(mpost$Beta)
  # effectiveSize(mpost$Beta) ## effective sample sizes are very close to the theoretical value of the actual number of samples, which is 20 (10 per chain)
  # 
  # ### Potential Scale Reduction Factors
  # ## Potential scale reduction factors close to one indicates that each chain give consistent results (confirm with visual inspection of the trace plot)
  # 
  # ## PSRF of the environmental effects
  # gelman.diag(mpost$Beta, multivariate=FALSE)$psrf 
  # 
  # ## PSRFs of the random levels, including species associations
  # gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf ## The top level of this list corresponds to the number of random levels. If the first random level was used to indicate species associations, then [[1]] will report the geolman diagonal for associations. the next level corresponds to the number of chains (one list of posteriors per chain)

  ##### Evaluate model fit #####
  ## Explanatory power
  preds <- computePredictedValues(m)
  r2 <- evaluateModelFit(hM = m, predY = preds)$R2 ##

  # ## Predictive power, 2-fold cross-evaluation
  # partition <- createPartition(m, nfolds = 2)
  # preds <- computePredictedValues(m, partition = partition)
  # evaluateModelFit(hM = m, predY = preds) ## as above

  # ## Conditional predictions, i.e. if we know the abundance of some other species
  # ## conceptually similar to using non-focal species as predictors.
  # preds = computePredictedValues(m, partition=partition,
                                 # partition.sp=c(1,2), mcmcStep=10) ## full leave-one-out cross validation across the five species. Uses MCMC sampling with 10 iterations.
  # ## recommend setting mcmcStep first e.g. to 10 and then to 100 to see if the results improve
  # evaluateModelFit(hM = m, predY = preds)
  # 
  # postBeta <- getPostEstimate(m, parName = "Beta")
  # plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95)
  # 
  # ##### Investigate relationships with environmental variables #####
  # ## Look at the parameter estimates
  # mpost = convertToCodaObject(m)
  # summary(mpost$Beta)
  # 
  # ## Plot relationships
  # Gradient <- constructGradient(m, focalVariable="x1", non.focalVariables=list(x2=list(1)))
  # predY <- predict(m, Gradient=Gradient, expected=T)
  # plotGradient(m, Gradient, predY, measure="S")
  # 
  ##### Investigate species-species associations #####
  ### Constant
  OmegaCor.con <- computeAssociations.regan(m, l="constant") ## list of association matrices corresponding to each random effect in the model and in the order that they were specified. Covariances are converted to correlations.
  supportLevel <- 0.95
  assoc.con <- OmegaCor.con[[1]]$mean[1,2] ## Correlation between the two species (assuming species is the first random level added above)
  supp.con <- OmegaCor.con[[1]]$support[1,2] ## Posterior probability for the correlation being negative or positive
  
  ### Covariate-dependent
  OmegaCor.cd <- computeAssociations.regan(m, l="cov-dependent") ## list of association matrices corresponding to each random effect in the model and in the order that they were specified. Covariances are converted to correlations.
  supportLevel <- 0.95
  assoc.cd <- OmegaCor.cd[[1]]$mean[1,2] ## Correlation between the two species (assuming species is the first random level added above)
  supp.cd <- OmegaCor.cd[[1]]$support[1,2] ## Posterior probability for the correlation being negative or positive
  
  ##### Make single plots for this species combo #####
  jpeg(paste0(wd.out, I, "_", NC, ".jpg"), width = 500, height=500)
  par(mfrow=c(2,4))
  
  ## Fits assumptions of link function?
  hist(nres.I, las = 1)
  plot(preds.mean.I, nres.I, las = 1)
  abline(a=0,b=0)
  
  hist(nres.NC, las = 1)
  plot(preds.mean.NC, nres.NC, las = 1)
  abline(a=0,b=0)
  
  # ## Convergence diagnostics
  # hist(effectiveSize(mpost$Beta), main="ess(beta)")
  # hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
  # hist(gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega - species effect)")
  
  ## Species associations - constant
  toPlot <- ((OmegaCor.con[[1]]$support>supportLevel) + 
               (OmegaCor.con[[1]]$support<(1-supportLevel))>0)*OmegaCor.con[[1]]$mean
  corrplot(toPlot, method = "color",
           col = colorRampPalette(c("blue","white","red"))(200),
           title = "Constant associations", mar=c(0,0,1,0))
  
  ## Species associations - dependent on selected environmental variable
  toPlot <- ((OmegaCor.cd[[1]]$support>supportLevel) + 
               (OmegaCor.cd[[1]]$support<(1-supportLevel))>0)*OmegaCor.cd[[1]]$mean
  corrplot(toPlot, method = "color",
           col = colorRampPalette(c("blue","white","red"))(200),
           title = "Associations dependent on human modification", mar=c(0,0,1,0))

  dev.off()
  
  ##### Output #####
  # *** Add metrics of convergence properties to this once can create mpost
  o <- c(I, NC, round(assoc.con, 5), round(supp.con, 5), round(assoc.cd, 5), round(supp.cd, 5), 
         round(r2[1], 5), round(r2[2], 5))
  # names(o) <- colnames(out)
  
  out[i,] <- o
  write.csv(out, paste0(wd.out,"IAS_vs_NC_pairwise_hm_npp_quad.csv"), row.names=F)
}
write.csv(out, paste0(wd.out,"IAS_vs_NC_pairwise_hm_npp_quad.csv"), row.names=F)

## I'm running pairwise associations. Is there a better procedure that could look at multiple species?


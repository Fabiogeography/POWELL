###############################################################################
##### Prelim models ########################
##### Written by: Regan Early ##########################################
##### Written on: 7th July 2022 #########################################
##### Modified on:  ##########################################
########################################################################

library(dplyr)
library(ggplot2)
library(Hmsc)
library(corrplot)
set.seed(0)

wd.out <- "~/POWELL/"

### Convergence rules
nChains = 2
test.run = F
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

##### Load data #####
dat.cmn <- read.csv(paste0(wd.out, "FULLDatabase_05272022_commonCov_env_IAS.csv"))
splist <- unique(dat.cmn$SpCode)

##### Create output file

##### Iterate over all species combinations #####
# sp <- splist[1]
for(sp in splist[1:3]) {
  
  d <- read.csv(paste0(wd.out, sp, "_occ.csv"))
  d$RelCov_I[is.na(d$RelCov_I)] <- 0 ## If invader cover is NA, set it to be 0
  d <- na.omit(d[,c("RelCov_I","PctCov_100", "ar", "hm", "Dataset")])
  # d <- d[d$PctCov_100>0,] ## Select only plots where focal native is present. Data 0-inflated. Doesn't help.
  colnames(d)[1:2] <- c("invader", "native")
  
  ##### Visualise data #####
  ## How often does the cover of the taget native and the invasives approach 100%? This might be an issue for interpreting the results.
  ggplot(d, aes(x= (native+invader))) + 
    geom_histogram(color="black", fill="white") +
    ggtitle(sp)
  
  ## Compare cover of invasives and target native
  d %>%
    ggplot (aes(invader, native)) + geom_point() +
    geom_smooth(method='loess', formula= y~x) +
    ggtitle(sp)
  
  ## Inspect data to see what kind of link function might work
  par(mfrow=c(2,2))
  plot(d$ar,d$invader, xlab="ar", main="I")
  plot(d$ar,d$native, xlab="ar", main="N")
  plot(d$hm,d$invader, xlab="hm", main="I")
  plot(d$hm,d$native, xlab="hm", main="N")
  
  ##### Make and test the model object #####
  Y <- as.matrix(d[,c("invader","native")]) ## species data
  XData = data.frame(ar=d$ar,hm=d$hm) ## environmental data
  
  ## Including species-species associations
  studyDesign <- as.data.frame(matrix(NA,nrow(d),2))
  studyDesign[,1] <- as.factor(d$Dataset)
  studyDesign[,2] <- as.factor(1:nrow(d))
  
  colnames(studyDesign) <- c("dataset", "species")
  
  rL.dataset <- HmscRandomLevel(units = studyDesign$dataset)
  rL.species <- HmscRandomLevel(units = studyDesign$species)
  
  ## Set the model formula
  # XFormula <- ~x1+x2 ## Simplest - additive effect
  XFormula <- ~ poly(ar, degree = 2, raw = TRUE) + poly(hm, degree = 2, raw = TRUE) ## implement polynomial effects
  
  m <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
            distr=c("normal", "normal"), ## Can have different observation models for the different 'species', i.e. binomial for native occurrence and normal for invader cover
            studyDesign = studyDesign, ranLevels = list("species"=rL.species, "dataset"=rL.dataset), YScale=T) ## the order of the random levels changes the order of the output of the computeAssociations and gelmanDiag functions (and others)
  
  m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, alignPost=FALSE) ## In order to implement covariate-dependent species associations alignPost=FALSE must be set. See https://github.com/hmsc-r/HMSC/issues/31. 
  
  ##### Check model assumptions are met #####
  preds <- computePredictedValues(m)
  colnames(preds) <- colnames(Y)
  
  ## Invader cover
  preds.mean.I <- apply(preds[,"invader",], FUN=mean, MARGIN=1)
  nres.I <- scale(Y[,"invader"]-preds.mean.I) ## residuals for the species combinations
  
  ## Focal native species
  preds.mean.N <- apply(preds[,"native",], FUN=mean, MARGIN=1)
  nres.N <- scale(Y[,"native"]-preds.mean.N) ## residuals for the species combinations
  
  ## Plot
  jpeg(paste0(wd.out, "resids_", sp,".jpg"))
  par(mfrow=c(2,2))
  
  hist(nres.I, las = 1)
  plot(preds.mean.I, nres.I, las = 1)
  abline(a=0,b=0)
  
  hist(nres.N, las = 1)
  plot(preds.mean.N, nres.N, las = 1)
  abline(a=0,b=0)
  dev.off()
  
  ##### Check MCMC convergence diagnostics. #####
  mpost <- convertToCodaObject(m) ## converts the posteriors into a named list of the posteriors of the various elements of the model
  # plot(mpost$Beta) ## trace plots
  effectiveSize(mpost$Beta) ## effective sample sizes should be very close to the theoretical value of the actual number of samples, which is nChains * samples (see values set above)
  
  ### Potential scale reduction factors (PSRF) of the environmental effects. Values close to one indicates that each chain give consistent results (confirm with visual inspection of the trace plot)
  gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
  
  ##### Evaluate model fit #####
  ## Explanatory power
  preds <- computePredictedValues(m)
  (expR2 <- evaluateModelFit(hM = m, predY = preds)$R2)
  
  ## Predictive power, 2-fold cross-evaluation
  partition <- createPartition(m, nfolds = 2)
  preds <- computePredictedValues(m, partition = partition)
  evaluateModelFit(hM = m, predY = preds) ## as above
  predR2 <- evaluateModelFit(hM = m, predY = preds)$R2
  
  ## Conditional predictions, i.e. if we know the abundance/occurrence of one species, what will the other be? Conceptually similar to using non-focal species as predictors.
  preds = computePredictedValues(m, partition=partition, partition.sp=c(1,2), mcmcStep=10) ## full leave-one-out cross validation across the five species. Uses MCMC sampling with 10 iterations.
  ## recommend setting mcmcStep first e.g. to 10 and then to 100 to see if the results improve
  evaluateModelFit(hM = m, predY = preds)
  
  postBeta <- getPostEstimate(m, parName = "Beta")
  # plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95)
  
  ##### Investigate relationships with environmental variables #####
  ## Look at the parameter estimates
  summary(mpost$Beta)
  
  ## Plot relationships
  jpeg(paste0(wd.out, sp, "_jSDM.jpg"))
  par(mfrow=c(2,2))
  
  Gradient <- constructGradient(m, focalVariable="ar", non.focalVariables=list(x2=list(1)))
  predY <- predict(m, Gradient=Gradient, expected=T)
  plotGradient(m, Gradient, predY, measure="S")
  
  Gradient <- constructGradient(m, focalVariable="hm", non.focalVariables=list(x2=list(1)))
  predY <- predict(m, Gradient=Gradient, expected=T)
  plotGradient(m, Gradient, predY, measure="S")
  
  ##### Investigate species-species associations #####
  OmegaCor <- computeAssociations(m) ## Extract associations from model object. Yields a list of association matrices corresponding to each random effect in the model and in the order that they were specified. This also converts the covariances to the more convenient scale of correlation (ranging from -1 to +1)
  supportLevel <- 0.95 ## plot only those associations for which the posterior probability for being negative or positive is at least 0.95
  
  (assoc <- OmegaCor[[1]]$mean[1,2]) ## Association: Correlation between the two species (assuming species is the first random level added above)
  (supp <- OmegaCor[[1]]$support[1,2]) ## Support: Posterior probability for the correlation being negative or positive
  
  toPlot = ((OmegaCor[[1]]$support>supportLevel)
            + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
  corrplot(toPlot, method = "color",
           col = colorRampPalette(c("blue","white","red"))(200),
           title = paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))
  
  plot(d$invader, d$native, main=sp)
  dev.off()
  
  ##### Write results to a data frame #####
  out[out$FocalNative==sp, ] <- c(sp, table((d$native + d$invader)>80)[2], expR2, predR2, ## How often summed cover is > 80% and the R2 values
                                  as.data.frame(summary(mpost$Beta)$statistics)$Mean, ## environmental parameter estimates
                                  assoc, supp, ## species associations
                                  as.data.frame(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf)[,1], ## psrf, 
                                  (nChains*samples) - effectiveSize(mpost$Beta) ## difference between effective sample size and theoretical value of the actual number of samples
  )
  
}

write.csv(out, paste0(wd.out, "AAA_Hmsc_output.csv"), row.names=F)


# Try  Goldfeld-Quandt test test for homoscedasticity: https://www.r-bloggers.com/2021/11/homoscedasticity-in-regression-analysis/. Only orks for linear models?
# It compares variances of two subgroups; one set of high values and one set of low values. If the variances differ, the test rejects the null hypothesis that the variances of the errors are not constant.
# Removing central one third of observations is optional, but recommended. 
###############################################################################
##### Prelim models ########################
##### Written by: Regan Early ##########################################
##### Written on: 7th July 2022 #########################################
##### Modified on:  ##########################################
########################################################################

.libPaths("D:/SOFTWARE/R-4.1.2/library")
library(dplyr)
library(ggplot2)
library(ggpubr) ## ggarrange
library(Hmsc)
library(corrplot)
set.seed(0)

wd.out <- "E:/NON_PROJECT/WORKSHOPS/POWELL/DATA/RESILIENT_NATIVES/"

model.type <- 5

if(model.type==1) {dn <- "normal"; ext <- "/NATIVE_NORMAL"} ## 1 is natives modelled with normal distribution
if(model.type==2) {dn <- "poisson"; ext <- "/NATIVE_POISSON"} ## 2 is natives modelled with poisson distribution
if(model.type==3) {dn <- "lognormal poisson"; ext <- "/NATIVE_LOGNORMAL_POISSON"} ## 3 is natives modelled with lognormal poisson distribution
if(model.type==4) {dn <- "normal"; ext <- "/NATIVE_INVADER_ABUNDANCE_ONLY_NORMAL"} ## 4 is native abundance modelled with normal distribution
if(model.type==5) {dn <- "normal"; ext <- "/NATIVE_INVADER_ABUNDANCE_ONLY_LOGTRANSFORM_NORMAL"} ## 5 is native abundance log-transformed and modelled with normal distribution
if(model.type==6) {dn <- "..."} ## 6 is native presence absence modelled with binomial distribution

### Convergence rules
nChains = 2
test.run = F
if (test.run){
  ## With this option, convergence runs fast but results are not reliable - use if doing a test run
  thin <- 1
  samples <- 3
  transient <- 5
  verbose <- 0
} else {
  ## Options used in Hmsc vignette 1. Lengthy.
  thin <- 5
  samples <- 1000
  transient <- 500*thin
  verbose <- 500*thin
}

##### Load data #####
dat.cmn <- read.csv(paste0(wd.out, "FULLDatabase_05272022_commonCov_env_IAS.csv")) ## Does not contain absences
splist <- unique(dat.cmn$SpCode)

##### Create output file
out <- as.data.frame(matrix(data=NA, nrow=length(splist), ncol=39))
colnames(out) <- c("FocalNative", "nplots", "exceed80pc", "expR2_I", "expR2_N", "predR2_I", "predR2_N",
                   "int_I", "ar_I", "ar2_I", "hm_I", "hm2_I", "int_N", "ar_N", "ar2_N", "hm_N", "hm2_N",
                   "assoc", "supp",
                   "psrf_int_I", "psrf_ar_I", "psrf_ar2_I", "psrf_hm_I", "psrf_hm2_I", "psrf_int_N", "psrf_ar_N", "psrf_ar2_N", "psrf_hm_N", "psrf_hm2_N", ## Potential scale reduction factors (PSRF) of the environmental effects
                   "diff_int_I", "diff_ar_I", "diff_ar2_I", "diff_hm_I", "diff_hm2_I", "diff_int_N", "diff_ar_N", "diff_ar2_N", "diff_hm_N", "diff_hm2_N") ## difference between effective sample size and theoretical value of the actual number of samples

out$FocalNative <- splist

##### Iterate over all species combinations #####
# sp <- splist[1]
for(sp in splist[1:20]) {
  
  d <- read.csv(paste0(wd.out, sp, "_occ.csv")) ## The PctCov_100 column contains 0s, which were calculated by randomly sampling plots without the target native from the same L4 ecoregion as the target native.
  d$RelCov_I[is.na(d$RelCov_I)] <- 0 ## If invader cover is NA, set it to be 0
  d <- na.omit(d[,c("RelCov_I","PctCov_100", "ar", "hm", "Dataset", "NA_L1NAME")]) 
  colnames(d)[1:2] <- c("invader", "native")
  
  ## If making abundance portion of hurdle model - replace 0s with NA and possibly log-transform
  if(model.type %in% c(4,5)) {
    d$invader[d$invader==0] <- NA
    d$native[d$native==0] <- NA
  }
  
  if(model.type==5) {
    d$native <- log(d$native)
  }
  
  ##### Visualise data #####
  i.plot <- ggplot(d, aes(x=invader)) + 
    geom_histogram(color="black", fill="white")
  
  n.plot <- ggplot(d, aes(x=native)) + 
    geom_histogram(color="black", fill="white")
  
  ## How often does the cover of the target native and the invasives approach 100%? This might be an issue for interpreting the results.
  both.plot <- ggplot(d, aes(x= (native+invader))) + 
    geom_histogram(color="black", fill="white")
  
  ## Compare cover of invasives and target native
  comp.plot <- d %>%
    ggplot (aes(invader, native)) + geom_point() +
    geom_smooth(method='loess', formula= y~x)

  ## Inspect data to see what kind of link function might work
  ar.i <- ggplot(d, aes(x=ar, y=invader)) + 
    geom_point()
  
  ar.n <- ggplot(d, aes(x=ar, y=native)) + 
    geom_point()
  
  hm.i <- ggplot(d, aes(x=hm, y=invader)) + 
    geom_point()
  
  hm.n <- ggplot(d, aes(x=hm, y=native)) + 
    geom_point()

  ## Plot everything
  ggarrange(i.plot, n.plot, both.plot, comp.plot, ar.i, ar.n, hm.i, hm.n,
            labels = c("A", "B", "C","D","E","F","G","H"),
            ncol = 2, nrow = 4) %>%
    ggexport(filename = paste0(wd.out, ext, "/", sp,"_data.jpg"), width=350)

  ##### Make and test the model object #####
  Y <- as.matrix(d[,c("invader","native")]) ## species data
  XData = data.frame(ar=d$ar,hm=d$hm) ## environmental data
  
  ## Including species-species associations
  studyDesign <- as.data.frame(matrix(NA,nrow(d),3))
  studyDesign[,1] <- as.factor(d$Dataset)
  studyDesign[,2] <- as.factor(1:nrow(d))
  studyDesign[,3] <- as.factor(d$NA_L1NAME)
  
  colnames(studyDesign) <- c("dataset", "species", "l1_ecoregion")
  
  rL.dataset <- HmscRandomLevel(units = studyDesign$dataset)
  rL.species <- HmscRandomLevel(units = studyDesign$species)
  rL.l1eco <- HmscRandomLevel(units = studyDesign$l1_ecoregion)
  
  ## Set the model formula
  # XFormula <- ~x1+x2 ## Simplest - additive effect
  XFormula <- ~ poly(ar, degree = 2, raw = TRUE) + poly(hm, degree = 2, raw = TRUE) ## implement polynomial effects
  
  m <- Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
            distr=c("normal", dn), ## Can have different observation models for the different 'species', i.e. binomial for native occurrence and normal for invader cover
            studyDesign = studyDesign, ranLevels = list("species"=rL.species, "dataset"=rL.dataset, "l1_ecoregion"=rL.l1eco), YScale=T) ## the order of the random levels changes the order of the output of the computeAssociations and gelmanDiag functions (and others)
  
  m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose, alignPost=FALSE) ## In order to implement covariate-dependent species associations alignPost=FALSE must be set. See https://github.com/hmsc-r/HMSC/issues/31. 
  
  ##### Check model assumptions are met #####
  preds <- computePredictedValues(m)
  colnames(preds) <- colnames(Y)
  
  ## Invader cover
  preds.mean.I <- apply(preds[,"invader",], FUN=mean, MARGIN=1) ## The predicted value at each observed value of invader cover, averaged over all the samples and chains
  nres.I <- scale(Y[,"invader"]-preds.mean.I) ## residuals for the species combinations
  
  ## Focal native species
  preds.mean.N <- apply(preds[,"native",], FUN=mean, MARGIN=1)
  nres.N <- scale(Y[,"native"]-preds.mean.N) ## residuals for the species combinations
  
  ## Plot
  jpeg(paste0(wd.out, ext, "/", sp,"_resids.jpg"))
  par(mfrow=c(2,2))
  
  hist(nres.I, las = 1, main="Invader residuals")
  plot(preds.mean.I, nres.I, las = 1, main="Invader resids against fitted")
  abline(a=0,b=0)
  
  hist(nres.N, las = 1, main="Native residuals")
  plot(preds.mean.N, nres.N, las = 1, main="Native resids against fitted")
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
  jpeg(paste0(wd.out, ext, "/", sp, "_jSDM.jpg"))
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
  out[out$FocalNative==sp, ] <- c(sp, nrow(d), table((d$native + d$invader)>80)[2], expR2, predR2, ## How often summed cover is > 80% and the R2 values
                                  as.data.frame(summary(mpost$Beta)$statistics)$Mean, ## environmental parameter estimates
                                  assoc, supp, ## species associations
                                  as.data.frame(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf)[,1], ## psrf, 
                                  (nChains*samples) - effectiveSize(mpost$Beta) ## difference between effective sample size and theoretical value of the actual number of samples
  )
  
}

write.csv(out, paste0(wd.out, ext, "/AAA_Hmsc_output.csv"), row.names=F)


# Try  Goldfeld-Quandt test test for homoscedasticity: https://www.r-bloggers.com/2021/11/homoscedasticity-in-regression-analysis/. Only orks for linear models?
# It compares variances of two subgroups; one set of high values and one set of low values. If the variances differ, the test rejects the null hypothesis that the variances of the errors are not constant.
# Removing central one third of observations is optional, but recommended. 
# Run manually: https://www.statisticshowto.com/goldfeld-quandt-test/
# Tests the ratio of mean square residual errors for the regressions on the two subsets of data.
# This corresponds to the F-Test for equality of variances. Both the one-tailed and two-tailed tests can be used.

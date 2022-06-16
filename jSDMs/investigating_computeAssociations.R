.libPaths("D:/SOFTWARE/R-4.1.1/library")
### Try to install version from github
install.packages("devtools") # if not yet installed

library(cli, lib="D:/SOFTWARE/R-4.1.1/library")
library(sessioninfo, lib="D:/SOFTWARE/R-4.1.1/library")
library(prettyunits, lib="D:/SOFTWARE/R-4.1.1/library")
library(memoise, lib="D:/SOFTWARE/R-4.1.1/library")
library(pkgbuild, lib="D:/SOFTWARE/R-4.1.1/library")
library(magrittr, lib="D:/SOFTWARE/R-4.1.1/library")
library(rlang, lib="D:/SOFTWARE/R-4.1.1/library")
library(usethis, lib="D:/SOFTWARE/R-4.1.1/library")
library(devtools, lib="D:/SOFTWARE/R-4.1.1/library")
install_github("hmsc-r/HMSC")

library(Hmsc) ## instlaled from github on 10th June 2022
studyDesign <- data.frame(sample = as.factor(1:50),
                          plot = as.factor(sample(1:20, 50, replace = TRUE)))
xData = data.frame(x1=rep(1, length(TD$X$x1)), x2=as.numeric(TD$X$x2=="c")) ## Samples where env var x2 has a value of o become 0 and where env var x2 has a value of c become 1.
rL1 <- HmscRandomLevel(xData = xData)
rL2 <- HmscRandomLevel(units = TD$studyDesign$plot)
rL3 <- HmscRandomLevel(units = TD$studyDesign$sample)
rL1$nfMin = 5
rL1$nfMax = 10 ## Adding these does not fix the problem. 

m.tik <- Hmsc(Y = TD$Y,XData = TD$X, XFormula = ~x1+x2,
          studyDesign = studyDesign, ranLevels = list("sample" = rL1, "plot" = rL2))
m.tik <- sampleMcmc(m.tik, samples = 100, alignPost=FALSE)

OmegaCor.tik <- computeAssociations.regan(m.tik) ##  
## list of association matrices corresponding to each random level in the model and in the order that they were specified. Covariances are cnverted to correlations.

m.simple <- Hmsc(Y = TD$Y,XData = TD$X, XFormula = ~x1+x2,
              studyDesign = studyDesign, ranLevels = list("plot" = rL2, "sample" = rL3))
m.simple <- sampleMcmc(m.simple, samples = 100, alignPost=FALSE)

OmegaCor.simple <- computeAssociations(m.simple) ## Computes species associations given each of the random effects, in the order that the random effects are entered into the model 

###################################################################
mpost = convertToCodaObject(ps)
## Fails: Error in dimnames(x) <- dn : 
# length of 'dimnames' [2] not equal to array extent   

mpost <- convertToCodaObject(ps, Lambda=F, Omega=F, Psi=F,Delta=F) ## fails at omega, Psi, Delta (slightly differently).

OmegaCor <- computeAssociations.regan(m) ## Fails: `colnames<-`(`*tmp*`, value = hM$spNames) m$spnames 
## list of association matrices corresponding to each random level in the model and in the order that they were specified. Covariances are cnverted to correlations.


m.simple <- Hmsc(Y = Y, XData = XData, XFormula = XFormula,
         studyDesign = studyDesign, ranLevels = list("species"=rL.species, "Dataset"=rL.Dataset), YScale=T) ## the order of the random levels changes the order of the output of the computeAssociations and gelmanDiag functions (and others)

m.simple <- sampleMcmc(m.simple, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose, alignPost=FALSE) ## In order to implement covariate-dependent species associations alignPost=FALSE must be set. See https://github.com/hmsc-r/HMSC/issues/31. 


## Can get beta characteristics, but that's cos they don't depend on the covariate-dependent random effect
plot(mpost$Beta)
effectiveSize(mpost$Beta)
# gelman.diag(mpost$Beta, multivariate=FALSE)$psrf 


r <- set to 1

## misc
xData.regan <- data.frame(x1=rep(1,length(XData$npp)),x2=XData$hm)





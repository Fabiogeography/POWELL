# hM <- m ## object produced by sampleMcmc and this random level formulation rL.species.hm = HmscRandomLevel(xData=data.frame(x1=rep(1,length(XData$npp)),x2=XData$hm)) ## see "working" example here: See https://github.com/hmsc-r/HMSC/issues/31
# hM <- m.simple ## 
# hM <- m.tik ## 4 species

# start <- 1
# thin <- 1
# r <- 1

## Am I using the right element/s of Lambda? 
## nr comes from Hmsc: hM$nr = ncol(hM$Pi). nr IS THE NUMBER OF RANDOM LEVELS

library(abind)
computeAssociations.regan <- function (hM, l=c("constant","cov-dependent"), start = 1, thin = 1) ## l is loadings
{
  OmegaCor <- vector("list", hM$nr)
  postList <- poolMcmcChains(hM$postList, start = start, thin = thin)
  # getOmegaCor.orig <- function(a, r = r) return(cov2cor(crossprod(a$Lambda[[r]]))) ## The original is not producing the right output
  # getOmegaCor = function(a, r = r) return(cov2cor(crossprod(a$Lambda[[1]][,,r]))) ## This seems to be producing the right kind of output!!!
  getOmegaCor <- function(a, r = r) {
    if(length(dim(a$Lambda[[r]]))==3) {
      if(l=="constant"){rr <- 1} ## select the first matrix
      if(l =="cov-dependent"){rr <- 2} ## select the second matrix
      ## This assumes the first matrix is the constant loadings, and the second is the loadings conditional on the environmental predictor
      return(cov2cor(crossprod(a$Lambda[[1]][,,rr])))
    } else {
      return(cov2cor(crossprod(a$Lambda[[r]]))) 
    }
  }
  
  # a <- postList[[1]]
  # a$Lambda[[1]] ## The first element (which should correspond to the covariate-dependent random effect of species) has two matrices, which seem quite similar to each other. Maybe these are what I need?
  ## The covariate-dependent random effect (rL1) gives an element of the Lambda array with two matrices in it, each with the number of columns and rows matching the number of species. 
  ## When Dataset random level is included as well, second in the code, the Lambda array has two elements [[1]] and [[2]], and the second one seems to correspond to the Dataset random effect.
  
  for (r in seq_len(hM$nr)) { ## Runs through the random levels
    OmegaCor1 <- lapply(postList, getOmegaCor, r = r) ## posteriors. Length is the number of chains * the number of samples, e.g. 22 for 2 chains and 11 samples. In the working version of the MCMC object each of the 22 levels produces a 2*2 matrix, seemingly showing the species associations in the sample. In Tikhonov's version, this is a list of 1s. 
    
    mOmegaCor1 <- apply(abind(OmegaCor1, along = 3), c(1, 2), mean) ## A species * species matrix, which is the mean of the matrices in the posterior lambda values
    OmegaCor2 <- lapply(OmegaCor1, function(a) return(a > 0))
    support1 <- apply(abind(OmegaCor2, along = 3), c(1, 2), mean)
    
    colnames(mOmegaCor1) = hM$spNames
    rownames(mOmegaCor1) = hM$spNames
    colnames(support1) = hM$spNames
    rownames(support1) = hM$spNames
    tmp = list()
    tmp$mean = mOmegaCor1
    tmp$support = support1
    OmegaCor[[r]] = tmp
  }
  return(OmegaCor)
}

## COuld two matrices correspond to lamda(x*i.) and lambda(x*i.)T, which are multiplied to get the omega matrix of covariances? 
## No, because T means transposed and the matrices are not transposed versions of each other. 
## The dot notation (.) indicates that the variable is a vector.
## x*i is the vector of predictors

## Could the first matrix be the constant latent variable loadings, and the second one be the loadings onto the latent factor conditional on the environmental predictor?
## This would mean that the two matrices should be summed to get the latent loadings (eqn 4 in Tikhonov 2017)
## THey would then need to be multiplied by the latent factors to get the linear predictor
## The latter seems unlikely, since the code above doesn't do that. 

##### message to github 
# Hi, I'm still struggling to get the covariate-dependent functionality to work.
# The code that @gtikhonov provoided on 28th Jan 2019 makes a model, but computeAssociations fails:

# Error in dimnames(x) <- dn : 
#   length of 'dimnames' [2] not equal to array extent

# This seems to be because in teh postList object created by poolMcmcChains in computeAssociations, 
# The element of Lambda that corresponds to the covariate-dependent random level (rL1) has two species * species matrices
postList[[1]]$Lambda

[[1]]
, , 1

[,1]          [,2]          [,3]          [,4]
[1,]  0.0262100640  0.0186711924  5.804635e-02  0.1077498535
[2,] -0.0136396407 -0.0261227143 -5.455211e-03 -0.0122868940
[3,]  0.0013262262  0.0002285267  6.567342e-03 -0.0011980602
[4,] -0.0004133883 -0.0007863265 -3.854753e-05  0.0004773631

, , 2

[,1]          [,2]         [,3]          [,4]
[1,] -0.2373420188 -0.0828967318 -0.066711818 -0.1553379441
[2,] -0.0478642459  0.0118648434 -0.032209240 -0.0578322325
[3,]  0.0001250210  0.0081198502 -0.003553992  0.0140523202
[4,] -0.0002287827 -0.0001740797  0.001141877  0.0001855759


[[2]]
[,1]        [,2]        [,3]        [,4]
[1,] 0.20153787 -0.10386425 -0.08746174 -0.04134730
[2,] 0.04805722  0.01782121  0.03143461  0.02981933

# I can update computeAssociations to work around this, but 
# I can't figure out how to use these two matrices to understand how the species associations are affected by the environmental variable. 

# convertToCodaObject also fails, I think for a similar reason.



# convert to coda 
# mpost = convertToCodaObject(ps)
# ## Fails: 
# Error in dimnames(x) <- dn : 
#  length of 'dimnames' [2] not equal to array extent    
# traceback()
# 2: `colnames<-`(`*tmp*`, value = sprintf("Lambda%d[%s, factor%s]", 
#        r, spNames[rep(1:ns, nfMax[r])], as.character(rep(1:nfMax[r], 
#            each = ns))))
# 
# The same issue happens not just with Lambda but with Delta, Psi, and Omega
# 
# I'm using the package installed from github (though the same happens when I use the version installed from CRAN).
# 
# The code that @gtikhonov provoided on 28th Jan 2019 makes a model but computeAssociations won't computeAssociations
# 
# > OmegaCor <- computeAssociations(m) ## Fails: `colnames<-`(`*tmp*`, value = hM$spNames) m$spnames 
# Error in dimnames(x) <- dn : 
#   length of 'dimnames' [2] not equal to array extent
# 
# In computeAssociations I updated the line:
# getOmegaCor = function(a, r = r) return(cov2cor(crossprod(a$Lambda[[r]]))) 
# with
# getOmegaCor = function(a, r = r) return(cov2cor(crossprod(a$Lambda[[1]][,,r])))

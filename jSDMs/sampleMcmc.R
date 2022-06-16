sampleMcmc
function (hM, samples, transient = 0, thin = 1, initPar = NULL, 
          verbose, adaptNf = rep(transient, hM$nr), nChains = 1, nParallel = 1, 
          useSocket = .Platform$OS.type == "windows", dataParList = NULL, 
          updater = list(), fromPrior = FALSE, alignPost = TRUE) 
{
  if (nParallel > 1 && .Platform$OS.type == "windows" && 
      !useSocket) {
    useSocket <- TRUE
    message("only socket clusters can be used in Windows: setting useSocket = TRUE")
  }
  if (missing(verbose)) {
    if (samples * thin <= 50) 
      verbose <- 1
    else verbose <- samples * thin/50
  }
  verbose <- as.integer(verbose)
  if (fromPrior) 
    nParallel = 1
  force(adaptNf)
  if (nParallel > nChains) {
    nParallel <- nChains
    message("using ", nParallel, " cores for ", 
            nChains, " chains")
  }
  if (any(adaptNf > transient)) 
    stop("transient parameter should be no less than any element of adaptNf parameter")
  X1 = hM$XScaled
  if (hM$ncsel > 0) {
    if (is.matrix(X1)) {
      X2 = X1
      X1 = list()
      for (j in 1:hM$ns) {
        X1[[j]] = X2
      }
    }
  }
  Tr = hM$TrScaled
  Y = hM$YScaled
  distr = hM$distr
  Pi = hM$Pi
  dfPi = hM$dfPi
  C = hM$C
  nr = hM$nr
  mGamma = hM$mGamma
  iUGamma = chol2inv(chol(hM$UGamma))
  V0 = hM$V0
  f0 = hM$f0
  aSigma = hM$aSigma
  bSigma = hM$bSigma
  rhopw = hM$rhopw
  if (is.null(dataParList)) 
    dataParList = computeDataParameters(hM)
  Qg = dataParList$Qg
  iQg = dataParList$iQg
  RQg = dataParList$RQg
  detQg = dataParList$detQg
  rLPar = dataParList$rLPar
  hM$postList = vector("list", nChains)
  hM$repList = vector("list", nChains)
  if (!exists(".Random.seed")) 
    runif(1)
  hM$randSeed <- .Random.seed
  initSeed = sample.int(.Machine$integer.max, nChains)
  EPS = 1e-06
  updaterWarningFlag = TRUE
  if (!identical(updater$Gamma2, FALSE) && any(abs(mGamma) > 
                                               EPS)) {
    updater$Gamma2 = FALSE
    if (updaterWarningFlag) 
      message("setting updater$Gamma2=FALSE due to non-zero mGamma")
  }
  if (!identical(updater$Gamma2, FALSE) && any(abs(iUGamma - 
                                                   kronecker(iUGamma[1:hM$nc, 1:hM$nc], diag(hM$nt))) > 
                                               EPS)) {
    updater$Gamma2 = FALSE
    if (updaterWarningFlag) 
      message("setting updater$Gamma2=FALSE due to non-kronecker structure of UGamma matrix")
  }
  if (!identical(updater$Gamma2, FALSE) && (!is.null(C))) {
    updater$Gamma2 = FALSE
    if (updaterWarningFlag) 
      message("setting updater$Gamma2=FALSE due to specified phylogeny matrix")
  }
  if (!identical(updater$GammaEta, FALSE) && any(abs(mGamma) > 
                                                 EPS)) {
    updater$GammaEta = FALSE
    if (updaterWarningFlag) 
      message("setting updater$GammaEta=FALSE due to non-zero mGamma")
  }
  if (!identical(updater$GammaEta, FALSE) && hM$nr == 0) {
    updater$GammaEta = FALSE
    if (updaterWarningFlag) 
      message("setting updater$GammaEta=FALSE due to absence of random effects included to the model")
  }
  if (!identical(updater$GammaEta, FALSE) && any(sapply(hM$rL, 
                                                        function(s) !is.null(s$spatialMethod) && s$spatialMethod %in% 
                                                        c("GPP", "NNGP")))) {
    updater$GammaEta = FALSE
    if (updaterWarningFlag) 
      message("setting updater$GammaEta=FALSE: not implemented for spatial methods 'GPP' and 'NNGP'")
  }
  if (identical(updater$latentLoadingOrderSwap, NULL)) {
    updater$latentLoadingOrderSwap = 0
    if (FALSE && updaterWarningFlag) 
      message("setting updater$latentLoadingOrderSwap=0 disabling full-conditional swapping of consecutive latent loadings")
  }
  sampleChain = function(chain) {
    if (nChains > 1) 
      cat(sprintf("Computing chain %d\n", chain))
    set.seed(initSeed[chain])
    parList = computeInitialParameters(hM, initPar)
    Gamma = parList$Gamma
    V = parList$V
    iV = chol2inv(chol(V))
    Beta = parList$Beta
    BetaSel = parList$BetaSel
    PsiRRR = parList$PsiRRR
    DeltaRRR = parList$DeltaRRR
    wRRR = parList$wRRR
    sigma = parList$sigma
    iSigma = 1/sigma
    Lambda = parList$Lambda
    Eta = parList$Eta
    Alpha = parList$Alpha
    Psi = parList$Psi
    Delta = parList$Delta
    rho = parList$rho
    Z = parList$Z
    X1A = X1
    if (hM$ncsel > 0) {
      for (i in 1:hM$ncsel) {
        XSel = hM$XSelect[[i]]
        for (spg in 1:length(XSel$q)) {
          if (!BetaSel[[i]][spg]) {
            fsp = which(XSel$spGroup == spg)
            for (j in fsp) {
              X1A[[j]][, XSel$covGroup] = 0
            }
          }
        }
      }
    }
    X = X1A
    if (hM$ncRRR > 0) {
      XB = hM$XRRRScaled %*% t(wRRR)
      if (is.matrix(X)) {
        X = cbind(X, XB)
      }
      else {
        for (j in 1:hM$ns) {
          X[[j]] = cbind(X[[j]], XB)
        }
      }
    }
    if (!identical(updater$Gamma2, FALSE) && (!is.matrix(X))) {
      updater$Gamma2 = FALSE
      if (updaterWarningFlag) 
        message("setting updater$Gamma2=FALSE due to X is not a matrix")
    }
    if (!identical(updater$GammaEta, FALSE) && (!is.matrix(X))) {
      updater$GammaEta = FALSE
      if (updaterWarningFlag) 
        message("setting updater$GammaEta=FALSE due to X is not a matrix")
    }
    postList = vector("list", samples)
    failed <- numeric(14)
    names(failed) <- c("Gamma2", "GammaEta", 
                       "BetaLambda", "wRRR", "BetaSel", 
                       "GammaV", "Rho", "LambdaPriors", 
                       "wRRRPriors", "Eta", "Alpha", "invSigma", 
                       "Nf", "LatentLoadingOrder")
    for (iter in seq_len(transient + samples * thin)) {
      if (!identical(updater$Gamma2, FALSE)) {
        out = try(updateGamma2(Z = Z, Gamma = Gamma, 
                               iV = iV, iSigma = iSigma, Eta = Eta, Lambda = Lambda, 
                               X = X, Pi = Pi, dfPi = dfPi, Tr = Tr, C = C, 
                               rL = hM$rL, iQg = iQg, mGamma = mGamma, iUGamma = iUGamma), 
                  silent = TRUE)
        if (!inherits(out, "try-error")) {
          Gamma <- out
        }
        else {
          failed["Gamma2"] <- failed["Gamma2"] + 
            1
        }
      }
      if (!identical(updater$GammaEta, FALSE)) {
        GammaEtaList = try(updateGammaEta(Z = Z, Gamma = Gamma, 
                                          V = chol2inv(chol(iV)), iV = iV, id = iSigma, 
                                          Eta = Eta, Lambda = Lambda, Alpha = Alpha, 
                                          X = X, Pi = Pi, dfPi = dfPi, Tr = Tr, rL = hM$rL, 
                                          rLPar = rLPar, Q = Qg[, , rho], iQ = iQg[, 
                                                                                   , rho], RQ = RQg[, , rho], mGamma = mGamma, 
                                          U = hM$UGamma, iU = iUGamma), silent = TRUE)
        if (!inherits(GammaEtaList, "try-error")) {
          Gamma = GammaEtaList$Gamma
          Eta = GammaEtaList$Eta
        }
        else {
          failed["GammaEta"] <- failed["GammaEta"] + 
            1
        }
      }
      if (!identical(updater$BetaLambda, FALSE)) {
        BetaLambdaList = try(updateBetaLambda(Y = Y, 
                                              Z = Z, Gamma = Gamma, iV = iV, iSigma = iSigma, 
                                              Eta = Eta, Psi = Psi, Delta = Delta, iQ = iQg[, 
                                                                                            , rho], X = X, Tr = Tr, Pi = Pi, dfPi = dfPi, 
                                              C = C, rL = hM$rL), silent = TRUE)
        if (!inherits(BetaLambdaList, "try-error")) {
          Beta = BetaLambdaList$Beta
          Lambda = BetaLambdaList$Lambda
        }
        else {
          failed["BetaLambda"] <- failed["BetaLambda"] + 
            1
        }
      }
      if (!identical(updater$wRRR, FALSE) && hM$ncRRR > 
          0) {
        wRRRXList = try(updatewRRR(Z = Z, Beta = Beta, 
                                   iSigma = iSigma, Eta = Eta, Lambda = Lambda, 
                                   X1A = X1A, XRRR = hM$XRRRScaled, Pi = Pi, dfPi = dfPi, 
                                   rL = hM$rL, PsiRRR = PsiRRR, DeltaRRR = DeltaRRR), 
                        silent = TRUE)
        if (!inherits(wRRRXList, "try-error")) {
          wRRR = wRRRXList$wRRR
          X = wRRRXList$X
        }
        else {
          failed["wRRR"] <- failed["wRRR"] + 
            1
        }
      }
      if (!identical(updater$BetaSel, FALSE) && hM$ncsel > 
          0) {
        BetaSelXList = try(updateBetaSel(Z = Z, XSelect = hM$XSelect, 
                                         BetaSel = BetaSel, Beta = Beta, iSigma = iSigma, 
                                         Lambda = Lambda, Eta = Eta, X1 = X1, Pi = Pi, 
                                         dfPi = dfPi, rL = hM$rL), silent = TRUE)
        if (!inherits(BetaSelXList, "try-error")) {
          BetaSel = BetaSelXList$BetaSel
          X = BetaSelXList$X
        }
        else {
          failed["BetaSel"] <- failed["BetaSel"] + 
            1
        }
      }
      if (!identical(updater$GammaV, FALSE)) {
        GammaVList = try(updateGammaV(Beta = Beta, Gamma = Gamma, 
                                      iV = iV, rho = rho, iQg = iQg, RQg = RQg, Tr = Tr, 
                                      C = C, mGamma = mGamma, iUGamma = iUGamma, 
                                      V0 = V0, f0 = f0), silent = TRUE)
        if (!inherits(GammaVList, "try-error")) {
          Gamma = GammaVList$Gamma
          iV = GammaVList$iV
        }
        else {
          failed["GammaV"] <- failed["GammaV"] + 
            1
        }
      }
      if (!is.null(hM$C) && !identical(updater$Rho, FALSE)) {
        out = try(updateRho(Beta = Beta, Gamma = Gamma, 
                            iV = iV, RQg = RQg, detQg = detQg, Tr = Tr, 
                            rhopw = rhopw), silent = TRUE)
        if (!inherits(out, "try-error")) 
          rho <- out
        else failed["Rho"] <- failed["Rho"] + 
            1
      }
      if (!identical(updater$LambdaPriors, FALSE)) {
        PsiDeltaList = try(updateLambdaPriors(Lambda = Lambda, 
                                              Delta = Delta, rL = hM$rL), silent = TRUE)
        if (!inherits(PsiDeltaList, "try-error")) {
          Psi = PsiDeltaList$Psi
          Delta = PsiDeltaList$Delta
        }
        else {
          failed["LambdaPriors"] <- failed["LambdaPriors"] + 
            1
        }
      }
      if (!identical(updater$wRRRPriors, FALSE) && hM$ncRRR > 
          0) {
        PsiDeltaList = try(updatewRRRPriors(wRRR = wRRR, 
                                            Delta = DeltaRRR, nu = hM$nuRRR, a1 = hM$a1RRR, 
                                            b1 = hM$b1RRR, a2 = hM$a2RRR, b2 = hM$b2RRR), 
                           silent = TRUE)
        if (!inherits(PsiDeltaList, "try-error")) {
          PsiRRR = PsiDeltaList$Psi
          DeltaRRR = PsiDeltaList$Delta
        }
        else {
          failed["wRRRPriors"] <- failed["wRRRPriors"] + 
            1
        }
      }
      if (!identical(updater$Eta, FALSE)) 
        out = try(updateEta(Y = Y, Z = Z, Beta = Beta, 
                            iSigma = iSigma, Eta = Eta, Lambda = Lambda, 
                            Alpha = Alpha, rLPar = rLPar, X = X, Pi = Pi, 
                            dfPi = dfPi, rL = hM$rL), silent = TRUE)
      if (!inherits(out, "try-error")) 
        Eta <- out
      else failed["Eta"] <- failed["Eta"] + 
          1
      if (!identical(updater$Alpha, FALSE)) 
        out = try(updateAlpha(Eta = Eta, rLPar = rLPar, 
                              rL = hM$rL), silent = TRUE)
      if (!inherits(out, "try-error")) 
        Alpha <- out
      else failed["Alpha"] <- failed["Alpha"] + 
          1
      if (!identical(updater$InvSigma, FALSE)) 
        out = try(updateInvSigma(Y = Y, Z = Z, Beta = Beta, 
                                 iSigma = iSigma, Eta = Eta, Lambda = Lambda, 
                                 distr = distr, X = X, Pi = Pi, dfPi = dfPi, 
                                 rL = hM$rL, aSigma = aSigma, bSigma = bSigma), 
                  silent = TRUE)
      if (!inherits(out, "try-error")) 
        iSigma <- out
      else failed["invSigma"] <- failed["invSigma"] + 
          1
      if (!identical(updater$Z, FALSE)) {
        Z = updateZ(Y = Y, Z = Z, Beta = Beta, iSigma = iSigma, 
                    Eta = Eta, Lambda = Lambda, X = X, Pi = Pi, 
                    dfPi = dfPi, distr = distr, rL = hM$rL)
      }
      for (r in seq_len(nr)) {
        if (iter <= adaptNf[r]) {
          listPar = try(updateNf(eta = Eta[[r]], lambda = Lambda[[r]], 
                                 alpha = Alpha[[r]], psi = Psi[[r]], delta = Delta[[r]], 
                                 rL = hM$rL[[r]], iter = iter), silent = TRUE)
          if (!inherits(listPar, "try-error")) {
            Lambda[[r]] = listPar$lambda
            Eta[[r]] = listPar$eta
            Alpha[[r]] = listPar$alpha
            Psi[[r]] = listPar$psi
            Delta[[r]] = listPar$delta
          }
          else {
            failed["Nf"] <- failed["Nf"] + 
              1
          }
        }
      }
      if (updater$latentLoadingOrderSwap > 0 && (iter%%updater$latentLoadingOrderSwap == 
                                                 0)) {
        for (r in seq_len(nr)) {
          listPar = try(updateLatentLoadingOrder(eta = Eta[[r]], 
                                                 lambda = Lambda[[r]], alpha = Alpha[[r]], 
                                                 delta = Delta[[r]], rL = hM$rL[[r]]), silent = TRUE)
          if (!inherits(listPar, "try-error")) {
            Lambda[[r]] = listPar$lambda
            Eta[[r]] = listPar$eta
            Alpha[[r]] = listPar$alpha
            Delta[[r]] = listPar$delta
          }
          else {
            failed["LatentLoadingOrder"] + failed["LatentLoadingOrder"] + 
              1
          }
        }
        PsiDeltaList = try(updateLambdaPriors(Lambda = Lambda, 
                                              Delta = Delta, rL = hM$rL))
        if (!inherits(PsiDeltaList, "try-error")) {
          Psi = PsiDeltaList$Psi
          Delta = PsiDeltaList$Delta
        }
        else {
          failed["PsiDelta"] <- failed["PsiDelta"] + 
            1
        }
      }
      if ((iter > transient) && ((iter - transient)%%thin == 
                                 0)) {
        postList[[(iter - transient)/thin]] = combineParameters(Beta = Beta, 
                                                                BetaSel = BetaSel, wRRR = wRRR, Gamma = Gamma, 
                                                                iV = iV, rho = rho, iSigma = iSigma, Eta = Eta, 
                                                                Lambda = Lambda, Alpha = Alpha, Psi = Psi, 
                                                                Delta = Delta, PsiRRR = PsiRRR, DeltaRRR = DeltaRRR, 
                                                                ncNRRR = hM$ncNRRR, ncRRR = hM$ncRRR, ncsel = hM$ncsel, 
                                                                XSelect = hM$XSelect, XScalePar = hM$XScalePar, 
                                                                XInterceptInd = hM$XInterceptInd, XRRRScalePar = hM$XRRRScalePar, 
                                                                nt = hM$nt, TrScalePar = hM$TrScalePar, TrInterceptInd = hM$TrInterceptInd, 
                                                                rhopw = rhopw)
      }
      postList$failedUpdates <- failed
      if ((verbose > 0) && (iter%%verbose == 0)) {
        if (iter > transient) {
          samplingStatusString = "sampling"
        }
        else {
          samplingStatusString = "transient"
        }
        cat(sprintf("Chain %d, iteration %d of %d (%s)\n", 
                    chain, iter, transient + samples * thin, samplingStatusString))
      }
    }
    return(postList)
  }
  if (nParallel > 1) {
    if (useSocket) {
      cl = makeCluster(nParallel)
      clusterExport(cl, c("hM", "nChains", 
                          "transient", "samples", "thin", 
                          "verbose", "adaptNf", "initSeed", 
                          "initPar", "updater", "X1", 
                          "Tr", "Y", "distr", "Pi", 
                          "C", "nr", "mGamma", "iUGamma", 
                          "V0", "f0", "aSigma", "bSigma", 
                          "rhopw", "Qg", "iQg", "RQg", 
                          "detQg", "rLPar"), envir = environment())
      clusterEvalQ(cl, {
        library(BayesLogit)
        library(MCMCpack)
        library(truncnorm)
        library(Matrix)
        library(abind)
        library(Hmsc)
      })
      hM$postList = clusterApplyLB(cl, 1:nChains, fun = sampleChain)
      stopCluster(cl)
    }
    else {
      verbose <- 0
      hM$postList <- mclapply(seq_len(nChains), function(i) sampleChain(i), 
                              mc.cores = nParallel)
    }
  }
  else {
    for (chain in 1:nChains) {
      if (fromPrior) {
        postList = vector("list", samples)
        for (iter in 1:samples) {
          postList[[iter]] = samplePrior(hM, dataParList = dataParList)
        }
        hM$postList[[chain]] = postList
      }
      else {
        hM$postList[[chain]] = sampleChain(chain)
      }
    }
  }
  for (chain in seq_len(nChains)) {
    ntries <- transient + samples * thin
    if (any(hM$postList[[chain]]$failedUpdates > 0)) {
      cat("Failed updaters and their counts in chain", 
          chain, " (", ntries, " attempts):\n")
      failures <- hM$postList[[chain]]$failedUpdates
      failures <- failures[failures > 0]
      print(failures)
    }
    attr(hM$postList[[chain]], "failedUpdates") <- hM$postList[[chain]]$failedUpdates
    hM$postList[[chain]]$failedUpdates <- NULL
  }
  hM$samples = samples
  hM$transient = transient
  hM$thin = thin
  hM$verbose = verbose
  hM$adaptNf = adaptNf
  if (alignPost) {
    for (i in 1:5) {
      hM = alignPosterior(hM)
    }
  }
  return(hM)
}
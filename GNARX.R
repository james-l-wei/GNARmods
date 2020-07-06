# This file contains code for the following functions:
# 1) GNARXfit: fits the GNARX model via OLS, allowing for +ve coefficient sign restriction
# 2) GNARXdesign: creates the design matrix for the OLS regression in GNARXfit
# 3) GNARXorder: computes the optimum GNARX model order based on greedy BIC optimisation
# 4) oneStepMSFE: computes the one-step ahead OOS MSFE for a given GNARX specification*
# 5) GNARXpredict: computes the one-step ahead prediction for a given GNARX specification*
# 6) VAR.oneStepMSFE: computes the one-step ahead OOS MSFE for a given VAR specification
# 7) inSampleError: computes the in-sample MSE for a given GNARX specification*2
# * Haven't yet been tested with a non-null xvts parameter
# *2 Doesn't yet work with positive coefficient restriction; need to remove first coefficient

library(matrixcalc)
library(GNAR)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

GNARXfit <- function (vts = GNAR::fiveVTS, xvts = NULL, xvts2 = NULL, net = GNAR::fiveNet, 
                      alphaOrder = 2, betaOrder = c(1, 1), lambdaOrder = NULL,
                      lambdaOrder2 = NULL, fact.var = NULL, globalalpha = TRUE, 
                      tvnets = NULL, netsstart = NULL, ErrorIfNoNei = TRUE, 
                      positiveCoef = FALSE) 
{
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if (!is.null(fact.var)) {
    stopifnot(length(fact.var) == length(net$edges))
  }
  stopifnot(is.matrix(vts))
  stopifnot(is.logical(globalalpha))
  if (!is.null(tvnets)) {
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  useNofNei <- 1
  if(is.null(lambdaOrder)){
    if(!is.null(xvts)){
      stop("lambdaOrder must be specified when xvts is provided")
    }else{
      maxOrder <- alphaOrder
    }
  }else{
    if((nrow(xvts) != nrow(vts))){
      stop("xvts and vts must have the same time range")
    }else{
      if(lambdaOrder%%1 != 0 || lambdaOrder < 0){
        stop("lambdaOrder must be a non-negative integer")
      }else{
        if(is.null(lambdaOrder2)){
          maxOrder <- max(alphaOrder, lambdaOrder)
        }else{
          if(lambdaOrder2%%1 != 0 || lambdaOrder2 < 0){
            stop("lambdaOrder2 must be a non-negative integer")
          }else{
            maxOrder <- max(alphaOrder, lambdaOrder, lambdaOrder2)
          }
        }
      }
    }
  }
  if(is.null(xvts) && !is.null(xvts2)){
    stop("xvts cannot be null when xvts2 is not null")
  }
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder, 
                betas.in = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, 
                final.in = vts[(nrow(vts) - maxOrder + 1):nrow(vts), ], 
                lambdas.in = lambdaOrder, positiveCoef = positiveCoef,
                lambdas2.in = lambdaOrder2)
  dmat <- GNARXdesign(vts = vts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart, 
                     xvts = xvts, lambdaOrder = lambdaOrder, xvts2 = xvts2, 
                     lambdaOrder2 = lambdaOrder2)
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 
                                            2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  predt <- nrow(vts) - maxOrder
  yvec <- NULL
  for (ii in 1:length(net$edges)) {
    yvec <- c(yvec, vts[((maxOrder + 1):(predt + maxOrder)), 
                        ii])
  }
  dmat2 <- dmat[complete.cases(dmat), ]
  yvec2 <- yvec[complete.cases(dmat)]
  if(sum(!is.na(yvec2)) > 0){
    if(ncol(dmat) != 1){
      dmat2 <- dmat2[!is.na(yvec2), ]
    }else{
      dmat2 <- dmat2[!is.na(yvec2)]
    }
    yvec2 <- yvec2[!is.na(yvec2)]
  }
  if(positiveCoef == FALSE){
    modNoIntercept <- lm(yvec2 ~ dmat2 + 0)
    modBIC <- BIC(modNoIntercept)
  }else{
    if(is.null(dim(dmat2))){
      modNoIntercept <- lm(yvec2 ~ dmat2 + 0)
      modBIC <- NA
    }else{
      modNoIntercept <- glmnet(dmat2, yvec2, lambda = 0, lower.limits = 0, intercept = FALSE)
      tLL <- modNoIntercept$nulldev - deviance(modNoIntercept)
      k <- modNoIntercept$df
      n <- modNoIntercept$nobs
      modBIC <- log(n) * k - tLL
    }
  }
  out <- list(mod = modNoIntercept, y = yvec, dd = dmat, frbic = frbic, BIC = modBIC)
  class(out) <- "GNARXfit"
  return(out)
}


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

GNARXdesign <- function (vts = GNAR::fiveVTS, xvts = NULL, xvts2 = NULL, net = GNAR::fiveNet, 
                         alphaOrder = 2, betaOrder = c(1, 1), lambdaOrder = NULL,
                         lambdaOrder2 = NULL, fact.var = NULL, globalalpha = TRUE, 
                         tvnets = NULL, netsstart = NULL) 
{
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if (!is.null(fact.var)) {
    stopifnot(length(fact.var) == length(net$edges))
    if (!globalalpha) {
      stop("Use factors OR individual alphas")
    }
  }
  stopifnot(is.matrix(vts))
  stopifnot(is.logical(globalalpha))
  if (!is.null(tvnets)) {
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  if(is.null(lambdaOrder)){
    if(!is.null(xvts)){
      stop("lambdaOrder must be specified when xvts is provided")
    }else{
      maxOrder <- alphaOrder
    }
  }else{
    if(nrow(xvts) != nrow(vts)){
      stop("xvts and vts must have the same time range")
    }else{
      if(lambdaOrder%%1 != 0 || lambdaOrder < 0){
        stop("lambdaOrder must be a non-negative integer")
      }else{
        if(is.null(lambdaOrder2)){
          maxOrder <- max(alphaOrder, lambdaOrder)
        }else{
          maxOrder <- max(alphaOrder, lambdaOrder, lambdaOrder2)
        }
      }
    }
  }
  netmat <- as.matrix(net, normalise = FALSE)
  if (!isSymmetric(netmat)) {
    net <- as.GNARnet(t(netmat))
  }
  parNames <- parLoc <- NULL
  for (jj in 1:alphaOrder) {
    if (globalalpha) {
      parNames <- c(parNames, paste("alpha", jj, sep = ""))
      parLoc <- c(parLoc, "a")
    }
    else {
      for (kk in 1:ncol(vts)) {
        parNames <- c(parNames, paste("alpha", jj, "node", 
                                      kk, sep = ""))
        parLoc <- c(parLoc, "a")
      }
    }
    if (betaOrder[jj] > 0) {
      for (kk in 1:betaOrder[jj]) {
        parNames <- c(parNames, paste("beta", jj, ".", 
                                      kk, sep = ""))
        parLoc <- c(parLoc, "b")
      }
    }
  }
  if(!is.null(lambdaOrder)){
    for(ii in 0:lambdaOrder){
      parNames <- c(parNames, paste("lambda.", ii, sep = ""))
      parLoc <- c(parLoc, "l1") 
    }
  }
  if(!is.null(lambdaOrder2)){
    for(ii in 0:lambdaOrder2){
      parNames <- c(parNames, paste("lambda2.", ii, sep = ""))
      parLoc <- c(parLoc, "l2") 
    }
  }
  predt <- nrow(vts) - maxOrder
  nnodes <- ncol(vts)
  if (globalalpha) {
    dmat <- matrix(0, nrow = predt * nnodes, ncol = length(parLoc), 
                   dimnames = list(NULL, parNames))
  }else {
    dmat <- matrix(0, nrow = predt * nnodes, ncol = length(parLoc), 
                   dimnames = list(NULL, parNames))
  }
  for (ii in 1:nnodes) {
    for (aa in 1:alphaOrder) {
      if (globalalpha) {
        alphaLoc <- which(parLoc == "a")[aa]
      } else {
        alphaLoc <- which(parLoc == "a")[nnodes * (aa - 1) + ii]
      }
      dmat[((predt * (ii - 1) + 1):(predt * ii)), alphaLoc] <- 
        vts[((maxOrder + 1 - aa):(predt + (maxOrder - aa))), ii]
    }
  }
  if (sum(betaOrder) > 0) {
    betaN <- NULL
    betaTimes <- rep(1:alphaOrder, betaOrder)
    for (jj in 1:alphaOrder) {
      if (betaOrder[jj] > 0) {
        betaN <- c(betaN, 1:betaOrder[jj])
      }
    }
    for (ii in 1:nnodes) {
      NofNei <- NofNeighbours(node = ii, stage = max(betaOrder), net = net)
      Nei <- NofNei$edges
      Wei <- NofNei$dist
      if ((!is.null(Nei)) & (length(Nei) > 0)) {
        if (!is.null(Nei[[1]]) & !is.na(Nei[[1]][1])) {
          Wei <- lapply(Wei, function(x) {
            1/(x * sum(1/x))
          })
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            if (length(Nei[[betaN[bb]]]) > 1) {
              vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):
                                (predt + (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]
              for (kk in 1:nrow(vts.cut)) {
                if (any(is.na(vts.cut[kk, ]))) {
                  if (all(is.na(vts.cut[kk, ]))) {
                    vts.cut[kk, ] <- 0
                  } else {
                    new.wei <- Wei[[betaN[bb]]][which(!is.na(vts.cut[kk, ]))]
                    new.wei <- new.wei/sum(new.wei)
                    sub.val <- vts.cut[kk, which(!is.na(vts.cut[kk, ]))] %*% new.wei
                    vts.cut[kk, which(is.na(vts.cut[kk, ]))] <- sub.val
                  }
                }
              }
              dmat[((predt * (ii - 1) + 1):(predt * ii)), 
                   betaLoc] <- vts.cut %*% Wei[[betaN[bb]]]
            } else {
              if ((length(Nei[[betaN[bb]]]) == 1) & (!is.na(Nei[[betaN[bb]]]))) {
                vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):
                                  (predt + (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]
                vts.cut[is.na(vts.cut)] <- 0
                dmat[((predt * (ii - 1) + 1):(predt * ii)), betaLoc] <- vts.cut * 
                  Wei[[betaN[bb]]]
              } else {
                dmat[((predt * (ii - 1) + 1):(predt * ii)), betaLoc] <- 0
              }
            }
          }
        } else {
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            dmat[((predt * (ii - 1) + 1):(predt * ii)), betaLoc] <- 0
          }
        }
      } else {
        for (bb in 1:sum(betaOrder)) {
          betaLoc <- which(parLoc == "b")[bb]
          dmat[((predt * (ii - 1) + 1):(predt * ii)), betaLoc] <- 0
        }
      }
    }
  }
  if(!is.null(xvts)){
    if(is.null(xvts2)){
      for(jj in 0:lambdaOrder){
        dmat[, ncol(dmat) - lambdaOrder + jj] <- vec(xvts[(maxOrder + 1 - jj):
                                                            (nrow(xvts) - jj), ])
      }
    }else{
      for(jj in 0:lambdaOrder){
        dmat[, ncol(dmat) - lambdaOrder2 - 1 - lambdaOrder + jj] <- vec(xvts[(maxOrder + 1 - jj):
                                                                               (nrow(xvts) - jj), ])
      }
      for(jj in 0:lambdaOrder2){
        dmat[, ncol(dmat) - lambdaOrder2 + jj] <- vec(xvts2[(maxOrder + 1 - jj):
                                                              (nrow(xvts2) - jj), ])
      }
    }
  }
  if (is.null(fact.var)) {
    return(dmat)
  } else {
    stop("Non-null fact.var has not been implemented for GNARX")
    facun <- unique(fact.var)
    if (length(facun) == 1) {
      return(dmat)
    } else {
      dmcol <- ncol(dmat)
      dmatex <- dmat
      exnames <- paste(colnames(dmat), " '", facun[1], "'", sep = "")
      for (ii in 2:length(facun)) {
        dmatex <- cbind(dmatex, dmat)
        exnames <- c(exnames, paste(colnames(dmat), " '", facun[ii], "'", sep = ""))
      }
      for (ii in 1:length(facun)) {
        dmatex[fact.var != facun[ii], ((ii - 1) * dmcol + (1:dmcol))] <- 0
      }
      colnames(dmatex) <- exnames
      return(dmatex)
    }
  }
}


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

GNARXOrder <- function(vts, xvts = NULL, xvts2 = NULL, net, alphaLim, 
                       lambdaLim = NULL, lambdaLim2 = NULL, globalalpha, 
                       positiveCoef = FALSE){
  if(is.null(xvts)){
    bestLambda <- NULL
  }else{
    bestLambda <- 0
  }
  if(is.null(xvts2)){
    bestLambda2 <- NULL
  }else{
    bestLambda2 <- 0
  }
  if(positiveCoef){
    bestAlpha <- 1
    newBIC <- GNARXfit(vts = vts, net = net, alphaOrder = 2, 
                       betaOrder = c(0,0), globalalpha = globalalpha, xvts = xvts, 
                       lambdaOrder = bestLambda, positiveCoef = positiveCoef, xvts2 = xvts2,
                       lambdaOrder2 = bestLambda2)$BIC
    bestBIC <- newBIC + 1
  }else{
    bestAlpha <- 0
    newBIC <- GNARXfit(vts = vts, net = net, alphaOrder = 1, 
                       betaOrder = 0, globalalpha = globalalpha, xvts = xvts, 
                       lambdaOrder = bestLambda, positiveCoef = positiveCoef, xvts2 = xvts2,
                       lambdaOrder2 = bestLambda2)$BIC
    bestBIC <- newBIC + 1
  }
  while(newBIC < bestBIC & bestAlpha < alphaLim){ # alpha search
    bestBIC <- newBIC
    bestAlpha <- bestAlpha + 1
    newBIC <- GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha + 1, 
                       betaOrder = rep(0, bestAlpha + 1), globalalpha = globalalpha, 
                       xvts = xvts, lambdaOrder = bestLambda, 
                       positiveCoef = positiveCoef, xvts2 = xvts2,
                       lambdaOrder2 = bestLambda2)$BIC
  }
  bestBeta <- rep(0, bestAlpha)
  bestBeta[1] <- -1 
  for(i in 1:bestAlpha){
    newBIC <- bestBIC
    bestBIC <- newBIC + 1
    while(newBIC < bestBIC){ # beta search
      bestBIC <- newBIC
      bestBeta[i] <- bestBeta[i] + 1
      newBeta <- bestBeta
      newBeta[i] <- bestBeta[i] + 1
      newBIC <- try(GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha, 
                             betaOrder = newBeta, globalalpha = globalalpha, 
                             xvts = xvts, lambdaOrder = bestLambda, 
                             positiveCoef = positiveCoef, xvts2 = xvts2,
                             lambdaOrder2 = bestLambda2)$BIC, silent = TRUE)
      if(inherits(newBIC, "try-error")){
        newBIC <- bestBIC + 1
      }
    }
    if(bestAlpha >= i+1){
      bestBeta[i+1] <- -1 
    }
  }
  if(is.null(xvts)){
    bestOrder <- list(bestBIC = bestBIC, bestAlpha = bestAlpha, bestBeta = bestBeta)
  }else{
    newBIC <- bestBIC
    bestBIC <- newBIC + 1
    bestLambda <- -1
    while(newBIC < bestBIC & bestLambda < lambdaLim){ # lambda search
      bestBIC <- newBIC
      bestLambda <- bestLambda + 1
      newBIC <- GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha, 
                         betaOrder = bestBeta, globalalpha = globalalpha, 
                         xvts = xvts, lambdaOrder = bestLambda + 1, 
                         positiveCoef = positiveCoef, xvts2 = xvts2,
                         lambdaOrder2 = bestLambda2)$BIC
    }
    if(is.null(xvts2)){
      bestOrder <- list(bestBIC = bestBIC, bestAlpha = bestAlpha, bestBeta = bestBeta, 
                        bestLambda = bestLambda)
    }else{
      newBIC <- bestBIC
      bestBIC <- newBIC + 1
      bestLambda2 <- -1
      while(newBIC < bestBIC & bestLambda2 < lambdaLim2){ # lambda search
        bestBIC <- newBIC
        bestLambda2 <- bestLambda2 + 1
        newBIC <- GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha, 
                           betaOrder = bestBeta, globalalpha = globalalpha, 
                           xvts = xvts, lambdaOrder = bestLambda, 
                           positiveCoef = positiveCoef, xvts2 = xvts2,
                           lambdaOrder2 = bestLambda2 + 1)$BIC
      }
      bestOrder <- list(bestBIC = bestBIC, bestAlpha = bestAlpha, bestBeta = bestBeta, 
                        bestLambda = bestLambda, bestLambda2 = bestLambda2)
    }
  }
  return(bestOrder)
}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

oneStepMSFE <- function(old.vts, old.xvts = NULL, old.xvts2 = NULL, new.vts, 
                        new.xvts = NULL, new.xvts2 = NULL, net, 
                        alphaOrder, betaOrder, lambdaOrder = NULL, lambdaOrder2 = NULL, 
                        globalalpha, positiveCoef = FALSE){
  noPreds <- nrow(new.vts) - nrow(old.vts)
  predMat <- matrix(NA, nrow = noPreds, ncol = ncol(old.vts))
  for(i in 1:noPreds){
    if(is.null(old.xvts)){
      fit <- GNARXfit(vts = new.vts[1:(nrow(old.vts) + i - 1),], net = net, 
                      alphaOrder = alphaOrder, betaOrder = betaOrder,
                      globalalpha = globalalpha, positiveCoef = positiveCoef)
      predMat[i, ] <- GNARXpredict(fit = fit, new.vts = new.vts[1:(nrow(old.vts) + i),]) 
    }else{
      if(is.null(old.xvts2)){
        fit <- GNARXfit(vts = new.vts[1:(nrow(old.vts) + i - 1),], 
                        xvts = new.xvts[1:(nrow(old.xvts) + i - 1),], 
                        net = net, alphaOrder = alphaOrder, betaOrder = betaOrder, 
                        lambdaOrder = lambdaOrder, globalalpha = globalalpha, 
                        positiveCoef = positiveCoef)
        predMat[i, ] <- GNARXpredict(fit = fit, new.vts = new.vts[1:(nrow(old.vts) + i),], 
                                     new.xvts = new.xvts[1:(nrow(old.xvts) + i),]) 
      }else{
        fit <- GNARXfit(vts = new.vts[1:(nrow(old.vts) + i - 1),], 
                        xvts = new.xvts[1:(nrow(old.xvts) + i - 1),],
                        xvts2 = new.xvts[1:(nrow(old.xvts2) + i - 1),],
                        net = net, alphaOrder = alphaOrder, betaOrder = betaOrder, 
                        lambdaOrder = lambdaOrder, lambdaOrder2 = lambdaOrder2,
                        globalalpha = globalalpha, positiveCoef = positiveCoef)
        predMat[i, ] <- GNARXpredict(fit = fit, new.vts = new.vts[1:(nrow(old.vts) + i),], 
                                     new.xvts = new.xvts[1:(nrow(old.xvts) + i),], 
                                     new.xvts2 = new.xvts2[1:(nrow(old.xvts2) + i),]) 
      }
    }
  }
  MSFE <- mean((predMat - new.vts[(nrow(old.vts) + 1):nrow(new.vts),])^2, na.rm = TRUE)
  return(list(MSFE = MSFE, predMat = predMat))
}


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

GNARXpredict <- function(fit, new.vts, new.xvts = NULL, new.xvts2 = NULL){
  nnodes <- fit$frbic$nnodes
  dmat <- GNARXdesign(vts = new.vts, xvts = new.xvts, xvts2 = new.xvts2,
                      net = fit$frbic$net.in, alphaOrder = fit$frbic$alphas.in, 
                      betaOrder = fit$frbic$betas.in, lambdaOrder = fit$frbic$lambdas.in, 
                      lambdaOrder2 = fit$frbic$lambdas2.in, fact.var = fit$frbic$fact.var, 
                      globalalpha = fit$frbic$globalalpha, tvnets = fit$frbic$tvnets, 
                      netsstart = fit$frbic$netsstart)
  gammaHat <- coef(fit$mod)
  if(fit$frbic$positiveCoef){
    gammaHat <- gammaHat[-1, ]
  }
  yvecHat <- dmat %*% gammaHat
  yvec <- matrix(yvecHat, nrow = length(yvecHat) / nnodes, ncol = nnodes)
  predictions <- tail(yvec, 1)
  return(predictions)
}


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

VAR.oneStepMSFE <- function(old.vts, new.vts, order){
  old.vts <- old.vts[complete.cases(old.vts), ]
  new.vts <- new.vts[complete.cases(new.vts), ]
  noPreds <- nrow(new.vts) - nrow(old.vts)
  predMat <- matrix(NA, nrow = noPreds, ncol = ncol(old.vts))
  for(i in 1:noPreds){
    fit <- VAR(new.vts[1:(nrow(old.vts) + i - 1),], p = order, include.mean = FALSE,
               output = FALSE)
    predMat[i, ] <- VARpred(fit, h=1)$pred
  }
  MSFE <- mean((predMat - new.vts[(nrow(old.vts) + 1):nrow(new.vts),])^2, na.rm = TRUE)
  return(list(MSFE = MSFE, predMat = predMat))
}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

inSampleError <- function(vts, xvts = NULL, xvts2 = NULL, net, alphaOrder, betaOrder, 
                          lambdaOrder = NULL, lambdaOrder2 = NULL, globalalpha, positiveCoef){
  if(positiveCoef){
    stop("positiveCoef = TRUE not supported yet")
  }
  nnodes <- ncol(vts)
  mod <- GNARXfit(vts = vts, xvts = xvts, xvts2 = xvts2, net = net, 
                  alphaOrder = alphaOrder, betaOrder = betaOrder, lambdaOrder = lambdaOrder,
                  lambdaOrder2 = lambdaOrder2, globalalpha = globalalpha, 
                  positiveCoef = positiveCoef)
  fittedVals <- mod$dd %*% coef(mod$mod)
  fittedVals <- matrix(fittedVals, ncol = nnodes)
  MSE <- mean((fittedVals - tail(vts, nrow(fittedVals)))^2, na.rm = TRUE)
}

library(matrixcalc)
library(GNAR)

# Line-by-line test conditions
# vts = PMIData
# net = NNTradeNet
# alphaOrder = 2
# betaOrder = c(2, 1)
# fact.var = NULL
# globalalpha = FALSE
# xvts = stringencyData
# lambdaOrder = 1

# Comparisons with GNARfit
# alphaOrder = 1, betaOrder = 0, globalalpha = TRUE
# alphaOrder = 2, betaOrder = c(0, 0), globalalpha = TRUE
# alphaOrder = 2, betaOrder = c(2, 1), globalalpha = TRUE
# alphaOrder = 1, betaOrder = 0, globalalpha = FALSE
# alphaOrder = 2, betaOrder = c(0, 0), globalalpha = FALSE
# alphaOrder = 2, betaOrder = c(2, 1), globalalpha = FALSE
# BIC.GNARXfit(GNARXfit(vts = PMIData, net = NNTradeNet, alphaOrder = 2, betaOrder = c(2,1), 
#                       globalalpha = FALSE)) - BIC(GNARfit(vts = PMIData, 
#                                                                   net = NNTradeNet, 
#                                                                   alphaOrder = 2, 
#                                                           betaOrder = c(2,1),
#                                                                   globalalpha = FALSE))



# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

GNARXfit <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                      betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
                      tvnets = NULL, netsstart = NULL, ErrorIfNoNei = TRUE, 
                      xvts = NULL, lambdaOrder = NULL) 
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
    if(nrow(xvts) != nrow(vts)){
      stop("xvts and vts must have the same time range")
    }else{
      if(lambdaOrder%%1 != 0 || lambdaOrder < 0){
        stop("lambdaOrder must be a non-negative integer")
      }else{
        maxOrder <- max(alphaOrder, lambdaOrder)
      }
    }
  }
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder, 
                betas.in = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, 
                final.in = vts[(nrow(vts) - maxOrder + 1):nrow(vts), ], 
                lambdas.in = lambdaOrder)
  dmat <- GNARXdesign(vts = vts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart, 
                     xvts = xvts, lambdaOrder = lambdaOrder)
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 
                                            2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  predt <- nrow(vts) - alphaOrder
  yvec <- NULL
  for (ii in 1:length(net$edges)) {
    yvec <- c(yvec, vts[((alphaOrder + 1):(predt + alphaOrder)), 
                        ii])
  }
  if (sum(is.na(yvec)) > 0) {
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec), ]
    modNoIntercept <- lm(yvec2 ~ dmat2 + 0)
  }
  else {
    modNoIntercept <- lm(yvec ~ dmat + 0)
  }
  out <- list(mod = modNoIntercept, y = yvec, dd = dmat, frbic = frbic)
  class(out) <- "GNARXfit"
  return(out)
}


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

GNARXdesign <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                         betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
                         tvnets = NULL, netsstart = NULL, xvts = NULL, lambdaOrder = NULL) 
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
        maxOrder <- max(alphaOrder, lambdaOrder)
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
      parNames <- c(parNames, paste("lambda", ii, sep = ""))
      parLoc <- c(parLoc, "l") 
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
    for(jj in 0:lambdaOrder){
        dmat[, ncol(dmat) - lambdaOrder + jj] <- vec(xvts[(maxOrder + 1 - jj):
                                                            (nrow(xvts) - jj), ])
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


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

BIC.GNARXfit <- function (object, ...) 
{
  stopifnot(is.GNARfit(object))
  nnodes.in <- object$frbic$nnodes
  alphas.in <- object$frbic$alphas.in
  betas.in <- object$frbic$betas.in
  lambdas.in <- object$frbic$lambdas.in
  fact.var <- object$frbic$fact.var
  tot.time <- object$frbic$time.in
  globalalpha <- object$frbic$globalalpha
  dotarg <- list(...)
  if (length(dotarg) != 0) {
    if (!is.null(names(dotarg))) {
      warning("... not used here, input(s) ", paste(names(dotarg), 
                                                    collapse = ", "), " ignored")
    }
    else {
      warning("... not used here, input(s) ", paste(dotarg, 
                                                    collapse = ", "), " ignored")
    }
  }
  if (!is.null(fact.var)) {
    f.in <- length(unique(fact.var))
  } else {
    f.in <- 1
  }
  stopifnot(is.logical(globalalpha))
  stopifnot(length(nnodes.in) == 1)
  stopifnot(floor(nnodes.in) == nnodes.in)
  stopifnot(tot.time != 0)
  tmp.resid <- residToMat(GNARobj = object, nnodes = nnodes.in)$resid
  tmp.resid[is.na(tmp.resid)] <- 0
  larg <- det((1/tot.time) * t(tmp.resid) %*% tmp.resid)
  stopifnot(larg != 0)
  tmp1 <- log(larg)
  if(!is.null(lambdas.in)){
    if (globalalpha) {
      tmp2 <- f.in * (alphas.in + sum(betas.in) + lambdas.in + 1) * log(tot.time)/tot.time
    }else {
      tmp2 <- (ncol(tmp.resid) * alphas.in + sum(betas.in) + lambdas.in + 1) * 
        log(tot.time)/tot.time
    }
  }else{
    if (globalalpha) {
      tmp2 <- f.in * (alphas.in + sum(betas.in)) * log(tot.time)/tot.time
    }else {
      tmp2 <- (ncol(tmp.resid) * alphas.in + sum(betas.in)) * 
        log(tot.time)/tot.time
    }
  }
  return(tmp1 + tmp2)
}
  
  
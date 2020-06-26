# This set of functions uses the restriction matrix to compute estimates for the 
# static-network GNARX model. Specifically, it regresses vec(Y) on (Z^T\otimes I_N)R.

library(GNAR)
library(matrixcalc)

# Test: comparison with GNARfit function (i.e. no exogenous regressors)

testFitGNARX <- GNARXfit(vts = fiveVTS, xvts = NULL, net = fiveNet, alphaOrder = 2, 
                        betaOrder = c(2,1), lambdaOrder = NULL, globalalpha = FALSE) 
testFitGNAR <- GNARfit(vts = fiveVTS, net = fiveNet, alphaOrder = 2, 
                         betaOrder = c(2,1), globalalpha = FALSE) 
summary(testFitGNARX$mod)
summary(testFitGNAR$mod)
print(testFitGNARX$BIC)
print(BIC(testFitGNAR))

testFitGNARX.loc <- GNARXfit(vts = fiveVTS, xvts = NULL, net = fiveNet, alphaOrder = 2, 
                         betaOrder = c(2, 1), lambdaOrder = NULL, globalalpha = TRUE) 
testFitGNAR.loc <- GNARfit(vts = fiveVTS, net = fiveNet, alphaOrder = 2, 
                       betaOrder = c(2, 1), globalalpha = TRUE) 
summary(testFitGNARX.loc$mod)
summary(testFitGNAR.loc$mod)
print(testFitGNARX.loc$BIC)
print(BIC(testFitGNAR.loc))

testFunction1 <- function(){
  GNARXfit(vts = fiveVTS, xvts = NULL, net = fiveNet, alphaOrder = 2, 
           betaOrder = c(2, 1), lambdaOrder = NULL, globalalpha = FALSE)
}
testFunction2 <- function(){
  GNARfit(vts = fiveVTS, net = fiveNet, alphaOrder = 2, 
          betaOrder = c(2, 1), globalalpha = FALSE)
}  
benchmark(replications=10, testFunction1(), testFunction2(),
          columns=c('test', 'elapsed', 'replications'))

# -------------------------------------------------------------------------------------

GNARXfit <- function (vts = GNAR::fiveVTS, xvts = NULL, net = GNAR::fiveNet, alphaOrder = 2, 
                     betaOrder = c(1, 1), lambdaOrder = NULL, fact.var = NULL, globalalpha = FALSE, 
                     tvnets = NULL, netsstart = NULL, ErrorIfNoNei = TRUE) 
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
  if(is.null(lambdaOrder) & (!is.null(xvts))){
    cat("When external regressors included, a non-null lambdaOrder is required")
  }
  if(is.null(lambdaOrder)){
    stopifnot(is.null(xvts))
  }
  if(!is.null(xvts)){
    if(nrow(vts) != nrow(xvts)){
      stop("vts and xvts must have the same time dimension")
    }
  }
  useNofNei <- 1
  nnodes = length(net$edges)
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder,
                betas.in = betaOrder, lambda.in = lambdaOrder,  fact.var = fact.var, 
                globalalpha = globalalpha, xtsp = tsp(vts), time.in = nrow(vts), 
                net.in = net, final.in = vts[(nrow(vts) - alphaOrder + 1):nrow(vts), ])
  designList <- GNARXdesign(vts = vts, xvts = xvts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, lambdaOrder = lambdaOrder, fact.var = fact.var, 
                     globalalpha = globalalpha, tvnets = tvnets, netsstart = netsstart)
  
  dmat <- designList[[1]]
  Rmat <- designList[[2]]
  Zmat <- designList[[3]]
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  if(!is.null(lambdaOrder)){
    maxOrder <- max(alphaOrder, lambdaOrder)
  }else{
    maxOrder <- alphaOrder
  }
  predt <- nrow(vts) - maxOrder
  ymat <- t(vts)[,(maxOrder + 1):(predt + maxOrder)]
  yvec <- vec(ymat)
  if (sum(is.na(yvec)) > 0) {
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec), ]
    modNoIntercept <- lm(yvec2 ~ dmat2 + 0)
  }else{
    modNoIntercept <- lm(yvec ~ dmat + 0) 
  }
  M <- length(coef(modNoIntercept))
  vecBhat <- Rmat %*% coef(modNoIntercept)
  Bhat <- matrix(vecBhat, nrow = nnodes)
  Uhat <- ymat - Bhat %*% Zmat
  Uhat[is.na(Uhat)] <- 0
  modBIC <- log(det((1/(nrow(vts))) * Uhat %*% t(Uhat))) + (1/(nrow(vts))) *
    M * log(nrow(vts))
  # modBIC <- M*log(length(yvec[!is.na(yvec)])) - 2 * logLik(modNoIntercept)
  out <- list(mod = modNoIntercept, y = yvec, dd = dmat, frbic = frbic, BIC = modBIC)
  class(out) <- "GNARXfit"
  return(out)
}

# --------------------------------------------------------------------------------------------

GNARXdesign <- function (vts = GNAR::fiveVTS, xvts = NULL, net = GNAR::fiveNet, alphaOrder = 2, 
                        betaOrder = c(1, 1), lambdaOrder = NULL, fact.var = NULL, globalalpha = FALSE, 
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
  if(is.null(xvts)){
    lambdaOrder <- NULL
  }else{
    if(is.null(lambdaOrder)){
      stop("lambdaOrder must be non-null if xvts is non-null")
    }else{
      if(lambdaOrder%%1 != 0 || lambdaOrder < 0){
        stop("lambdaOrder must be a non-negative integer")
      }
    }
  }
  if(!is.null(xvts)){
    if(nrow(vts) != nrow(xvts)){
      stop("vts and xvts must have the same time dimension")
    }
  }
  parNames <- parLoc <- NULL
  for (jj in 1:alphaOrder) {                                  # For j in 1:p
    if (globalalpha) {                                        
      parNames <- c(parNames, paste("alpha", jj, sep = ""))
      parLoc <- c(parLoc, "a")
    }
    else {
      for (kk in 1:ncol(vts)) {                               # For i in 1:N
        parNames <- c(parNames, paste("alpha", jj, "node",    
                                      kk, sep = ""))
        parLoc <- c(parLoc, "a")                              
      }
    }
    if (betaOrder[jj] > 0) {                                  # If s_j>0
      for (kk in 1:betaOrder[jj]) {                           # For r in 1:s_j
        parNames <- c(parNames, paste("beta", jj, ".",        
                                      kk, sep = ""))          # E.g. parNames = c("alpha1node1", ..., "alpha1node5", "beta1.1", "alpha2node1", ..., "alpha2node5", "beta2.1")
        parLoc <- c(parLoc, "b")                              # E.g. parLoc = c("a", ..., "a", "b", "a", ..., "a", "b")
      }
    }
  }
  maxOrder <- alphaOrder
  if(!is.null(lambdaOrder)){
    for(ii in 0:lambdaOrder){
      parNames <- c(parNames, paste("lambda", ii, sep = ""))
      parLoc <- c(parLoc, "l") 
    }
    maxOrder <- max(alphaOrder, lambdaOrder)
  }
  predt <- nrow(vts) - maxOrder                             
  nnodes <- ncol(vts)
  Zmat <- t(vts[maxOrder:(nrow(vts)-1),])
  if(alphaOrder > 1){
    for(ii in 1:(alphaOrder-1)){
      Zmat <- rbind(Zmat, t(vts[(maxOrder-ii):(nrow(vts)-1-ii),]))
    }
  }
  if(!is.null(lambdaOrder)){
    for(jj in 1:(lambdaOrder + 1)){
      Zmat <- rbind(Zmat, t(xvts[(maxOrder+2-jj):(nrow(xvts)+1-jj),]))
    }
  }
  Lmat <- matrix(0, ncol=1, nrow=nnodes^2)
  ii2 <- 1
  for(ii in 1:nnodes){
    Lmat[ii2,1] <- 1
    ii2 <- ii2 + nnodes + 1
  }
  if(globalalpha == FALSE){
    Amat <- matrix(0, ncol=nnodes, nrow=1)
    Amat[1, 1] <- 1
    for(ii in 2:nnodes){
      Amat <- rbind(Amat, matrix(0, nnodes, nnodes))
      Amat2 <- matrix(0, ncol=nnodes, nrow=1)
      Amat2[1, ii] <- 1
      Amat <- rbind(Amat, Amat2)
    }
  }else{
    Amat <- Lmat
  }
  Rkmatlist <- list()
  for(ii in 1:alphaOrder){
    if(betaOrder[ii]!=0){
      Rkmat <- vec(as.matrix(net, stage = 1, normalise = TRUE))
      
      
      
      if(betaOrder[ii] > 1){
        for(jj in 2:betaOrder[ii]){
          Rkmat <- cbind(Rkmat, vec(as.matrix(net, stage = jj, normalise = TRUE)))
        }
      }
      Rkmatlist[[ii]] <- Rkmat
    }
  }
  if(!is.null(lambdaOrder)){
    Rmat <- matrix(0, nrow = (nnodes^2 * (alphaOrder + lambdaOrder + 1)), ncol = length(parLoc))
  }else{
    Rmat <- matrix(0, nrow = (nnodes^2 * alphaOrder), ncol = length(parLoc))
  }
  rowStart <- 1
  rowEnd <- nnodes^2
  colStart <- 1
  if(globalalpha == FALSE){
    colEnd <- nnodes + betaOrder[1]
    if(betaOrder[1] == 0){
      Rmat[rowStart:rowEnd, colStart:colEnd] <- Amat
    }else{
      Rmat[rowStart:rowEnd, colStart:colEnd] <- cbind(Amat, Rkmatlist[[1]])
    }
    if(alphaOrder > 1){
      for(ii in 2:alphaOrder){
        rowStart <- rowEnd + 1
        rowEnd <- rowEnd + nnodes^2
        colStart <- colEnd + 1
        colEnd <- colEnd + nnodes + betaOrder[ii]
        if(betaOrder[ii] == 0){
          Rmat[rowStart:rowEnd, colStart:colEnd] <- Amat
        }else{
          Rmat[rowStart:rowEnd, colStart:colEnd] <- cbind(Amat, Rkmatlist[[ii]])
        }
      }
    }
    if(!is.null(lambdaOrder)){
      rowStart <- rowEnd + 1
      rowEnd <- rowStart + nnodes^2 - 1
      colStart <- colEnd + 1
      Rmat[rowStart:rowEnd, colStart] <- Lmat
      if(lambdaOrder > 0){
        for(jj in 1:lambdaOrder){
          rowStart <- rowEnd + 1
          rowEnd <- rowStart + nnodes^2 - 1
          colStart <- colStart + 1
          Rmat[rowStart:rowEnd, colStart] <- Lmat
        }
      }
    }
  }else{
    colEnd <- 1 + betaOrder[1]
    if(betaOrder[1] == 0){
      Rmat[rowStart:rowEnd, colStart:colEnd] <- Amat
    }else{
      Rmat[rowStart:rowEnd, colStart:colEnd] <- cbind(Amat, Rkmatlist[[1]])
    }
    if(alphaOrder > 1){
      for(ii in 2:alphaOrder){
        rowStart <- rowEnd + 1
        rowEnd <- rowEnd + nnodes^2
        colStart <- colEnd + 1
        colEnd <- colEnd + 1 + betaOrder[ii]
        if(betaOrder[ii] == 0){
          Rmat[rowStart:rowEnd, colStart:colEnd] <- Amat
        }else{
          Rmat[rowStart:rowEnd, colStart:colEnd] <- cbind(Amat, Rkmatlist[[ii]])
        }
      }
    }
    if(!is.null(lambdaOrder)){
      rowStart <- rowEnd + 1
      rowEnd <- rowStart + nnodes^2 - 1
      colStart <- colEnd + 1
      Rmat[rowStart:rowEnd, colStart] <- Lmat
      if(lambdaOrder > 0){
        for(jj in 1:lambdaOrder){
          rowStart <- rowEnd + 1
          rowEnd <- rowStart + nnodes^2 - 1
          colStart <- colStart + 1
          Rmat[rowStart:rowEnd, colStart] <- Lmat
        }
      }
    }
  }
  dmat <- (t(Zmat) %x% diag(nnodes)) %*% Rmat
  colnames(dmat) <- parNames
  return(list(dmat, Rmat, Zmat))
}
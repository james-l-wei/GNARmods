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
  if (globalalpha == TRUE) {
    cat("Global alpha not yet supported")
  }
  stopifnot(globalalpha == FALSE)
  if(is.null(lambdaOrder) & (!is.null(xvts))){
    cat("When external regressors included, a non-null lambdaOrder is required")
  }
  if(is.null(lambdaOrder)){
    stopifnot(is.null(xvts))
  }
  useNofNei <- 1
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder,
                betas.in = betaOrder, lambda.in = lambdaOrder,  fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, final.in = vts[(nrow(vts) - 
                                                                                      alphaOrder + 1):nrow(vts), ])
  dmat <- GNARdesign(vts = vts, xvts = xvts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, lambdaOrder = lambdaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart)
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  predt <- nrow(vts) - alphaOrder
  yvec <- NULL
  ymat <- t(vts)
  for (ii in (alphaOrder + 1):(predt + alphaOrder)) {
    yvec <- c(yvec, ymat[, ii])
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
  if (globalalpha == TRUE) {
    cat("Global alpha not yet supported")
  }
  stopifnot(globalalpha == FALSE)
  if(is.null(lambdaOrder) & (!is.null(xvts))){
    cat("When external regressors included, a non-null lambdaOrder is required")
  }
  if(is.null(lambdaOrder)){
    stopifnot(is.null(xvts))
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
  if(!is.null(lambdaOrder)){
    for(ii in 0:lambdaOrder){
      parNames <- c(parNames, paste("lambda", ii, sep = ""))
      parLoc <- c(parLoc, "l") 
    }
  }
  maxOrder <- max(alphaOrder, lambdaOrder)                    # p*
  predt <- nrow(vts) - maxOrder                               # T-p*
  nnodes <- ncol(vts)                                         # N
  # if (globalalpha) {
  #   dmat <- matrix(0, nrow = predt * nnodes, 
  #                  ncol = sum(c(alphaOrder, betaOrder)), dimnames = list(NULL, parNames))
  # } else {
  #   dmat <- matrix(0, nrow = predt * nnodes, 
  #                  ncol = sum(c(alphaOrder * nnodes, betaOrder)), 
  #                  dimnames = list(NULL, parNames))  
  # }
  Zmat <- matrix(0, nrow = nnodes * (alphaOrder + lambdaOrder + 1), ncol = predt)
  Zmat <- t(vts[maxOrder:(nrow(vts)-1),])
  for(ii in 1:(alphaOrder-1)){
    Zmat <- rbind(Zmat, t(vts[(maxOrder-ii):(nrow(vts)-1-ii),]))
  }
  for(jj in 1:(lambdaOrder + 1)){
    Zmat <- rbind(Zmat, t(xvts[(maxOrder+2-jj):(nrow(xvts)+1-jj),]))
  }
  Amat <- matrix(0, ncol=nnodes, nrow=1)
  Amat[1, 1] <- 1
  for(ii in 2:nnodes){
    Amat <- rbind(Amat, matrix(0, nnodes, nnodes))
    Amat2 <- matrix(0, ncol=nnodes, nrow=1)
    Amat2[1, ii] <- 1
    Amat <- rbind(Amat, Amat2)
  }
  Lmat <- matrix(0, ncol=1, nrow=nnodes^2)
  ii2 <- 1
  for(ii in 1:nnodes){
    Lmat[ii2,1] <- 1
    ii2 <- ii2 + nnodes + 1
  }
  Rkmatlist <- list()
  for(ii in 1:alphaOrder){ # What about if betaorder=c(1,0)?
    if(betaOrder[[ii]]!=0){
      Rkmat <- vec(as.matrix(net, stage = 1))
      if(betaOrder[[ii]] > 1){
        for(jj in 2:betaOrder[[ii]]){
          Rkmat <- cbind(Rkmat, vec(as.matrix(net, stage = jj)))
        }
      }
      Rkmatlist[[ii]] <- Rkmat
    }
  }
  
  
  Rmat <- xxx
  dmat <- (t(Zmat) %x% diag(nnodes)) %*% Rmat
  
  
  
  for (ii in 1:nnodes) {                                     # For i in 1:N
    for (aa in 1:alphaOrder) {                               # For j in 1:p
      if (globalalpha) {
        alphaLoc <- which(parLoc == "a")[aa]
      }
      else {
        alphaLoc <- which(parLoc == "a")[nnodes * (aa - 1) + ii]
      }
      dmat[((predt * (ii - 1) + 1):(predt * ii)), alphaLoc] <- 
        vts[((maxOrder + 1 - aa):(predt + (maxOrder - aa))), ii] 
    }
  }
  if (sum(betaOrder) > 0) {
    betaN <- NULL
    betaTimes <- rep(1:alphaOrder, betaOrder)               # 1:p with each element repeated s_j times, e.g. c(1, 2)
    for (jj in 1:alphaOrder) {                              # For j in 1:p
      if (betaOrder[jj] > 0) {                              
        betaN <- c(betaN, 1:betaOrder[jj])                  # betaN is the vector of concatenated 1:s_j's, e.g. c(1, 1)
      }
    }
    for (ii in 1:nnodes) {                                  # For i in 1:N
      NofNei <- NofNeighbours(node = ii, stage = max(betaOrder),  # Neighbours of node i of stage max(s_j); [[1]] is vector of neighbouring nodes, [[2]] is vector of corresponding distances
                              net = net)
      Nei <- NofNei$edges                                   # Vector of edges connected to node i at stage max(s_j), e.g. for node 1, [[1]] c(4, 5)
      Wei <- NofNei$dist                                    # Vector of distances connected to node i at stage max(s_j), e.g. for node 1, [[1]] c(1,1)
      if ((!is.null(Nei)) & (length(Nei) > 0)) {            # If node i has neighbours
        if (!is.null(Nei[[1]])) {                           # If node i has neighbours
          Wei <- lapply(Wei, function(x) {
            1/(x * sum(1/x))
          })                                                # Note that Wei is a list, hence lapply; rescaled so weights sum to one
          for (bb in 1:sum(betaOrder)) {                    # For b in 1:\sum_j{s_j}
            betaLoc <- which(parLoc == "b")[bb]             # Position of the b-th beta parameter, e.g. for bb=1, betaLoc=6 
            if (length(Nei[[betaN[bb]]]) > 1) {             # E.g. length(c(4, 5))>1
              vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):(predt + 
                                                                (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]   # X_{ p+1-betaTimes[b] : T+(p-betaTimes[b]) , Nei[[beta[b]]]},
              # e.g. X_{2:199, c(4, 5)}; entries of the neighbouring nodes
              for (kk in 1:nrow(vts.cut)) {                 # For t' in T', where T' denotes time period of neighbour nodes  
                if (any(is.na(vts.cut[kk, ]))) {            # If there are any NAs at time t' in the neighbour nodes
                  if (all(is.na(vts.cut[kk, ]))) {          
                    vts.cut[kk, ] <- 0                      # If all neighbours have NAs at time t', change NA values to zeros 
                  }
                  else {
                    new.wei <- Wei[[betaN[bb]]][which(!is.na(vts.cut[kk,            # Remove weights corresponding to neighbours with NA values at time t'
                                                                     ]))]
                    new.wei <- new.wei/sum(new.wei)                                 # Weights rescaled so they sum to one
                    sub.val <- vts.cut[kk, which(!is.na(vts.cut[kk,                 # sub.val = X'_{t', which(not NA(X'_{t'}))} %*% new weights (i.e. total effect of non-NA neighbours)
                                                                ]))] %*% new.wei
                    vts.cut[kk, which(is.na(vts.cut[kk,                             # Replace NA values of neighbours by sub.val, so almost as if 'NA nodes' don't exist
                                                    ]))] <- sub.val
                  }
                }
              }
              dmat[((predt * (ii - 1) + 1):(predt * ii)), 
                   betaLoc] <- vts.cut %*% Wei[[betaN[bb]]]
            }
            else {
              if ((length(Nei[[betaN[bb]]]) == 1) & (!is.na(Nei[[betaN[bb]]]))) {
                vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):(predt + 
                                                                  (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]
                vts.cut[is.na(vts.cut)] <- 0
                dmat[((predt * (ii - 1) + 1):(predt * 
                                                ii)), betaLoc] <- vts.cut * Wei[[betaN[bb]]]
              }
              else {
                dmat[((predt * (ii - 1) + 1):(predt * 
                                                ii)), betaLoc] <- 0
              }
            }
          }
        }
        else {
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            dmat[((predt * (ii - 1) + 1):(predt * ii)), 
                 betaLoc] <- 0
          }
        }
      }
      else {
        for (bb in 1:sum(betaOrder)) {
          betaLoc <- which(parLoc == "b")[bb]
          dmat[((predt * (ii - 1) + 1):(predt * ii)), 
               betaLoc] <- 0
        }
      }
    }
  }
  if (is.null(fact.var)) {
    return(dmat)
  }
  else {
    facun <- unique(fact.var)
    if (length(facun) == 1) {
      return(dmat)
    }
    else {
      dmcol <- ncol(dmat)
      dmatex <- dmat
      exnames <- paste(colnames(dmat), " '", facun[1], 
                       "'", sep = "")
      for (ii in 2:length(facun)) {
        dmatex <- cbind(dmatex, dmat)
        exnames <- c(exnames, paste(colnames(dmat), " '", 
                                    facun[ii], "'", sep = ""))
      }
      for (ii in 1:length(facun)) {
        dmatex[fact.var != facun[ii], ((ii - 1) * dmcol + 
                                         (1:dmcol))] <- 0
      }
      colnames(dmatex) <- exnames
      return(dmatex)
    }
  }
}
# GNAR package source code ----------------------------------------------------------------------------------

GNARfit <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                     betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
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
  useNofNei <- 1
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder, 
                betas.in = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, final.in = vts[(nrow(vts) - 
                                                                                      alphaOrder + 1):nrow(vts), ])
  dmat <- GNARdesign(vts = vts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart)
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
  class(out) <- "GNARfit"
  return(out)
}





vts = GNAR::fiveVTS
net = GNAR::fiveNet
alphaOrder = 2 
betaOrder = c(1, 1)
fact.var = NULL
globalalpha = FALSE 
tvnets = NULL
netsstart = NULL

GNARdesign <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                        betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
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
  maxOrder <- alphaOrder                                      # p, e.g. 2
  predt <- nrow(vts) - maxOrder                               # T, e.g. 198
  nnodes <- ncol(vts)                                         # N, e.g. 5
  if (globalalpha) {
    dmat <- matrix(0, nrow = predt * nnodes, ncol = sum(c(alphaOrder, 
                                                          betaOrder)), dimnames = list(NULL, parNames))
  } else {
    dmat <- matrix(0, nrow = predt * nnodes, ncol = sum(c(alphaOrder * 
                                                            nnodes, betaOrder)), dimnames = list(NULL, parNames))  # NT\times M matrix of zeros
  }
  for (ii in 1:nnodes) {                                     # For i in 1:N
    for (aa in 1:alphaOrder) {                               # For j in 1:p
      if (globalalpha) {
        alphaLoc <- which(parLoc == "a")[aa]
      }
      else {
        alphaLoc <- which(parLoc == "a")[nnodes * (aa -           # E.g. which(parLoc == "a") = c(1, 2, ..., 11)
                                                     1) + ii]     # alphaloc is the N(j-1)+i entry of which(parLoc == "a"), e.g. for i=5 and j=2, 10th entry gives alphaLoc=11
      }
      dmat[((predt * (ii - 1) + 1):(predt * ii)), alphaLoc] <- vts[((maxOrder + 
                                                                       1 - aa):(predt + (maxOrder - aa))), ii] # Filled up alpha columns with values, D_{ T*(i-1)+1 : T*i , alphaloc} = 
      # X_{  (p+1-j) : T+(p-j), i }, T\times1 block assignment for each i and j
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
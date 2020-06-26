setwd("/home/james/covidProject")
PMIData <- as.matrix(read.csv("PMIData.csv")[,-1])
stringencyData <- as.matrix(read.csv("StringencyData.csv")[,-1])
tradeData <- read.csv("tradeData.csv")[,-1]
tradeNet <- matrixtoGNAR(tradeData)
NNTradeData <- matrix(0, dim(tradeData)[1], dim(tradeData)[2])
for(i in 1:nrow(NNTradeData)){
  NNTradeData[i, which.max(tradeData[i,])] <- 1
}
NNTradeNet <- matrixtoGNAR(NNTradeData)
plot(NNTradeNet, vertex.label = colnames(tradeData))

tradeData2 <- tradeData
for(i in 1:nrow(tradeData)){
  tradeData2[i,] <- tradeData[i,]/sum(tradeData[i,])
}

# alphaOrder* in 1:10
# betaOrder*[k] in 0:alphaOrder* 
# lambdaOrder in 1:10
# where * denotes optimum values based on BIC
# Note that we can't have beta orders above 1 for full trade network

# Test settings:
# test1: just PMIData; global alpha; full trade network
# test2: just PMIData; local alpha; full trade network
# test3: just PMIData; global alpha; nearest neighbour network
# test4: just PMIData; local alpha; nearest neighbour network
# test5: test1 with stringencyData
# test6: test2 with stringencyData
# test7: test3 with stringencyData
# test8: test4 with stringencyData

# Results:
# test1: BIC: 14014.5, alphaOrder: 1, betaOrder: 1
# test2: BIC: 12928.4, alphaOrder: 3, betaOrder: c(1,1,1)
# test3: BIC: 14417.0, alphaOrder: 1, betaOrder: 0
# test4: BIC: 13226.4, alphaOrder: 3, betaOrder: c(0,0,0)
# test5: BIC: 12000.3, alphaOrder: 5, betaOrder: c(1,1,1,0,1), lambdaOrder: 2
# test6: BIC: 11547.0, alphaOrder: 4, betaOrder: c(1,1,1,1), lambdaOrder: 2
# test7: BIC: 12124.9, alphaOrder: 5, betaOrder: c(1,1,1,0,0), lambdaOrder: 1
# test8: BIC: 11733.9, alphaOrder: 4, betaOrder: c(0,0,0,0), lambdaOrder: 2

test1 <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = TRUE, xvts = NULL, 
                    alphaLim = 4, lambdaLim = 4)
test2 <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = FALSE, xvts = NULL, 
                    alphaLim = 4, lambdaLim = 4)
test3 <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = TRUE, xvts = NULL, 
                    alphaLim = 4, lambdaLim = 4)
test4 <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = FALSE, xvts = NULL, 
                    alphaLim = 4, lambdaLim = 4)
test5 <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = TRUE, xvts = stringencyData, 
                    alphaLim = 4, lambdaLim = 4)
test6 <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = FALSE, xvts = stringencyData, 
                    alphaLim = 4, lambdaLim = 4)
test7 <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = TRUE, xvts = stringencyData, 
                    alphaLim = 4, lambdaLim = 4)
test8 <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = FALSE, xvts = stringencyData, 
                    alphaLim = 4, lambdaLim = 4)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------


# Line-by-line test
vts = PMIData
net = tradeNet
globalalpha = TRUE
xvts = NULL
alphaLim = 20
lambdaLim = 20

GNARXOrder <- function(vts, net, globalalpha, xvts, alphaLim, lambdaLim){
  bestAlpha <- 0
  if(is.null(xvts)){
    bestLambda <- NULL
  }else{
    bestLambda <- 0
  }
  newBIC <- BIC(GNARXfit(vts = vts, net = net, alphaOrder = 1, 
                         betaOrder = 0, globalalpha = globalalpha, xvts = xvts, 
                         lambdaOrder = bestLambda)$mod)
  bestBIC <- newBIC + 1
  while(newBIC < bestBIC & bestAlpha <= alphaLim){ # alpha search
    bestBIC <- newBIC
    bestAlpha <- bestAlpha + 1
    newBIC <- BIC(GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha + 1, 
                           betaOrder = rep(0, bestAlpha + 1), globalalpha = globalalpha, 
                           xvts = xvts, lambdaOrder = bestLambda)$mod)
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
      newBIC <- try(BIC(GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha, 
                             betaOrder = newBeta, globalalpha = globalalpha, 
                             xvts = xvts, lambdaOrder = bestLambda)$mod), silent = TRUE)
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
    while(newBIC < bestBIC & bestLambda <= lambdaLim){ # lambda search
      bestBIC <- newBIC
      bestLambda <- bestLambda + 1
      newBIC <- BIC(GNARXfit(vts = vts, net = net, alphaOrder = bestAlpha, 
                             betaOrder = bestBeta, globalalpha = globalalpha, 
                             xvts = xvts, lambdaOrder = bestLambda + 1)$mod)
    }
    bestOrder <- list(bestBIC = bestBIC, bestAlpha = bestAlpha, bestBeta = bestBeta, 
                      bestLambda = bestLambda)
  }
  return(bestOrder)
}


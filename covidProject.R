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

# alphaOrder* in 1:5
# betaOrder*[k] in 0:alphaOrder* 
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
# test1: BIC -6.326, alphaOrder 5, betaOrder (0, 0, 1, 1, 0)
# test2: BIC -7.069, alphaOrder 1, betaOrder 1
# test3: BIC -5.986, alphaOrder 5, betaOrder (0, 0, 0, 3, 0)
# test4: BIC -7.936, alphaOrder 1, betaOrder 4
# test5: 

# ------------------------------------------------------------------------------------

test1.alpha <- rep(NA, 5)
for(i in 1:length(test1.alpha)){
  test1.alpha[i] <- GNARXfit(vts = PMIDataNoNA, xvts = NULL, net = tradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), lambdaOrder = NULL, 
                             globalalpha = TRUE)$BIC
}
which.min(test1.alpha) # Best order is 5

test1.alpha <- rep(NA, 5)
for(i in 1:length(test1.alpha)){
  test1.alpha[i] <- BIC(GNARfit(vts = PMIDataNoNA, net = tradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), 
                             globalalpha = TRUE))
}
which.min(test1.alpha) # Best order is 5









test1.beta1 <- rep(NA, 6)
for(i in 0:length(test1.beta1)){
  test1.beta1[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 5, 
                             betaOrder = c(i,0,0,0,0), lambdaOrder = NULL, 
                             globalalpha = TRUE)$BIC
}
which.min(test1.beta1) - 1 # Best order is 0

test1.beta2 <- rep(NA, 6)
for(i in 0:length(test1.beta2)){
  test1.beta2[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,i,0,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test1.beta2) - 1 # Best order is 0

test1.beta3 <- rep(NA, 6)
for(i in 0:length(test1.beta3)){
  test1.beta3[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,i,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test1.beta3) - 1 # Best order is 1

test1.beta4 <- rep(NA, 6)
for(i in 0:length(test1.beta4)){
  test1.beta4[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,1,i,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test1.beta4) - 1 # Best order is 1

test1.beta5 <- rep(NA, 6)
for(i in 0:length(test1.beta5)){
  test1.beta5[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,1,1,i), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test1.beta5) - 1 # Best order is 0

# Hence test1: alphaOrder 5, betaOrder (0, 0, 1, 1, 0)

# ------------------------------------------------------------------------------------

test2.alpha <- rep(NA, 5)
for(i in 1:length(test2.alpha)){
  test2.alpha[i] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), lambdaOrder = NULL, 
                             globalalpha = FALSE)$BIC
}
which.min(test2.alpha) # Best order is 1

test2.beta1 <- rep(NA, 6)
for(i in 0:length(test2.beta1)){
  test2.beta1[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = tradeNet, alphaOrder = 1, 
                                 betaOrder = i, lambdaOrder = NULL, 
                                 globalalpha = FALSE)$BIC
}
which.min(test2.beta1) - 1 # Best order is 1

# ------------------------------------------------------------------------------------

test3.alpha <- rep(NA, 5)
for(i in 1:length(test3.alpha)){
  test3.alpha[i] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), lambdaOrder = NULL, 
                             globalalpha = TRUE)$BIC
}
which.min(test3.alpha) # Best order is 5

test3.beta1 <- rep(NA, 6)
for(i in 0:length(test3.beta1)){
  test3.beta1[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 5, 
                                 betaOrder = c(i,0,0,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test3.beta1) - 1 # Best order is 0

test3.beta2 <- rep(NA, 6)
for(i in 0:length(test3.beta2)){
  test3.beta2[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,i,0,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test3.beta2) - 1 # Best order is 0

test3.beta3 <- rep(NA, 6)
for(i in 0:length(test3.beta3)){
  test3.beta3[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,i,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test3.beta3) - 1 # Best order is 0

test3.beta4 <- rep(NA, 6)
for(i in 0:length(test3.beta4)){
  test3.beta4[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,0,i,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test3.beta4) - 1 # Best order is 3

test3.beta5 <- rep(NA, 6)
for(i in 0:length(test3.beta5)){
  test3.beta5[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,0,3,i), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test3.beta5) - 1 # Best order is 0

# Hence test3: alphaOrder 5, betaOrder (0, 0, 0, 3, 0)

# ------------------------------------------------------------------------------------

test4.alpha <- rep(NA, 5)
for(i in 1:length(test4.alpha)){
  test4.alpha[i] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), lambdaOrder = NULL, 
                             globalalpha = FALSE)$BIC
}
which.min(test4.alpha) # Best order is 1

test4.beta1 <- rep(NA, 6)
for(i in 0:length(test4.beta1)){
  test4.beta1[i + 1] <- GNARXfit(vts = PMIData, xvts = NULL, net = NNTradeNet, alphaOrder = 1, 
                                 betaOrder = i, lambdaOrder = NULL, 
                                 globalalpha = FALSE)$BIC
}
which.min(test4.beta1) - 1 # Best order is 4

# ------------------------------------------------------------------------------------

test5.alpha <- rep(NA, 5)
for(i in 1:length(test5.alpha)){
  test5.alpha[i] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = i, 
                             betaOrder = rep(0, i), lambdaOrder = 0, 
                             globalalpha = TRUE)$BIC
}
which.min(test5.alpha) # Best order is 5

vts = PMIData
xvts = stringencyData
net = tradeNet
alphaOrder = 1
betaOrder = 0
lambdaOrder = 0
globalalpha = TRUE

test5.beta1 <- rep(NA, 6)
for(i in 0:length(test5.beta1)){
  test5.beta1[i + 1] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(i,0,0,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test5.beta1) - 1 # Best order is 0

test5.beta2 <- rep(NA, 6)
for(i in 0:length(test5.beta2)){
  test5.beta2[i + 1] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,i,0,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test5.beta2) - 1 # Best order is 0

test5.beta3 <- rep(NA, 6)
for(i in 0:length(test5.beta3)){
  test5.beta3[i + 1] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,i,0,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test5.beta3) - 1 # Best order is 1

test5.beta4 <- rep(NA, 6)
for(i in 0:length(test5.beta4)){
  test5.beta4[i + 1] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,1,i,0), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test5.beta4) - 1 # Best order is 1

test5.beta5 <- rep(NA, 6)
for(i in 0:length(test5.beta5)){
  test5.beta5[i + 1] <- GNARXfit(vts = PMIData, xvts = stringencyData, net = tradeNet, alphaOrder = 5, 
                                 betaOrder = c(0,0,1,1,i), lambdaOrder = NULL, 
                                 globalalpha = TRUE)$BIC
}
which.min(test5.beta5) - 1 # Best order is 0

# Hence test5: alphaOrder 5, betaOrder (0, 0, 1, 1, 0)

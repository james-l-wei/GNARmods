setwd("/home/james/covidProject")
PMIData <- as.matrix(read.csv("PMIData.csv")[,-1]) - 50
stringencyData <- as.matrix(read.csv("StringencyData.csv")[,-1])
tradeData <- read.csv("tradeData.csv")[,-1]
tradeNet <- matrixtoGNAR(tradeData)
NNTradeData <- matrix(0, dim(tradeData)[1], dim(tradeData)[2])
for(i in 1:nrow(NNTradeData)){
  NNTradeData[i, which.max(tradeData[i,])] <- 1
}
NNTradeNet <- matrixtoGNAR(NNTradeData)
externalData <- rbind(rep(0, ncol(tradeData)), diff(stringencyData, lag=1))

# Without stringency indices

aGNAR1order <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 2, lambdaLim = 2, positiveCoef = FALSE)
aGNAR2order <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 2, lambdaLim = 2, positiveCoef = FALSE)
aGNAR3order <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 2, lambdaLim = 2, positiveCoef = FALSE)
aGNAR4order <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 2, lambdaLim = 2, positiveCoef = FALSE)

# --------------------------------------------------------------------------------------
# With stringency indices 

bGNAR1order <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = TRUE, 
                          xvts = externalData, alphaLim = 2, lambdaLim = 2, 
                          positiveCoef = FALSE)

bGNAR2order <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = TRUE, 
                          xvts = externalData, alphaLim = 2, lambdaLim = 2, 
                          positiveCoef = FALSE)

bGNAR3order <- GNARXOrder(vts = PMIData, net = tradeNet, globalalpha = FALSE, 
                          xvts = externalData, alphaLim = 2, lambdaLim = 2, 
                          positiveCoef = FALSE)
bGNAR4order <- GNARXOrder(vts = PMIData, net = NNTradeNet, globalalpha = FALSE, 
                          xvts = externalData, alphaLim = 2, lambdaLim = 2, 
                          positiveCoef = FALSE)

# --------------------------------------------------------------------------------------

# MSE calculations

aGNAR1MSE <- inSampleError(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1, 1),
                           globalalpha = TRUE, xvts = NULL, lambdaOrder = NULL, 
                           positiveCoef = FALSE)
aGNAR2MSE <- inSampleError(vts = PMIData, net = NNTradeNet, alphaOrder = 2, betaOrder = c(0, 0),
                           globalalpha = TRUE, xvts = NULL, lambdaOrder = NULL, 
                           positiveCoef = FALSE)
aGNAR3MSE <- inSampleError(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1, 1),
                           globalalpha = FALSE, xvts = NULL, lambdaOrder = NULL, 
                           positiveCoef = FALSE)
aGNAR4MSE <- inSampleError(vts = PMIData, net = NNTradeNet, alphaOrder = 2, betaOrder = c(0, 0),
                           globalalpha = FALSE, xvts = NULL, lambdaOrder = NULL, 
                           positiveCoef = FALSE)

bGNAR1MSE <- inSampleError(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1, 0),
                           globalalpha = TRUE, xvts = externalData, lambdaOrder = 2, 
                           positiveCoef = FALSE)
bGNAR2MSE <- inSampleError(vts = PMIData, net = NNTradeNet, alphaOrder = 2, betaOrder = c(0, 0),
                           globalalpha = TRUE, xvts = externalData, lambdaOrder = 2, 
                           positiveCoef = FALSE)
bGNAR3MSE <- inSampleError(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1, 0),
                           globalalpha = FALSE, xvts = externalData, lambdaOrder = 2, 
                           positiveCoef = FALSE)
bGNAR4MSE <- inSampleError(vts = PMIData, net = NNTradeNet, alphaOrder = 2, betaOrder = c(0, 0),
                           globalalpha = FALSE, xvts = externalData, lambdaOrder = 2, 
                           positiveCoef = FALSE)

bGNAR3fit <- GNARXfit(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1, 0),
                      globalalpha = TRUE, xvts = externalData, lambdaOrder = 0, 
                      positiveCoef = FALSE)
summary(bGNAR3fit$mod)

# --------------------------------------------------------------------------------------

# Fit charts

GNARmod <- GNARXfit(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1,1),
                globalalpha = TRUE, xvts = NULL, lambdaOrder = NULL,
                positiveCoef = FALSE)
GNARXmod <- GNARXfit(vts = PMIData, net = tradeNet, alphaOrder = 2, betaOrder = c(1,0),
                    globalalpha = TRUE, xvts = externalData, lambdaOrder = 0,
                    positiveCoef = FALSE)
GNARFittedVals <- GNARmod$dd %*% coef(GNARmod$mod)
GNARFittedVals <- matrix(GNARFittedVals, ncol = ncol(tradeData))
GNARXFittedVals <- GNARXmod$dd %*% coef(GNARXmod$mod)
GNARXFittedVals <- matrix(GNARXFittedVals, ncol = ncol(tradeData))


GNARFittedVals <- GNARFittedVals + 50
GNARXFittedVals <- GNARXFittedVals + 50
TrueSeries <- (PMIData + 50)[-(1:2),]
dateVector <- seq(from = as.Date("1998/1/31"), to = as.Date("2020/5/31"), 
                  by= "month")[-(1:2)]

countryNo <- 13
experiment2Frame1 <- data.frame(date = dateVector, trueSeries = TrueSeries[,countryNo], 
                                GNARSeries = GNARFittedVals[,countryNo], 
                                GNARXSeries = GNARXFittedVals[,countryNo])
experiment2Frame1 <- experiment2Frame1[complete.cases(experiment2Frame1),]
experiment2Frame1 <- tail(experiment2Frame1, 60)
experiment2Frame1.long <- melt(experiment2Frame1, id.var = "date")
ggplot(experiment2Frame1.long, aes(x = date, y = value, colour = variable)) + 
  geom_point() + 
  geom_line() +
  xlab(' ') +
  ylab('Composite PMI') + 
  theme(legend.position="bottom") +
  scale_color_discrete(name="", labels = c("Observations", "GNAR Fit", "GNARX Fit"))

# --------------------------------------------------------------------------------------

# Generate Forecasts

# First create new UK stringency series
# Then create code to use one-period forward estimates as next true observation of PMI

constantStringencyData <- stringencyData
for(i in 1:6){
  constantStringencyData <- rbind(constantStringencyData, tail(stringencyData, 1))
}
easingStringencyData <- constantStringencyData
easingDiff <- -(1/6)*easingStringencyData[269,2]
for(i in 1:6){
  easingStringencyData[269 + i,2] <- easingStringencyData[269,2] + i * easingDiff
}
tighteningStringencyData <- constantStringencyData
tighteningDiff <- (100 - tighteningStringencyData[269,2])/6
for(i in 1:6){
  tighteningStringencyData[269 + i,2] <- tighteningStringencyData[269,2] + i * tighteningDiff
}
constantExternalData <- rbind(rep(0, ncol(tradeData)), diff(constantStringencyData, lag=1))
easingExternalData <- rbind(rep(0, ncol(tradeData)), diff(easingStringencyData, lag=1))
tighteningExternalData <- rbind(rep(0, ncol(tradeData)), diff(tighteningStringencyData, lag=1))
  
constantPMIData <- PMIData
nObs <- nrow(PMIData)
nNodes <- ncol(PMIData)
for(i in 1:6){
  fit <- GNARXfit(vts = constantPMIData[1:(nObs + i - 1), ], net = tradeNet, 
                  alphaOrder = 2, betaOrder = c(1,0),
                  globalalpha = TRUE, xvts = constantExternalData[1:(nObs + i - 1), ], 
                  lambdaOrder = 0, positiveCoef = positiveCoef)
  newPMIPreds <- GNARXpredict(fit = fit, new.vts = rbind(constantPMIData[1:(nObs + i - 1), ], 
                                                         rep(NA, nNodes)), 
                              new.xvts = constantExternalData[1:(nObs + i), ]) 
  constantPMIData <- rbind(constantPMIData, newPMIPreds)
}

easingPMIData <- PMIData
nObs <- nrow(PMIData)
nNodes <- ncol(PMIData)
for(i in 1:6){
  fit <- GNARXfit(vts = easingPMIData[1:(nObs + i - 1), ], net = tradeNet, 
                  alphaOrder = 2, betaOrder = c(1,0),
                  globalalpha = TRUE, xvts = easingExternalData[1:(nObs + i - 1), ], 
                  lambdaOrder = 0, positiveCoef = positiveCoef)
  newPMIPreds <- GNARXpredict(fit = fit, new.vts = rbind(easingPMIData[1:(nObs + i - 1), ], 
                                                         rep(NA, nNodes)), 
                              new.xvts = easingExternalData[1:(nObs + i), ]) 
  easingPMIData <- rbind(easingPMIData, newPMIPreds)
}

tighteningPMIData <- PMIData
nObs <- nrow(PMIData)
nNodes <- ncol(PMIData)
for(i in 1:6){
  fit <- GNARXfit(vts = tighteningPMIData[1:(nObs + i - 1), ], net = tradeNet, 
                  alphaOrder = 2, betaOrder = c(1,0),
                  globalalpha = TRUE, xvts = tighteningExternalData[1:(nObs + i - 1), ], 
                  lambdaOrder = 0, positiveCoef = positiveCoef)
  newPMIPreds <- GNARXpredict(fit = fit, new.vts = rbind(tighteningPMIData[1:(nObs + i - 1), ], 
                                                         rep(NA, nNodes)), 
                              new.xvts = tighteningExternalData[1:(nObs + i), ]) 
  tighteningPMIData <- rbind(tighteningPMIData, newPMIPreds)
}


constantPMIVals <- constantPMIData + 50
constantPMIHi
constantPMILo

easingPMIVals <- easingPMIData + 50
tighteningPMIVals <- tighteningPMIData + 50
PMIVals <- constantPMIVals
PMIVals[270:275,] <- NA
dateVector <- seq(from = as.Date("1998/1/31"), to = as.Date("2020/12/31"), 
                  by= "month")[-1]
experiment2Frame2 <- data.frame(date = dateVector, constantPMIVals = constantPMIVals[,2], 
                                easingPMIVals = easingPMIVals[,2], 
                                tighteningPMIVals = tighteningPMIVals[,2],
                                PMIVals = PMIVals[,2])
experiment2Frame2 <- tail(experiment2Frame2, 60)
experiment2Frame2.long <- melt(experiment2Frame2, id.var = "date")
ggplot(experiment2Frame2.long, aes(x = date, y = value, colour = variable)) + 
  geom_point() + 
  geom_line() +
  xlab(' ') +
  ylab('Composite PMI') + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))+
  scale_color_discrete(name="", labels = c("Constant Stringency Path", "Easing Stringency Path", 
                                           "Tightening Stringency Path", "Observed Values"))


# 
# # Get the data from the web !
# CC <- read.table("http://www.sr.bham.ac.uk/~ajrs/papers/sanderson06/mean_Tprofile-CC.txt" ,  header=TRUE)
# nCC <- read.table("http://www.sr.bham.ac.uk/~ajrs/papers/sanderson06/mean_Tprofile-nCC.txt" , header=TRUE)
# CC$type <- "Cool core"
# nCC$type <- "Non-cool core"
# A <- rbind(CC, nCC)
# 
# 
# # Make the plot
# ggplot(data=A, aes(x=r.r500, y=sckT, ymin=sckT.lo, ymax=sckT.up, fill=type, linetype=type)) + 
#   geom_line() + 
#   geom_ribbon(alpha=0.5) + 
#   scale_x_log10() + 
#   scale_y_log10() + 
#   xlab(as.expression(expression( paste("Radius (", R[500], ")") ))) + 
#   ylab("Scaled Temperature")
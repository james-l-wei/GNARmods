library(MTS)
library(tseries)
library(ggplot2)
library(dplyr)
library(reshape2)

# Data sets

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

PMIDataNoAUS <- PMIData[,-12]
PMIDataNoAUSInSample <- PMIDataNoAUS[1:240,]

tradeDataNoAUS <- tradeData[,-12]
tradeDataNoAUS <- tradeDataNoAUS[-12,]
tradeNetNoAUS <- matrixtoGNAR(tradeDataNoAUS)

NNTradeDataNoAUS <- NNTradeData[,-12]
NNTradeDataNoAUS <- NNTradeDataNoAUS[-12,]
NNTradeNetNoAUS <- matrixtoGNAR(NNTradeDataNoAUS)

stringencyDataNoAUS <- stringencyData[,-12]

# Model settings:
# 1) net = tradeNetNoAUS, globalalpha = TRUE
# 2) net = NNTradeNetNoAUS, globalalpha = TRUE
# 3) net = tradeNetNoAUS, globalalpha = FALSE
# 4) net = NNTradeNetNoAUS, globalalpha = FALSE
# Model order selection methods:
# a) Greedy BIC-based method
# b) PACF based method (alphaOrder = 1, betarOrder from greed BIC-base method)
# c) Positive coefficients method

aGNAR1order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = TRUE, 
                         xvts = NULL, alphaLim = 12, lambdaLim = 3)
aGNAR2order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = TRUE, 
                         xvts = NULL, alphaLim = 12, lambdaLim = 3)
aGNAR3order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = FALSE, 
                         xvts = NULL, alphaLim = 12, lambdaLim = 3)
aGNAR4order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = FALSE, 
                         xvts = NULL, alphaLim = 12, lambdaLim = 3)

pacf(na.omit(PMIDataNoAUSInSample[,1]))
adf.test(na.omit(PMIDataNoAUSInSample[,1]))
bGNAR1order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 1, lambdaLim = 3)
bGNAR2order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 1, lambdaLim = 3)
bGNAR3order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 1, lambdaLim = 3)
bGNAR4order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 1, lambdaLim = 3)

cGNAR1order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 12, lambdaLim = 3, positiveCoef = TRUE)
cGNAR2order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = TRUE, 
                          xvts = NULL, alphaLim = 12, lambdaLim = 3, positiveCoef = TRUE)
cGNAR3order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 12, lambdaLim = 3, positiveCoef = TRUE)
cGNAR4order <- GNARXOrder(vts = PMIDataNoAUSInSample, net = NNTradeNetNoAUS, globalalpha = FALSE, 
                          xvts = NULL, alphaLim = 12, lambdaLim = 3, positiveCoef = TRUE)

VAR1order <- VARorder(x = PMIDataNoAUSInSample[complete.cases(PMIDataNoAUSInSample),], maxp = 5,
                      output = F)

aGNAR1MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                        new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                        globalalpha = TRUE, alphaOrder = 12, 
                        betaOrder = c(1,1,1,1,0,0,0,0,0,0,0,0), lambdaOrder = NULL, 
                        positiveCoef = FALSE)
aGNAR2MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net =NNTradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 12, 
                          betaOrder = c(0,0,0,0,0,0,0,0,0,0,0,0), lambdaOrder = NULL, 
                          positiveCoef = FALSE)
aGNAR3MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 1, 
                          betaOrder = 1, lambdaOrder = NULL, 
                          positiveCoef = FALSE) 
aGNAR4MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = NNTradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 1, 
                          betaOrder = 1, lambdaOrder = NULL, 
                          positiveCoef = FALSE) 

bGNAR1MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 1, 
                          betaOrder = 1, lambdaOrder = NULL, 
                          positiveCoef = FALSE)
bGNAR2MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net =NNTradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 1, 
                          betaOrder = 0, lambdaOrder = NULL, 
                          positiveCoef = FALSE)
bGNAR3MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 1, 
                          betaOrder = 1, lambdaOrder = NULL, 
                          positiveCoef = FALSE) 
bGNAR4MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = NNTradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 1, 
                          betaOrder = 1, lambdaOrder = NULL, 
                          positiveCoef = FALSE) 

cGNAR1MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 2, 
                          betaOrder = c(1,0), lambdaOrder = NULL, 
                          positiveCoef = TRUE)
cGNAR2MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net =NNTradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 2, 
                          betaOrder = c(1,0), lambdaOrder = NULL, 
                          positiveCoef = TRUE)
cGNAR3MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 2, 
                          betaOrder = c(1,0), lambdaOrder = NULL, 
                          positiveCoef = TRUE) 
cGNAR4MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = NNTradeNetNoAUS, 
                          globalalpha = FALSE, alphaOrder = 2, 
                          betaOrder = c(1,0), lambdaOrder = NULL, 
                          positiveCoef = TRUE)

VARMSFE <- VAR.oneStepMSFE(old.vts = PMIDataNoAUSInSample, new.vts = PMIDataNoAUS, order = 2)

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

GNARseries <- cGNAR1MSFE[[2]] + 50
VARseries <- VARMSFE[[2]] + 50
TrueSeries <- PMIDataNoAUS[(nrow(PMIDataNoAUSInSample)+1):nrow(PMIDataNoAUS),] + 50
dateVector <- seq(from = as.Date("2018/1/31"), to = as.Date("2020/5/31"), 
                  by= "month")

experiment1Frame1 <- data.frame(date = dateVector, trueSeries = TrueSeries[,12], 
                                GNARSeries = GNARseries[,12], VARSeries = VARseries[,12])
experiment1Frame1 <- experiment1Frame1[complete.cases(experiment1Frame1),]
experiment1Frame1.long <- melt(experiment1Frame1, id.var = "date")
ggplot(experiment1Frame1.long, aes(x = date, y = value, colour = variable)) + 
  geom_point() + 
  geom_line() +
  xlab(' ') +
  ylab('Composite PMI') + 
  theme(legend.position="bottom") +
  scale_color_discrete(name="", labels = c("Observations", "GNAR Forecast", "VAR Forecast"))


# Naive Forecast -------------------------------------------------

TrueSeriesLagged <- PMIDataNoAUS[(nrow(PMIDataNoAUSInSample)):(nrow(PMIDataNoAUS)-1),] + 50
mean((TrueSeries - TrueSeriesLagged)^2, na.rm = TRUE)

# Regression summary table -------------------------------------------------

GNARRegTable <- summary(GNARXfit(vts = PMIDataNoAUSInSample, net = tradeNetNoAUS, alphaOrder = 2, 
                                 betaOrder = c(1,0), positiveCoef = TRUE)$mod)

cGNAR1MSFE <- oneStepMSFE(old.vts = PMIDataNoAUSInSample, old.xvts = NULL, 
                          new.vts = PMIDataNoAUS, new.xvts = NULL, net = tradeNetNoAUS, 
                          globalalpha = TRUE, alphaOrder = 2, 
                          betaOrder = c(1,0), lambdaOrder = NULL, 
                          positiveCoef = TRUE)
    
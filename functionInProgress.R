library(casebase)
library(survival)
library(riskRegression)
library(reprex)
library(reshape2)
library(ggplot2)
library(prodlim)
library(data.table)


brierScoreKMCens <- function(predSurvs, times, newData) {
  # inverse probability censoring weights
  # Note that we take probabilities right before the drops.
  
  #DIFFERENCES START HERE
  fitCens <- prodlim::prodlim(Hist(time, event != 0) ~ 1, newData, reverse = TRUE)
  IPCW.subject.times <- prodlim::predictSurvIndividual(fitCens, lag = 1) # G(t-|X)
  #DIFFERENCES END HERE
  
  # Empty matrix that will be filled in with the following loop
  Score <- matrix(NA, nrow(predSurvs), ncol(predSurvs))
  IPCWMatrix <- matrix(NA, nrow(predSurvs), ncol(predSurvs))
  # for each point in time we have predicted
  for (i in 1:length(times)) {
    # get number of censored individuals so long as their survival time is less than times[i]
    # these individuals do not have an effect on right-censored brier score.
    CensBefore <- newData$event == 0 & newData$time < times[i]
    # y encompasses all survival times larger than t[i] with a 1, 0 otherwise
    y <- drop(t(newData$time > times[i]))
    # above permits the two parts of right-censored brier score to be calculated, without IPCW, in one line
    Score[i, ] <- (y - predSurvs[i, ])^2

    # Generate IPCWMatrix
    IPCWMatrix[i, y == 0] <- IPCW.subject.times[y == 0] # G(t-|X) filled in corresponding positions
    IPCW.time <- predict(fitCens, newdata = newData, times = times[i], level.chaos = 1, mode = "matrix", type = "surv", lag = 1)
    IPCWMatrix[i, y == 1] <- IPCW.time # G(t) filled, same value, for remaining positions.
    # above calculated individuals who don't have an effect on brier score at specfied time. they're scores are set to 0.
    Score[i, CensBefore] <- 0
  }

  # apply IPCW to all scores
  Err <- Score / IPCWMatrix

  # Average curve demonstrating right-censored brier averaged over test-set, for each time of interest.
  Err <- apply(Err, 1, mean)
  return(Err)
}

brierScoreweibCens <- function(predSurvs= absRiskcb, times=times, newDataX=newDataX,newDataY=newDataY) {
  # inverse probability censoring weights
  
  #DIFFERENCES START HERE
  yCens <- newDataY
  yCens[,1] <- yCens[,1]==0
  fitCens <- fitSmoothHazard(event~log(time) ,data= as.data.frame(yCens), time="time",ratio=100)
  IPCW.subject.times <- absoluteRisk(object = fitCens, time = times, s = "lambda.1se")
  IPCW.subject.times[, -c(1)] <- 1 - IPCW.subject.times[, -c(1)] # make it survival probabilities
  IPCW.subject.times <- rowMeans(IPCW.subject.times[-c(1), -c(1)]) # remove extra time at 0. Also remove the first column (times) as it matches perfectly now
  #DIFFERENCES END HERE
  
  # Empty matrix that will be filled in with the following loop
  Score <- matrix(NA, nrow(predSurvs), ncol(predSurvs))
  IPCWMatrix <- matrix(NA, nrow(predSurvs), ncol(predSurvs))
  # for each point in time we have predicted
  for (i in 1:length(times)) {
    # get number of censored individuals so long as their survival time is less than times[i]
    # these individuals do not have an effect on right-censored brier score.
    CensBefore <- newData$event == 0 & newData$time < times[i]
    # y encompasses all survival times larger than t[i] with a 1, 0 otherwise
    y <- drop(t(newData$time > times[i]))
    # above permits the two parts of right-censored brier score to be calculated, without IPCW, in one line
    Score[i, ] <- (y - predSurvs[i, ])^2
    
    # Generate IPCWMatrix
    IPCWMatrix[i, y == 0] <- IPCW.subject.times[y == 0] # G(t-|X) filled in corresponding positions
    IPCW.time <- IPCW.subject.times[i]
    IPCWMatrix[i, y == 1] <- IPCW.time # G(t) filled, same value, for remaining positions.
    # above calculated individuals who don't have an effect on brier score at specfied time. they're scores are set to 0.
    Score[i, CensBefore] <- 0
  }
  
  # apply IPCW to all scores
  Err <- Score / IPCWMatrix
  
  # Average curve demonstrating right-censored brier averaged over test-set, for each time of interest.
  Err <- apply(Err, 1, mean)
  return(Err)
}


#Simulate/prepare data
set.seed(18)

astrain <- simActiveSurveillance(278)
astrain$event <- astrain$event != 0 #take care of competing risks
data.table::setorder(astrain, time, -event)
astrainX <- model.matrix(event ~ . - time, data = astrain)[, -c(1)] #covariates
astrainY <- data.matrix(astrain[, c(9, 8)]) #response

newData <- simActiveSurveillance(208)
newData$event <- newData$event != 0 #take care of competing risks
data.table::setorder(newData, time, -event)
newDataX <- model.matrix(event ~ . - time, data = newData)[, -c(1)] #covariates
newDataY <- data.matrix(newData[, c(9, 8)]) #response

times <- sort(unique(newData$time))


# Fit a cox model
coxfit <- coxph(Surv(time, event != 0) ~ ., data = astrain, x = TRUE)
# Get Cox survivals
coxPred <- summary(survfit(coxfit, newdata = newData), times = times)$surv
# calculate the IPA/brier scores with riskRegression
X2 <- Score(list("PredictionModel" = coxfit), data = newData, formula = Surv(time, event != 0) ~ 1, summary = "ipa", se.fit = 0L, metrics = "brier", contrasts = FALSE, times = times)


#Fit Weibull casebase model
cbModel <- fitSmoothHazard.fit(astrainX, astrainY, family = "glmnet", time = "time", event = "event", alpha = 0,formula_time = ~log(time), ratio = 100, standardize = TRUE)
#we require a specific format for the probabilities
absRiskcb <- absoluteRisk(object = cbModel, newdata = newDataX, time = times, s = "lambda.1se")
absRiskcb[, -c(1)] <- 1 - absRiskcb[, -c(1)] # make it survival probabilities
absRiskcb <- absRiskcb[-c(1), -c(1)] # remove extra time at 0. Also remove the first column (times) as it matches perfectly now


#data is prepared to calculate brier scores
kmCensoredBrierScores <- brierScoreKMCens(predSurvs = absRiskcb, times = times, newData = newData)
weibullCensoredBrierScores<-brierScoreweibCens(predSurvs=absRiskcb, times, newDataX,newDataY)

# restructuring to make plotting easier
resultsCB<-rbind(rbind(data.frame(brierScore=kmCensoredBrierScores,times=times,model="KMCensoring")),rbind(data.frame(brierScore=weibullCensoredBrierScores,times=times,model="WeibullCensoring")),data.frame(brierScore = X2$Brier$score$Brier[X2$Brier$score$model == "PredictionModel"],times=times,model="riskRegressionCOX"))




#brierScorePlot
ggplot(data = resultsCB, mapping = aes(x = times, y = brierScore, col = model)) +
  geom_line() +
  geom_rug(sides = "b") +
  xlab("Times") +
  ylab("Brier-Score")+
  ggtitle("KMcensoring vs. WeibullCensoring")




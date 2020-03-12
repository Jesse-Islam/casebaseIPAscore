library(casebase)
# See example usage at http://sahirbhatnagar.com/casebase/

library(survival)

library(riskRegression)
#> riskRegression version 2019.11.03



library(reprex)
library(reshape2)
library(ggplot2)
library(prodlim)
library(data.table)


brierScore <- function(predSurvs,times,newData){
  #inverse probability censoring weights
  #note that we remove the last time, as this would give us 0. I shift all times by 1.
  #IPCW=summary(survfit(Surv(time, event==0)~1, newData),c(min(times),head(times,-1)))$surv
  fitCens=prodlim::prodlim(Hist(time, event!=0)~1,newData,reverse = TRUE)
  IPCW.subject.times=prodlim::predictSurvIndividual(fitCens,lag=1)#G(t-|X)
  
  #Empty matrix that will be filled in with the following loop
  Score <- matrix(NA, nrow(predSurvs), ncol(predSurvs))
  IPCWMatrix<-matrix(NA, nrow(predSurvs), ncol(predSurvs))
  #for each point in time we have predicted
  for (i in 1:length(times)) {
    #get number of censored individuals so long as their survival time is less than times[i]
    #these individuals do not have an effect on right-censored brier score.
    CensBefore <- newData$event==0 & newData$time < times[i]
    #y encompasses all survival times larger than t[i] with a 1, 0 otherwise
    y <- drop(t(newData$time > times[i]))
    #above permits the two parts of right-censored brier score to be calculated, without IPCW, in one line
    Score[i,] <- (y - predSurvs[i,])^2
    
    #Generate IPCWMatrix
    IPCWMatrix[i,y==0]<-IPCW.subject.times[y==0]#G(t-|X) filled in corresponding positions
    IPCW.time <- predict(fitCens,newdata=newData,times=times[i],level.chaos=1,mode="matrix",type="surv",lag=1)
    IPCWMatrix[i,y==1]<-IPCW.time#G(t) filled, same value, for remaining positions.
    # above calculated individuals who don't have an effect on brier score at specfied time. they're scores are set to 0.
    Score[i, CensBefore] <- 0
  }
  
  #apply IPCW to all scores
  Err<-Score/IPCWMatrix
  
  #Average curve demonstrating right-censored brier averaged over test-set, for each time of interest.
  Err<-apply(Err,1,mean)
  return(Err)
}

set.seed(18)
astrain <- simActiveSurveillance(278)
data.table::setorder(astrain,time,-event)

newData <- simActiveSurveillance(208)
data.table::setorder(newData,time,-event)
times=sort(unique(newData$time))

cbModel <- fitSmoothHazard(event ~ ., data = astrain, ratio = 100)
absRiskcb <- absoluteRisk(object = cbModel,newdata=newData, time = times)

#> 
#> Attaching package: 'data.table'
#> The following objects are masked from 'package:reshape2':
#> 
#>     dcast, melt

#set seed and simulate data

#fit a cox model
coxfit <- coxph(Surv(time,event!=0)~.,data=astrain,x=TRUE)

#specify prediction times of interest

# calculate the IPA/brier scores

X2 <- Score(list("PredictionModel"=coxfit),data=newData,formula=Surv(time,event!=0)~1,summary="ipa",se.fit=0L,metrics="brier",contrasts=FALSE,times=times)



#predictions using cox model
predSurvs=summary(survfit(coxfit,newdata=newData),times=times)$surv

#restructuring to make plotting easier
results=data.frame(riskRegression=X2$Brier$score$Brier[X2$Brier$score$model=="PredictionModel"],CustomBrier=Err,times=times)
results <- melt(results, id.vars="times")
#> Warning in melt(results, id.vars = "times"): The melt generic in data.table has
#> been passed a data.frame and will attempt to redirect to the relevant reshape2
#> method; please note that reshape2 is deprecated, and this redirection is now
#> deprecated as well. To continue using melt methods from reshape2 while both
#> libraries are attached, e.g. melt.list, you can prepend the namespace like
#> reshape2::melt(results). In the next version, this warning will become an error.

ggplot(data=results,mapping=aes(x=times,y=value,col=variable))+
  geom_line()+
  geom_rug(sides = 'b')+
  xlab("Times")+
  ylab("Brier-Score")
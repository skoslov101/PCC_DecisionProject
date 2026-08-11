### Behavioral plots for WoF and RecDec for manuscript
rm(list=ls())

library(doBy)
library(ggplot2)
library(reshape2)
require(gridExtra)
library(lme4)
library(lmerTest)
library(RColorBrewer)

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}


nanmean=function(x)(mean(x,na.rm=T))

### You'll need to indicate the folder/set the path to where the behavioral data is stored
behDir='<path2folder>/data/behavior'
dat=read.csv(sprintf('%s/behavior_valueBasedDM.csv',behDir),header=T)
head(dat)
str(dat)
dat$Subject=as.factor(dat$Subject)

wofdat2=summaryBy(wasRisky~Subject*trial*blockNum*riskyEV,data=dat)
names(wofdat2)=c('SubName','trial','block','riskyEV','wasRisky')
head(wofdat2)

#For each participant perform logistic regression to be able to predict the indecision point and response across risky-expected values for visualization.
subList=unique(wofdat2$SubName)
wofdat2$block=as.factor(wofdat2$block)
head(wofdat2)
decPlots=data.frame(subject=character(),riskyEV=numeric(),prop=numeric())
for(subI in 1:length(subList)){
  subValue=subList[subI]
  thisSub=wofdat2[wofdat2$SubName==subValue,]
  locVal=thisSub$locName[1]
  glmSub=glmer(wasRisky~riskyEV+(1|block),data=thisSub,family='binomial')
  predframe <- data.frame(riskyEV=seq(2, 27, length.out = 100))
  predframe$prop <- predict(glmSub,newdata=predframe,type="response",re.form=NA)
  predframe$subject=subValue
  decPlots=rbind(decPlots,predframe)
  rm(predframe, subValue, thisSub)
}

result <- do.call(rbind, lapply(split(decPlots, decPlots$subject), function(subdf) {
  subdf[which.min(abs(subdf$prop - 0.5)), c("riskyEV", "prop")]
}))
mean(result$riskyEV) #Average point of 50% risky decision across participants fit independently

#Fit the model for the reporting statistic
allMod <- glmer(wasRisky ~ riskyEV + (1 + riskyEV | SubName),
                data = wofdat2,
                family = binomial)
confint(allMod, method = "Wald")
summary(allMod) #Reporting stats
predframe <- data.frame(riskyEV=seq(2, 27, length.out = 100))
predframe$prop <- predict(allMod,newdata=predframe,type="response",re.form=NA)
summary(allMod)


wofChoicePlot1=ggplot() + theme_bw() + geom_line(data=decPlots,aes(x=riskyEV,y=prop,group=subject),linewidth=1.25,color='grey70') + theme(legend.position="") + scale_x_continuous(limits=c(0,30),breaks=c(0,5,10,15,20,25,30)) +geom_line(data=predframe,aes(x=riskyEV,y=prop),color='black',linewidth=2,linetype='solid')
wofChoicePlot1
wofBeh1=sprintf('%s/figure1d_base.pdf',behDir)
ggsave(wofBeh1,height=7,width=10,wofChoicePlot1)

wofChoicePlot2=ggplot() + geom_line(data=decPlots,aes(x=riskyEV,y=prop,group=subject),linewidth=1.25,color='grey60') + theme(legend.position="") + scale_x_continuous(limits=c(2,27),breaks=c(0,5,10,15,20,25)) +geom_line(data=predframe,aes(x=riskyEV,y=prop),color='black',linewidth=3,linetype='solid') + theme(
  panel.background = element_blank(),  # remove gray panel
  panel.grid       = element_blank(),  # remove grid lines
  panel.border = element_blank(),
  axis.line = element_line(),
  panel.grid.major = element_line(color='grey80')
)
wofChoicePlot2
wofBeh2=sprintf('%s/figure1d_base2.pdf',behDir)
ggsave(wofBeh2,height=7,width=10,wofChoicePlot2)

#################
#################
## RDD section
#################
#################
rm(list=ls())
behDir='<path2folder>/data/behavior'
rddDat=read.csv(sprintf('%s/behavior_memoryBasedDM.csv',behDir))
head(rddDat)
str(rddDat)

rddDat$Subject=as.factor(rddDat$Subject)
rdddat2=summaryBy(wasRisky~Subject*trial*block*numExp,data=rddDat)
head(rdddat2)

rddSubs=unique(rddDat$Subject)

rdddat3=summaryBy(wasRisky.mean~Subject*numExp,data=rdddat2,FUN=c(mean))
names(rdddat3)=c('SubName','expNum','wasRisky')
str(rdddat3)

rddAllMod <- lmer(wasRisky ~ expNum + (1  | SubName),
                  data = rdddat3)
summary(rddAllMod)
confint(rddAllMod, method = "Wald") 
rddExpSum=summaryBy(wasRisky~expNum,data=rdddat3,FUN=c(mean))
names(rddExpSum)=c('expNum','wasRisky')

rddPlot1=ggplot() + theme_bw() + geom_line(data=rdddat3,aes(x=expNum,y=wasRisky,group=SubName,color=SubName),color="grey70",lwd=1.25) + geom_line(data=rddExpSum,aes(x=expNum,y=wasRisky),color="black",lwd=2,linetype='solid') + scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))
rddPlot1
rddBeh1=sprintf('%s/figure1d_bottom_base.pdf',behDir)
ggsave(rddBeh1,height=7,width=10,rddPlot1)

##############
##############
# Delayed Recognition
##############
##############
keepList=c('rddSubs','behDir')
rm(list = setdiff(ls(), c("rddSubs", "behDir")))

itemDat=read.csv(sprintf('%s/behavior_delayedRecognition.csv',behDir))
names(itemDat)[names(itemDat) == 'nRDD'] <- 'nRiskyRDD'
head(itemDat)
str(itemDat)

itemDat2 = itemDat[itemDat$Subject %in% rddSubs,]
unique(itemDat2$Subject)


# Create hit indicator (1 if hit, 0 if miss)
oldDat=itemDat2[(itemDat2$totExp>0) & !is.na(itemDat2$totExp),]
hrArr <- summaryBy(acc2 ~ Subject + totExp, data=oldDat, FUN=c(mean, length))
hrArr$numHit=hrArr$acc2.mean*hrArr$acc2.length
colnames(hrArr)=c("Subject","totExp","hitRate","totTrials","numHits")
hrArr$adjHR=(hrArr$numHits+0.5)/(hrArr$totTrials+1)

hrArr2=hrArr[hrArr$totExp!=5,]
hrArr2$totExp2=NaN
#Reframe old-safe
hrArr2$totExp2=ifelse(hrArr2$totExp==15, 5,hrArr2$totExp)

#New data
newDat=itemDat2[itemDat2$wasRDD==0,]
faRate <- summaryBy(acc2 ~ Subject + totExp, data=newDat, FUN=c(mean,length)) #So inverse of this would be false alarms
names(faRate)=c("Subject","totExp","CR","totNews")
faRate$fas=1-faRate$CR
faRate$numFAs=faRate$fas*faRate$totNews
faRate$adjFA <- (faRate$numFAs + 0.5)/(faRate$totNews + 1)

dDat=merge(hrArr2,faRate,by='Subject')

dDat$zHit=qnorm(dDat$adjHR)
dDat$zFA=qnorm(dDat$adjFA)

dDat$dPrime=dDat$zHit-dDat$zFA

dDat2=summaryBy(dPrime~totExp2,data=dDat,FUN=c(mean))
names(dDat2)=c("totExp2","dPrime")

itemPlot1=ggplot() + theme_bw() + geom_line(data=dDat,aes(x=totExp2,y=dPrime,group=Subject),lwd=1.5,color='grey70') + scale_x_continuous(breaks=1:5,labels=c('1','2','3','4','>5')) + scale_y_continuous(limits=c(-1,3),breaks=c(-1,0,1,2,3)) + geom_line(data=dDat2,aes(x=totExp2,y=dPrime),color='black',lwd=3,linetype='solid')
itemPlot1
itemName1=sprintf('%s/itemRec_dPrimeXnumExp_acrossSubs_pccOnly_032526_v1.pdf',behDir)
ggsave(itemName1,height=7,width=10,itemPlot1)

#Stats for recognition memory across exposure number
lmPrime1=lmer(dPrime~totExp2+(1|Subject),data=dDat)
summary(lmPrime1)
confint(lmPrime1, method = "Wald") 

#I also want the overall d-Prime for all old stimuli
allHR <- summaryBy(acc2 ~ Subject , data=oldDat, FUN=c(mean, length))
allHR$numHit=allHR$acc2.mean*allHR$acc2.length
colnames(allHR)=c("Subject","hitRate","totTrials","numHits")
allHR$adjHR=(allHR$numHits+0.5)/(allHR$totTrials+1)

allMerge1=merge(allHR,faRate,by='Subject')
allMerge1$zHit=qnorm(allMerge1$adjHR)
allMerge1$zFA=qnorm(allMerge1$adjFA)

allMerge1$dPrime=allMerge1$zHit-allMerge1$zFA
t.test(allMerge1$dPrime,mu=0)

########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
###               Apply to timeseries data              ####
###                   MemBased - Choice                 ####
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data. Recording sites that did not pass quality assessment (e.g. frequent interictal spiking activity) have been removed already.
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/memoryDM/hippocampus'

hippDat <- fread(sprintf('%s/hippMacro_memoryBasedDM_choiceLocked.csv',baseDir), data.table = FALSE)
head(hippDat)

regRDD1=summaryBy(psc~time*elecID*blockValue*locName,data=hippDat,FUN=c(mean))
str(regRDD1)

regRDD2=summaryBy(psc.mean~time*elecID*locName,data=regRDD1,FUN=c(mean))
names(regRDD2)=c('time','elecID','locName','psc')
regRDD2$elecID=as.factor(regRDD2$elecID)
str(regRDD2)


#Select for choice window and compute stats
choicePeriod=regRDD2[(regRDD2$time>1700) & (regRDD2$time<2200),]
choiceSum2=summaryBy(psc~elecID*locName,data=choicePeriod,FUN=c(mean))
names(choiceSum2)=c('elecID','locName','psc')
head(choiceSum2)
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2$locName=as.factor(choiceSum2$locName)
str(choiceSum2)

lmC1=lmer(psc~locName+(1|Subject),data=choiceSum2)
summary(lmC1)
confint(lmC1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
lmCp=lmer(psc~1+(1|Subject),data=pDat)
summary(lmCp)
confint(lmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
lmCa=lmer(psc~1+(1|Subject),data=aDat)
summary(lmCa) 
confint(lmCa,method='Wald')

############
############
### Risk Modulation Hippocampus Choice - Memory Task
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippDat")
rm(list = setdiff(ls(), to_keep))

riskMem1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippDat,FUN=c(mean))
str(riskMem1)

#Smooth for plotting
smoothRisk=data.frame()
smoothRisk <- do.call(rbind, lapply(split(riskMem1,with(riskMem1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))

smoothRisk <- smoothRisk[order(smoothRisk$elecID, 
                               smoothRisk$blockValue, 
                               smoothRisk$wasRisky), ]
rownames(smoothRisk) <- NULL


#Now average across blocks and then probes
riskMem2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk,FUN=c(mean))
names(riskMem2)=c('time','elecID','wasRisky','locName','psc')
riskMem2$elecID=as.factor(riskMem2$elecID)
str(riskMem2)


head(riskMem2)
riskMem3=summaryBy(psc~time*locName*wasRisky,data=riskMem2,FUN=c(nanmean,se))
names(riskMem3)=c('time','locName','wasRisky','psc','se')
head(riskMem3)

riskMem3$wasRisky=as.factor(riskMem3$wasRisky)
riskMem3$locName=as.factor(riskMem3$locName)
str(riskMem3)

riskMem3b=riskMem3[riskMem3$time>1000 & riskMem3$time<=3000,]

rddRiskyPlot1=ggplot() + theme_bw() + geom_line(data=riskMem3b,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=riskMem3b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black')  + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + ggtitle('RDD - ChoiceOnset: Risky vs. safe x ROI') + theme(panel.grid.minor = element_blank())
rddRiskyPlot1
plotName=sprintf('%s/hippMacro_memoryDM_byRisk_choice.pdf',outDir)
ggsave(plotName,height=7,width=10,rddRiskyPlot1)

# Now summarize/collect data for risky/safe diff plots
choicePeriod=riskMem2[(riskMem2$time>1700) & (riskMem2$time<2200),]
choiceSum1=summaryBy(psc~elecID*locName*wasRisky,data=choicePeriod,FUN=c(mean))
choiceSum1$timeperiod='choice'

choiceSum2=dcast(data=choiceSum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(choiceSum2)=c('elecID','locName','timeperiod','safe','risky')
choiceSum2$diff=choiceSum2$risky-choiceSum2$safe
head(choiceSum2)

choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2$locName=as.factor(choiceSum2$locName)

cLm1=lmer(diff~locName+(1|Subject),data=choiceSum2)
summary(cLm1)
confint(cLm1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
lmCp=lmer(diff~1+(1|Subject),data=pDat)
summary(lmCp)
confint(lmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
lmCa=lmer(diff~1+(1|Subject),data=aDat)
summary(lmCa)
confint(lmCa,method='Wald')

############
##################
########################
## Delay
#######################
##################
############
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/memoryDM/hippocampus'

hippDat <- fread(sprintf('%s/hippMacro_memoryBasedDM_respLocked.csv',baseDir), data.table = FALSE)
head(hippDat)

regRDD1=summaryBy(psc~time*elecID*blockValue*locName,data=hippDat,FUN=c(mean))

regRDD2=summaryBy(psc.mean~time*elecID*locName,data=regRDD1,FUN=c(mean))
names(regRDD2)=c('time','elecID','locName','psc')
regRDD2$elecID=as.factor(regRDD2$elecID)

#Compare against baseline stats
choicePeriod=regRDD2[(regRDD2$time>1500) & (regRDD2$time<2000),]
choiceSum2=summaryBy(psc~elecID*locName,data=choicePeriod,FUN=c(mean))
names(choiceSum2)=c('elecID','locName','psc')
head(choiceSum2)
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2$locName=as.factor(choiceSum2$locName)
str(choiceSum2)

lmC1=lmer(psc~locName+(1|Subject),data=choiceSum2)
summary(lmC1) #
confint(lmC1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
lmCp=lmer(psc~1+(1|Subject),data=pDat)
summary(lmCp) # 
confint(lmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
lmCa=lmer(psc~1+(1|Subject),data=aDat)
summary(lmCa) #  
confint(lmCa,method='Wald')


############
############
### Risk Modulation Hippocampus Response - Memory Task
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippDat")
rm(list = setdiff(ls(), to_keep))

regRDD1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippDat,FUN=c(mean))
str(regRDD1)

#Smooth for plotting
smoothRisk=data.frame()
smoothRisk <- do.call(rbind, lapply(split(regRDD1,with(regRDD1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))

smoothRisk <- smoothRisk[order(smoothRisk$elecID, 
                               smoothRisk$blockValue, 
                               smoothRisk$wasRisky), ]
rownames(smoothRisk) <- NULL

regRDD2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk,FUN=c(mean))
names(regRDD2)=c('time','elecID','wasRisky','locName','psc')
regRDD2$elecID=as.factor(regRDD2$elecID)
str(regRDD2)


head(regRDD2)
rddRisk3=summaryBy(psc~time*locName*wasRisky,data=regRDD2,FUN=c(nanmean,se))
names(rddRisk3)=c('time','locName','wasRisky','psc','se')
head(rddRisk3)

rddRisk3$wasRisky=as.factor(rddRisk3$wasRisky)
rddRisk3$locName=as.factor(rddRisk3$locName)
str(rddRisk3)

rddRisk3$wasRisky=as.factor(rddRisk3$wasRisky)
rddRisk3$locName=as.factor(rddRisk3$locName)
str(rddRisk3)

rddRisk3b=rddRisk3[rddRisk3$time>500 & rddRisk3$time<=2500,]
rddRiskyPlot2=ggplot() + theme_bw() + geom_line(data=rddRisk3b,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=rddRisk3b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black')  + scale_y_continuous(limits=c(-4,8),breaks=c(-4,0,4,8)) + ggtitle('RDD - ChoiceOnset: Risky vs. safe x ROI') + theme(panel.grid.minor = element_blank())
rddRiskyPlot2
plotName=sprintf('%s/hippMacro_memoryDM_byRisk_response.pdf',outDir)
ggsave(plotName,height=7,width=10,rddRiskyPlot2)


#Response period risk modulation stats
delayPeriod=regRDD2[(regRDD2$time>1500) & (regRDD2$time<2000),]
delaySum1=summaryBy(psc~elecID*locName*wasRisky,data=delayPeriod,FUN=c(mean))
delaySum1$timeperiod='delay'

delaySum2=dcast(data=delaySum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(delaySum2)=c('elecID','locName','timeperiod','safe','risky')
delaySum2$diff=delaySum2$risky-delaySum2$safe
delaySum2$Subject=substr(delaySum2$elecID,1,6)
delaySum2$Subject=as.factor(delaySum2$Subject)
delaySum2$locName=as.factor(delaySum2$locName)

dLm1=lmer(diff~locName+(1|Subject),data=delaySum2)
summary(dLm1)
confint(dLm1,method='Wald')

pDatD=delaySum2[delaySum2$locName=='postHipp',]
lmDp=lmer(diff~1+(1|Subject),data=pDatD)
summary(lmDp) 
confint(lmDp,method='Wald')

aDatD=delaySum2[delaySum2$locName=='antHipp',]
lmDa=lmer(diff~1+(1|Subject),data=aDatD)
summary(lmDa)
confint(lmDa,method='Wald')

############
##################
########################
## Feedback
#######################
##################
############
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

#Load the data
outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/memoryDM/hippocampus'

hippDat <- fread(sprintf('%s/hippMacro_memoryBasedDM_feedLocked.csv',baseDir), data.table = FALSE)
head(hippDat)

regRDD1=summaryBy(psc~time*elecID*blockValue*locName,data=hippDat,FUN=c(mean))
regRDD2=summaryBy(psc.mean~time*elecID*locName,data=regRDD1,FUN=c(mean))
names(regRDD2)=c('time','elecID','locName','psc')
regRDD2$elecID=as.factor(regRDD2$elecID)

choicePeriod=regRDD2[(regRDD2$time>1500) & (regRDD2$time<2000),]
choiceSum2=summaryBy(psc~elecID*locName,data=choicePeriod,FUN=c(mean))
names(choiceSum2)=c('elecID','locName','psc')
head(choiceSum2)
choiceSum2$Subject=substr(choiceSum2$elecID,1,6)
choiceSum2$Subject=as.factor(choiceSum2$Subject)
choiceSum2$locName=as.factor(choiceSum2$locName)
str(choiceSum2)

lmC1=lmer(psc~locName+(1|Subject),data=choiceSum2)
summary(lmC1) #
confint(lmC1,method='Wald')

pDat=choiceSum2[choiceSum2$locName=='postHipp',]
lmCp=lmer(psc~1+(1|Subject),data=pDat)
summary(lmCp) # 
confint(lmCp,method='Wald')

aDat=choiceSum2[choiceSum2$locName=='antHipp',]
lmCa=lmer(psc~1+(1|Subject),data=aDat)
summary(lmCa) #  
confint(lmCa,method='Wald')

############
############
### Risk Modulation Hippocampus Feedback - Memory Task
############
############
to_keep <- c("outDir", "baseDir","nanmean","se","hippDat")
rm(list = setdiff(ls(), to_keep))

regRDD1=summaryBy(psc~time*elecID*blockValue*wasRisky*locName,data=hippDat,FUN=c(mean))
str(regRDD1)

#Smooth for plotting
smoothRisk=data.frame()
smoothRisk <- do.call(rbind, lapply(split(regRDD1,with(regRDD1, paste(elecID, blockValue, wasRisky, sep = "_"))),function(chunk) {
  chunk <- chunk[order(chunk$time), ]
  chunk$smoothDat <- ksmooth(chunk$time, chunk$psc.mean, kernel = "normal", bandwidth = 100)$y
  chunk
}))

smoothRisk <- smoothRisk[order(smoothRisk$elecID, 
                               smoothRisk$blockValue, 
                               smoothRisk$wasRisky), ]
rownames(smoothRisk) <- NULL


regRDD2=summaryBy(smoothDat~time*elecID*wasRisky*locName,data=smoothRisk,FUN=c(mean))
names(regRDD2)=c('time','elecID','wasRisky','locName','psc')
regRDD2$elecID=as.factor(regRDD2$elecID)
str(regRDD2)


head(regRDD2)
rddRisk3=summaryBy(psc~time*locName*wasRisky,data=regRDD2,FUN=c(nanmean,se))
names(rddRisk3)=c('time','locName','wasRisky','psc','se')
head(rddRisk3)

rddRisk3$wasRisky=as.factor(rddRisk3$wasRisky)
rddRisk3$locName=as.factor(rddRisk3$locName)
str(rddRisk3)

rddRisk3$wasRisky=as.factor(rddRisk3$wasRisky)
rddRisk3$locName=as.factor(rddRisk3$locName)
str(rddRisk3)

rddRisk3b=rddRisk3[rddRisk3$time>500 & rddRisk3$time<=2500,]

rddRiskyPlot3=ggplot() + theme_bw() + geom_line(data=rddRisk3b,aes(x=time,y=psc,color=locName,linetype=wasRisky),lwd=2) + geom_ribbon(data=rddRisk3b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName,linetype=wasRisky),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black')  + scale_y_continuous(limits=c(-4,8.1),breaks=c(-4,0,4,8)) + ggtitle('RDD - Feedback: Risky vs. safe x ROI') + theme(panel.grid.minor = element_blank())
rddRiskyPlot3
plotName=sprintf('%s/hippMacro_memoryDM_byRisk_feedback.pdf',outDir)
ggsave(plotName,height=7,width=10,rddRiskyPlot3)


#Stats for risk modulation post feedback
feedPeriod=regRDD2[(regRDD2$time>1500) & (regRDD2$time<2000),]
feedSum1=summaryBy(psc~elecID*locName*wasRisky,data=feedPeriod,FUN=c(mean))
feedSum1$timeperiod='feedback'

feedSum2=dcast(data=feedSum1,elecID+locName+timeperiod~wasRisky,value.var="psc.mean")
names(feedSum2)=c('elecID','locName','timeperiod','safe','risky')
feedSum2$diff=feedSum2$risky-feedSum2$safe
head(feedSum2)
feedSum2$Subject=substr(feedSum2$elecID,1,6)
feedSum2$Subject=as.factor(feedSum2$Subject)
feedSum2$locName=as.factor(feedSum2$locName)

fLm1=lmer(diff~locName+(1|Subject),data=feedSum2)
summary(fLm1)
confint(fLm1,method='Wald')

pDatF=feedSum2[feedSum2$locName=='postHipp',]
lmFp=lmer(diff~1+(1|Subject),data=pDatF)
summary(lmFp)
confint(lmFp,method='Wald')

aDatF=feedSum2[feedSum2$locName=='antHipp',]
lmFa=lmer(diff~1+(1|Subject),data=aDatF)
summary(lmFa)
confint(lmFa,method='Wald')

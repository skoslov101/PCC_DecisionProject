rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(data.table)

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/delayRec'

delayRecDat <- fread(sprintf('%s/delayRec_hippAndPCC_trialOnset_preproc.csv',baseDir), data.table = FALSE)
head(delayRecDat)

drPeriod=delayRecDat[(delayRecDat$time>1700) & (delayRecDat$time<=2200),]
drSum1=summaryBy(psc~elecID*locName,data=drPeriod,FUN=c(mean))
drSum1$Subject=substr(drSum1$elecID,1,6)
drSum1$Subject=as.factor(drSum1$Subject)
drSum1$locName=as.factor(drSum1$locName)
drSum1$locName <- relevel(drSum1$locName, ref = "dPCC")

lmO1=lmer(psc.mean~locName+(1|Subject),data=drSum1)
summary(lmO1)
confint(lmO1,method='Wald')

#Compare to baseline of 0
pDat=drSum1[drSum1$locName=='postHipp',]
lmOp=lmer(psc.mean~1+(1|Subject),data=pDat)
summary(lmOp)
confint(lmOp,method='Wald')

aDat=drSum1[drSum1$locName=='antHipp',]
lmOa=lmer(psc.mean~1+(1|Subject),data=aDat)
summary(lmOa)
confint(lmOa,method='Wald')

dDat=drSum1[drSum1$locName=='dPCC',]
lmOd=lmer(psc.mean~1+(1|Subject),data=dDat)
summary(lmOd)
confint(lmOd,method='Wald')

vDat=drSum1[drSum1$locName=='vPCC',]
lmOv=lmer(psc.mean~1+(1|Subject),data=vDat)
summary(lmOv)
confint(lmOv,method='Wald')

#Plot
delayRecDat
delayRecDat2=summaryBy(psc~time*locName,data=delayRecDat,FUN=c(mean,se))
names(delayRecDat2)=c('time','locName','psc','se')
head(delayRecDat2)
delayRecDat2b=delayRecDat2[(delayRecDat2$time>=1000) & (delayRecDat2$time<=3000),]

basePlotBoth=ggplot() + theme_bw() + geom_line(data=delayRecDat2b,aes(x=time,y=psc,color=locName)) + geom_ribbon(data=delayRecDat2b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-5,15),breaks=c(-5,0,5,10,15)) + ggtitle('DelayRec Onset: Baseline Resp') + theme(panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(color='black'),axis.ticks.length = unit(4, "pt")) + coord_cartesian(expand = FALSE)
basePlotBoth
plotName=sprintf('%s/delayRec_trialOnset_allMacro.pdf',outDir)
ggsave(plotName,height=7,width=10,basePlotBoth)

####################
##########
#### Confidence response period
##########
####################
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/delayRec'

confDat <- fread(sprintf('%s/delayRec_hippAndPCC_confindence_preproc.csv',baseDir), data.table = FALSE)
head(confDat)

confTrim=confDat[(confDat$time>1500) & (confDat$time<=2000),]
confTrimSum=summaryBy(psc~elecID*locName,data=confTrim,FUN=c(mean))
confTrimSum$Subject=substr(confTrimSum$elecID,1,6)
confTrimSum$Subject=as.factor(confTrimSum$Subject)
confTrimSum$locName=as.factor(confTrimSum$locName)

confTrimSum$locName=relevel(confTrimSum$locName,ref='dPCC')
lmO1c2=lmer(psc.mean~locName+(1|Subject),data=confTrimSum)
summary(lmO1c2)
confint(lmO1c2,method='Wald')

#Compare to baseline (0)
pDat=confTrimSum[confTrimSum$locName=='postHipp',]
lmOp=lmer(psc.mean~1+(1|Subject),data=pDat)
summary(lmOp) 
confint(lmOp,method='Wald')

aDat=confTrimSum[confTrimSum$locName=='antHipp',]
lmOa=lmer(psc.mean~1+(1|Subject),data=aDat)
summary(lmOa)
confint(lmOa,method='Wald')

dDat=confTrimSum[confTrimSum$locName=='dPCC',]
lmOd=lmer(psc.mean~1+(1|Subject),data=dDat)
summary(lmOd)
confint(lmOd,method='Wald')

vDat=confTrimSum[confTrimSum$locName=='vPCC',]
lmOv=lmer(psc.mean~1+(1|Subject),data=vDat)
summary(lmOv)
confint(lmOv,method='Wald')

#Plot
confDat2=summaryBy(psc~time*locName,data=confDat,FUN=c(nanmean,se))
names(confDat2)=c('time','locName','psc','se')
confDat2b=confDat2[confDat2$time>=500 & confDat2$time<=2500,]


confPlot=ggplot() + theme_bw() + geom_line(data=confDat2b,aes(x=time,y=psc,color=locName)) + geom_ribbon(data=confDat2b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName),alpha=.3) + geom_vline(xintercept=1500,linetype='dashed',color='black') + scale_y_continuous(limits=c(-5,15),breaks=c(-5,0,5,10,15)) + ggtitle('ItemRec Choice Onset: Baseline Resp') + theme(panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(color='black'),axis.ticks.length = unit(4, "pt")) + coord_cartesian(expand = FALSE)
confPlot
plotName=sprintf('%s/delayRec_confidence_allMacro.pdf',outDir)
ggsave(plotName,height=7,width=10,confPlot)

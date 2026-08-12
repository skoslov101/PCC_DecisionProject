rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(data.table)

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/encoding'

encOnsetDat <- fread(sprintf('%s/encoding_hippAndPCC_trialOnset_preproc.csv',baseDir), data.table = FALSE)
head(encOnsetDat)

#Concat to time after trial onset and run stats
onsetPeriod=encOnsetDat[(encOnsetDat$time>1200) & (encOnsetDat$time<=1700),]
onsetSum1=summaryBy(psc~elecID*locName,data=onsetPeriod,FUN=c(mean))
onsetSum1$Subject=substr(onsetSum1$elecID,1,6)
onsetSum1$Subject=as.factor(onsetSum1$Subject)
onsetSum1$locName=as.factor(onsetSum1$locName)
onsetSum1$locName <- relevel(onsetSum1$locName, ref = "dPCC")

lmO1=lmer(psc.mean~locName+(1|Subject),data=onsetSum1)
summary(lmO1)
confint(lmO1,method='Wald')

#Compare against baseline of 0
pDat=onsetSum1[onsetSum1$locName=='postHipp',]
lmOp=lmer(psc.mean~1+(1|Subject),data=pDat)
summary(lmOp)
confint(lmOp,method='Wald')

aDat=onsetSum1[onsetSum1$locName=='antHipp',]
lmOa=lmer(psc.mean~1+(1|Subject),data=aDat)
summary(lmOa)
confint(lmOa,method='Wald')

dDat=onsetSum1[onsetSum1$locName=='dPCC',]
lmOd=lmer(psc.mean~1+(1|Subject),data=dDat)
summary(lmOd) 
confint(lmOd,method='Wald')

vDat=onsetSum1[onsetSum1$locName=='vPCC',]
lmOv=lmer(psc.mean~1+(1|Subject),data=vDat)
summary(lmOv) 
confint(lmOv,method='Wald')

#Plot
encDat2=summaryBy(psc~time*locName,data=encOnsetDat,FUN=c(mean,se))
names(encDat2)=c('time','locName','psc','se')
str(encDat2)

encDat2b=encDat2[encDat2$time>=500,]
basePlotBoth=ggplot() + theme_bw() + geom_line(data=encDat2b,aes(x=time,y=psc,color=locName)) + geom_ribbon(data=encDat2b,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName),alpha=.3) + geom_vline(xintercept=1000,linetype='dashed',color='black') + scale_y_continuous(limits=c(-5,20),breaks=c(-5,0,5,10,15,20)) + ggtitle('ItemRec Choice Onset: Baseline Resp') + theme(panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(color='black'),axis.ticks.length = unit(4, "pt")) + coord_cartesian(expand = FALSE)
basePlotBoth
plotName=sprintf('%s/encoding_trialOnset_allMacro.pdf',outDir)
ggsave(plotName,height=7,width=10,basePlotBoth)



###############################
################
### By Target-Responses
################
###############################
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

outDir='<path2folder>/figures/Figure4'
baseDir='<path2folder>/data/encoding'

encTargDat <- fread(sprintf('%s/encoding_hippAndPCC_target_preproc.csv',baseDir), data.table = FALSE)
head(encTargDat)

#Subset to period around target responses
targetPeriod=encTargDat[(encTargDat$time>750) & (encTargDat$time<=1250),]
targetSum1=summaryBy(psc~elecID*locName,data=targetPeriod,FUN=c(mean))
targetSum1$Subject=substr(targetSum1$elecID,1,6)
targetSum1$Subject=as.factor(targetSum1$Subject)
targetSum1$locName=as.factor(targetSum1$locName)
targetSum1$locName <- relevel(targetSum1$locName, ref = "dPCC")

lmT1=lmer(psc.mean~locName+(1|Subject),data=targetSum1)
summary(lmT1)
confint(lmT1,method='Wald')

#Compare against baseline of 0
pDat=targetSum1[targetSum1$locName=='postHipp',]
lmFp=lmer(psc.mean~1+(1|Subject),data=pDat)
summary(lmFp)
confint(lmFp,method='Wald')

aDat=targetSum1[targetSum1$locName=='antHipp',]
lmFa=lmer(psc.mean~1+(1|Subject),data=aDat)
summary(lmFa)
confint(lmFa,method='Wald')

dDat=targetSum1[targetSum1$locName=='dPCC',]
lmFd=lmer(psc.mean~1+(1|Subject),data=dDat)
summary(lmFd)
confint(lmFd,method='Wald')

vDat=targetSum1[targetSum1$locName=='vPCC',]
lmFv=lmer(psc.mean~1+(1|Subject),data=vDat)
summary(lmFv)
confint(lmFv,method='Wald')

#Plot
targDat2=summaryBy(psc~time*locName,data=encTargDat,FUN=c(mean,se))
names(targDat2)=c('time','locName','psc','se')

targPlot=ggplot() + theme_bw() + geom_line(data=targDat2,aes(x=time,y=psc,color=locName)) + geom_ribbon(data=targDat2,aes(x=time,ymin=psc-se,ymax=psc+se,fill=locName),alpha=.3) + geom_vline(xintercept=1000,linetype='dashed',color='black') + scale_y_continuous(limits=c(-5,20),breaks=c(-5,0,5,10,15,20)) + ggtitle('ItemRec Choice Onset: Baseline Resp') + theme(panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(color='black'),axis.ticks.length = unit(4, "pt")) + coord_cartesian(expand = FALSE)
targPlot
plotName=sprintf('%s/encoding_targetResp_allMacro.pdf',outDir)
ggsave(plotName,height=7,width=10,targPlot)
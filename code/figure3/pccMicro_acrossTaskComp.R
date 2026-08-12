### Code for computing the congruent/incongruent stats
rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(lme4)
library(lmerTest)
library(MuMIn)

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

valChoice <- fread(sprintf('%s/pccMicro_valueBasedDM_choiceLocked_smooth_byBlock.csv',baseDir), data.table = FALSE)
head(valChoice)

valSum=summaryBy(smoothDat~time*wire*cluster,data=valChoice,FUN=c(mean))
names(valSum)=c('time','wire','cluster','fr')

memDir='<path2folder>/data/memoryDM'
memChoice <- fread(sprintf('%s/pccMicro_memBasedDM_choiceLocked_smooth_byBlock.csv',memDir), data.table = FALSE)
head(memChoice)

memSum=summaryBy(smoothDat~time*wire*cluster,data=memChoice,FUN=c(mean))
names(memSum)=c('time','wire','cluster','fr')

#Now concat for just the main choice period
memSumChoice=memSum[(memSum$time>1700) & (memSum$time<=2201),]
memSumChoice$task='mem'

#Truncate times to be just about the feedback period
valSumChoice=valSum[(valSum$time>1700) & (valSum$time<=2201),]
valSumChoice$task='val'

#Add in cluster info
valSumChoice2=summaryBy(fr~wire*cluster,data=valSumChoice,FUN=c(mean))
memSumChoice2=summaryBy(fr~wire*cluster,data=memSumChoice,FUN=c(mean))

shared_elec <- intersect(valSumChoice2$wire, memSumChoice2$wire)

valSumChoice2s=valSumChoice2[valSumChoice2$wire %in% shared_elec, ]
memSumChoice2s=memSumChoice2[memSumChoice2$wire %in% shared_elec, ]

mergeChoice=merge(valSumChoice2s,memSumChoice2s,by='wire')
names(mergeChoice)=c('wire','valClust','valFR','memClust','memFR')


mergeChoice$valType=NaN
mergeChoice$valType=ifelse(mergeChoice$memClust==1 & mergeChoice$valClust<3,'same',mergeChoice$valType)
mergeChoice$valType=ifelse(mergeChoice$memClust==2 & mergeChoice$valClust<3,'flip',mergeChoice$valType)
mergeChoice$valType=ifelse(mergeChoice$memClust==1 & mergeChoice$valClust==3,'flip',mergeChoice$valType)
mergeChoice$valType=ifelse(mergeChoice$memClust==2 & mergeChoice$valClust==3,'same',mergeChoice$valType)
mergeChoice$valType=ifelse(mergeChoice$memClust>=3,'other',mergeChoice$valType)


mergeChoice$valType=as.factor(mergeChoice$valType)
mergeChoice$Subject=substr(mergeChoice$wire,1,6)
mergeChoice$Subject=as.factor(mergeChoice$Subject)
mergeChoice$timeperiod='choice'

#Number of each type
table(mergeChoice$valType)

#Now compute stats for onset period
choiceSame=mergeChoice[mergeChoice$valType=='same',]
lmcs=lmer(memFR~valFR+(1|Subject),data=choiceSame)
summary(lmcs)
confint(lmcs,method='Wald')
r.squaredGLMM(lmcs)

choiceFlip=mergeChoice[mergeChoice$valType=='flip',]
lmcf=lmer(memFR~valFR+(1|Subject),data=choiceFlip)
summary(lmcf)
confint(lmcf,method='Wald')
r.squaredGLMM(lmcf)

choiceOther=mergeChoice[mergeChoice$valType=='other',]
lmco=lmer(memFR~valFR+(1|Subject),data=choiceOther)
summary(lmco)
confint(lmco,method='Wald')
r.squaredGLMM(lmco)


#### Response
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

valChoice <- fread(sprintf('%s/pccMicro_valueBasedDM_respLocked_smooth_byBlock.csv',baseDir), data.table = FALSE)
head(valChoice)

valSum=summaryBy(smoothDat~time*wire*cluster,data=valChoice,FUN=c(mean))
names(valSum)=c('time','wire','cluster','fr')

memDir='<path2folder>/data/memoryDM'
memChoice <- fread(sprintf('%s/pccMicro_memBasedDM_respLocked_smooth_byBlock.csv',memDir), data.table = FALSE)
head(memChoice)

memSum=summaryBy(smoothDat~time*wire*cluster,data=memChoice,FUN=c(mean))
names(memSum)=c('time','wire','cluster','fr')

#Now concat for just the main choice period
memSumDelay=memSum[(memSum$time>1500) & (memSum$time<=2001),]
memSumDelay$task='mem'

#Truncate times to be just about the feedback period
valSumDelay=valSum[(valSum$time>1500) & (valSum$time<=2001),]
valSumDelay$task='val'

#Add in cluster info
valSumDelay2=summaryBy(fr~wire*cluster,data=valSumDelay,FUN=c(mean))
memSumDelay2=summaryBy(fr~wire*cluster,data=memSumDelay,FUN=c(mean))

shared_elec <- intersect(valSumDelay2$wire, memSumDelay2$wire)

valSumDelay2s=valSumDelay2[valSumDelay2$wire %in% shared_elec, ]
memSumDelay2s=memSumDelay2[memSumDelay2$wire %in% shared_elec, ]

mergeDelay=merge(valSumDelay2s,memSumDelay2s,by='wire')
names(mergeDelay)=c('wire','valClust','valFR','memClust','memFR')


mergeDelay$valType=NaN
mergeDelay$valType=ifelse(mergeDelay$memClust==1 & mergeDelay$valClust<3,'same',mergeDelay$valType)
mergeDelay$valType=ifelse(mergeDelay$memClust==2 & mergeDelay$valClust<3,'flip',mergeDelay$valType)
mergeDelay$valType=ifelse(mergeDelay$memClust==1 & mergeDelay$valClust==3,'flip',mergeDelay$valType)
mergeDelay$valType=ifelse(mergeDelay$memClust==2 & mergeDelay$valClust==3,'same',mergeDelay$valType)
mergeDelay$valType=ifelse(mergeDelay$memClust>=3,'other',mergeDelay$valType)

mergeDelay$valType=as.factor(mergeDelay$valType)
mergeDelay$Subject=substr(mergeDelay$wire,1,6)
mergeDelay$Subject=as.factor(mergeDelay$Subject)
mergeDelay$timeperiod='Response'


respSame=mergeDelay[mergeDelay$valType=='same',]
lmrs=lmer(memFR~valFR+(1|Subject),data=respSame)
summary(lmrs)
confint(lmrs,method='Wald')
r.squaredGLMM(lmrs)

respFlip=mergeDelay[mergeDelay$valType=='flip',]
lmrf=lmer(memFR~valFR+(1|Subject),data=respFlip)
summary(lmrf)
confint(lmrf,method='Wald')
r.squaredGLMM(lmrf)

respOther=mergeDelay[mergeDelay$valType=='other',]
lmro=lmer(memFR~valFR+(1|Subject),data=respOther)
summary(lmro)
confint(lmro,method='Wald')
r.squaredGLMM(lmro)


#### Feedback
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

valChoice <- fread(sprintf('%s/pccMicro_valueBasedDM_feedbackLocked_smooth_byBlock.csv',baseDir), data.table = FALSE)
head(valChoice)

valSum=summaryBy(smoothDat~time*wire*cluster,data=valChoice,FUN=c(mean))
names(valSum)=c('time','wire','cluster','fr')

memDir='<path2folder>/data/memoryDM'
memChoice <- fread(sprintf('%s/pccMicro_memBasedDM_feedLocked_smooth_byBlock.csv',memDir), data.table = FALSE)
head(memChoice)

memSum=summaryBy(smoothDat~time*wire*cluster,data=memChoice,FUN=c(mean))
names(memSum)=c('time','wire','cluster','fr')

#Now concat for just the main choice period
memSumFeed=memSum[(memSum$time>1500) & (memSum$time<=2001),]
memSumFeed$task='mem'

#Truncate times to be just about the feedback period
valSumFeed=valSum[(valSum$time>1500) & (valSum$time<=2001),]
valSumFeed$task='val'

#Add in cluster info
valSumFeed2=summaryBy(fr~wire*cluster,data=valSumFeed,FUN=c(mean))
memSumFeed2=summaryBy(fr~wire*cluster,data=memSumFeed,FUN=c(mean))

shared_elec <- intersect(valSumFeed2$wire, memSumFeed2$wire)

valSumFeed2s=valSumFeed2[valSumFeed2$wire %in% shared_elec, ]
memSumFeed2s=memSumFeed2[memSumFeed2$wire %in% shared_elec, ]

mergeFeed=merge(valSumFeed2s,memSumFeed2s,by='wire')
names(mergeFeed)=c('wire','valClust','valFR','memClust','memFR')


mergeFeed$valType=NaN
mergeFeed$valType=ifelse(mergeFeed$memClust==1 & mergeFeed$valClust<3,'same',mergeFeed$valType)
mergeFeed$valType=ifelse(mergeFeed$memClust==2 & mergeFeed$valClust<3,'flip',mergeFeed$valType)
mergeFeed$valType=ifelse(mergeFeed$memClust==1 & mergeFeed$valClust==3,'flip',mergeFeed$valType)
mergeFeed$valType=ifelse(mergeFeed$memClust==2 & mergeFeed$valClust==3,'same',mergeFeed$valType)
mergeFeed$valType=ifelse(mergeFeed$memClust>=3,'other',mergeFeed$valType)

mergeFeed$valType=as.factor(mergeFeed$valType)
mergeFeed$Subject=substr(mergeFeed$wire,1,6)
mergeFeed$Subject=as.factor(mergeFeed$Subject)
mergeFeed$timeperiod='feedback'


feedSame=mergeFeed[mergeFeed$valType=='same',]
lmfs=lmer(memFR~valFR+(1|Subject),data=feedSame)
summary(lmfs)
confint(lmfs,method='Wald')
r.squaredGLMM(lmfs)

feedFlip=mergeFeed[mergeFeed$valType=='flip',]
lmff=lmer(memFR~valFR+(1|Subject),data=feedFlip)
summary(lmff)
confint(lmff,method='Wald')
r.squaredGLMM(lmff)

feedOther=mergeFeed[mergeFeed$valType=='other',]
lmfo=lmer(memFR~valFR+(1|Subject),data=feedOther)
summary(lmfo)
confint(lmfo,method='Wald')
r.squaredGLMM(lmfo)

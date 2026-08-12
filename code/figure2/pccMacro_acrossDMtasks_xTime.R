##### Script to compute stats for cross task PCC macro electrode correlations and to make Figure 2f.
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
###                         Choice                      ####
###                        Across Task                  ####
########### ########## ############ ########### ############
########### ########## ############ ########### ############
########### ########## ############ ########### ############
rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(zoo)
library(MuMIn)

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/'


#Load the Value-Based choice data
pccDat=read.csv(sprintf('%s/valueDM/pccMacro_valueBasedDM_onsetLocked.csv',baseDir),header=T)
head(pccDat)

#Load the Memory-Based choice data
pccMem=read.csv(sprintf('%s/memoryDM/pccMacro_memoryBasedDM_onsetLocked.csv',baseDir),header=T)
head(pccMem)

#Base resp merged
regMem=summaryBy(psc~time*elecID*blockValue*locName*dendClust,data=pccMem,FUN=c(mean))
regMemb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regMem,FUN=c(mean))
names(regMemb)=c('time','elecID','locName','clustVal','psc')

regWof=summaryBy(psc~time*elecID*blockVal*locName*dendClust,data=pccDat,FUN=c(mean))
regWofb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regWof,FUN=c(mean))
names(regWofb)=c('time','elecID','locName','clustVal','psc')

regWof2=regWofb[regWofb$time>1500 & regWofb$time<=2000,]
regRDD2=regMemb[regMemb$time>1500 & regMemb$time<=2000,]

regRDD3=summaryBy(psc~elecID*locName*clustVal,data=regRDD2,FUN=c(mean))
regWof3=summaryBy(psc~elecID*locName*clustVal,data=regWof2,FUN=c(mean))

merge3 <- merge(regWof3, regRDD3[,c('elecID','psc.mean')], by = c("elecID"))
names(merge3)=c('elecID','locName','clustVal','valPSC','memPSC')
merge3$Subject=substr(merge3$elecID,1,6)
merge3$Subject=as.factor(merge3$Subject)
str(merge3)
merge3$clustVal=as.factor(merge3$clustVal)
merge3=merge3[merge3$clustVal=='1' | merge3$clustVal=='2',]

lmAlld=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='1',])
summary(lmAlld)
confint(lmAlld,method='Wald')
r.squaredGLMM(lmAlld)

lmAllv=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='2',])
summary(lmAllv)
confint(lmAllv,method='Wald')
r.squaredGLMM(lmAllv)

################################
################################
######                    ######
###### Delay              ######
######                    ######
################################
################################
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/'

#Load the Value-Based response data
pccDat=read.csv(sprintf('%s/valueDM/pccMacro_valueBasedDM_responseLocked.csv',baseDir),header=T)
head(pccDat)


#Load the Memory-Based response data
pccMem=read.csv(sprintf('%s/memoryDM/pccMacro_memoryBasedDM_responseLocked.csv',baseDir),header=T)
head(pccMem)

#Base resp merged
regMem=summaryBy(psc~time*elecID*blockVal*locName*dendClust,data=pccMem,FUN=c(mean))
regMemb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regMem,FUN=c(mean))
names(regMemb)=c('time','elecID','locName','clustVal','psc')

regWof=summaryBy(psc~time*elecID*blockVal*locName*dendClust,data=pccDat,FUN=c(mean))
regWofb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regWof,FUN=c(mean))
names(regWofb)=c('time','elecID','locName','clustVal','psc')

regWof2=regWofb[regWofb$time>1500 & regWofb$time<=2000,]
regRDD2=regMemb[regMemb$time>1500 & regMemb$time<=2000,]

regRDD3=summaryBy(psc~elecID*locName*clustVal,data=regRDD2,FUN=c(mean))
regWof3=summaryBy(psc~elecID*locName*clustVal,data=regWof2,FUN=c(mean))

merge3 <- merge(regWof3, regRDD3[,c('elecID','psc.mean')], by = c("elecID"))
names(merge3)=c('elecID','locName','clustVal','valPSC','memPSC')
merge3$Subject=substr(merge3$elecID,1,6)
merge3$Subject=as.factor(merge3$Subject)
str(merge3)
merge3$clustVal=as.factor(merge3$clustVal)
merge3=merge3[merge3$clustVal=='1' | merge3$clustVal=='2',]

lmAlld=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='1',])
summary(lmAlld)
confint(lmAlld,method='Wald')
r.squaredGLMM(lmAlld)

lmAllv=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='2',])
summary(lmAllv)
confint(lmAllv,method='Wald')
r.squaredGLMM(lmAllv)


#Plot relationship
respPlot1=ggplot() + theme_bw() + geom_point(data=merge3,aes(x=memPSC,y=valPSC,color=clustVal),size=5,alpha=.5) + geom_smooth(data=merge3,aes(x=memPSC,y=valPSC,color=clustVal),se=F,method='glm') + scale_y_continuous(limits=c(-10,20),breaks=c(-10,-5,0,5,10,15,20)) + scale_x_continuous(limits=c(-10,15),breaks=c(-10,-5,0,5,10,15)) + theme(panel.grid.minor = element_blank()) + theme(legend.position='') + xlab('Memory Amplitude (%d)') + ylab('Value Amplitude (%d)') + ggtitle('Post Response Task Correlation')
respPlot1
corrPlotName=sprintf('%s/macroPCC_memXval_resp_xClust.pdf',outDir)
ggsave(corrPlotName,height=7,width=10,respPlot1)


################################
################################
######                    ######
###### Feedback           ######
######                    ######
################################
################################
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}

outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/'

#Load the Value-Based feedback data
pccDat=read.csv(sprintf('%s/valueDM/pccMacro_valueBasedDM_feedbackLocked.csv',baseDir),header=T)
head(pccDat)

#Load the Memory-Based feedback data
pccMem=read.csv(sprintf('%s/memoryDM/pccMacro_memoryBasedDM_feedbackLocked.csv',baseDir),header=T)
head(pccMem)

#Base resp merged
regMem=summaryBy(psc~time*elecID*blockVal*locName*dendClust,data=pccMem,FUN=c(mean))
regMemb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regMem,FUN=c(mean))
names(regMemb)=c('time','elecID','locName','clustVal','psc')

regWof=summaryBy(psc~time*elecID*blockVal*locName*dendClust,data=pccDat,FUN=c(mean))
regWofb=summaryBy(psc.mean~time*elecID*locName*dendClust,data=regWof,FUN=c(mean))
names(regWofb)=c('time','elecID','locName','clustVal','psc')

regWof2=regWofb[regWofb$time>1500 & regWofb$time<=2000,]
regRDD2=regMemb[regMemb$time>1500 & regMemb$time<=2000,]

regRDD3=summaryBy(psc~elecID*locName*clustVal,data=regRDD2,FUN=c(mean))
regWof3=summaryBy(psc~elecID*locName*clustVal,data=regWof2,FUN=c(mean))

merge3 <- merge(regWof3, regRDD3[,c('elecID','psc.mean')], by = c("elecID"))
names(merge3)=c('elecID','locName','clustVal','valPSC','memPSC')
merge3$Subject=substr(merge3$elecID,1,6)
merge3$Subject=as.factor(merge3$Subject)
str(merge3)
merge3$clustVal=as.factor(merge3$clustVal)
merge3=merge3[merge3$clustVal=='1' | merge3$clustVal=='2',]

lmAlld=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='1',])
summary(lmAlld)
confint(lmAlld,method='Wald')
r.squaredGLMM(lmAlld)

lmAllv=lmer(valPSC~memPSC+(1|Subject),data=merge3[merge3$clustVal=='2',])
summary(lmAllv)
confint(lmAllv,method='Wald')
r.squaredGLMM(lmAllv)
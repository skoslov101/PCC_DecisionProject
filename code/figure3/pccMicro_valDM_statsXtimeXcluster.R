# ── Initiate Script ───────────────────────────────────────────────
rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(data.table)
library(lme4)
library(lmerTest)

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

choiceMicro <- fread(sprintf('%s/pccMicro_valueBasedDM_choiceLocked.csv',baseDir), data.table = FALSE)
head(choiceMicro)

#Collapse across trials for each block for each wire
choiceProbe=summaryBy(firingRate~time*wireVal2*Block*hcCluster,data=choiceMicro,FUN=c(nanmean))
names(choiceProbe)=c('time','wire','block','cluster','fr')
head(choiceProbe)

#Smooth
str(choiceProbe)
choiceProbe$time=as.integer(choiceProbe$time)
choiceProbe$wire=as.factor(choiceProbe$wire)
wireList=unique(choiceProbe$wire)
smoothDat=data.frame()
for(clustI in wireList){
  clustDatC=choiceProbe[choiceProbe$wire==clustI,]
  blockList=unique(clustDatC$block)
  for(blockI in 1:length(blockList)){
    blockVal=blockList[blockI]
    blockDat=clustDatC[clustDatC$block==blockVal,]
    blockDat$smoothDat <- ksmooth(blockDat$time, blockDat$fr, kernel = "normal", bandwidth = 50)$y
    smoothDat=rbind(smoothDat,blockDat)
  }
}

head(smoothDat)


#Stats for response vs. baseline for unit types during onset period
choiceProbe2b=summaryBy(smoothDat~time*wire*cluster,data=smoothDat,FUN=c(mean))
names(choiceProbe2b)=c('time','wire','cluster','fr')

choicePeriod=choiceProbe2b[(choiceProbe2b$time>1700) & (choiceProbe2b$time<2200),]
cpSum1=summaryBy(fr~wire*cluster,data=choicePeriod,FUN=c(mean))

str(cpSum1)
cpSum1$Subject=substr(cpSum1$wire,1,6)
cpSum1$Subject=as.factor(cpSum1$Subject)
cpSum1$cluster=as.factor(cpSum1$cluster)
lmerChoice1=lmer(fr.mean~cluster+(1|Subject),data=cpSum1)
summary(lmerChoice1)

cpSum1$cluster=relevel(cpSum1$cluster,ref='2')
lmerChoice1=lmer(fr.mean~cluster+(1|Subject),data=cpSum1)
summary(lmerChoice1)
confint(lmerChoice1,method='Wald')

cpSum11=cpSum1[cpSum1$cluster=='1',]
c1choice=lmer(fr.mean~1+(1|Subject),data=cpSum11)
summary(c1choice)
confint(c1choice,method='Wald')

cpSum12=cpSum1[cpSum1$cluster=='2',]
c2choice=lmer(fr.mean~1+(1|Subject),data=cpSum12)
summary(c2choice)
confint(c2choice,method='Wald')

cpSum13=cpSum1[cpSum1$cluster=='3',]
c3choice=lmer(fr.mean~1+(1|Subject),data=cpSum13)
summary(c3choice)
confint(c3choice,method='Wald')


#Plot for each cluster type
choiceSum2=summaryBy(fr~time*cluster,data=choiceProbe2b,FUN=c(mean,se,length))
head(choiceSum2)
str(choiceSum2)
names(choiceSum2)=c('time','cluster','zFR','se','len')
choiceSum2$cluster=as.factor(choiceSum2$cluster)

choiceSumTrim=choiceSum2[choiceSum2$time>=1000 & choiceSum2$time<=3000,]
clusPlot1=ggplot() + theme_bw() + geom_line(data=choiceSumTrim,aes(x=time,y=zFR,color=cluster)) + geom_ribbon(data=choiceSumTrim,aes(x=time,ymin=zFR-se,ymax=zFR+se,fill=cluster),alpha=.3) + geom_vline(xintercept=1500,color='black',linetype='solid') + scale_y_continuous(limits=c(-2,3),breaks=c(-2,-1,0,1,2,3)) + theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())
clusPlot1
plotName=sprintf('%s/pccMicro_valDM_choice_byClust.pdf',outDir)
ggsave(plotName,height=7,width=10,clusPlot1)

# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  Clusters x wasRisky
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────

#load the data
cname=sprintf('%s/pccMicro_valueBasedDM_byRisk_choiceLocked_smooth_byWire.csv',baseDir)
choiceSum2 <- fread(cname, data.table = FALSE)

rsChoice=choiceSum2[(choiceSum2$time>1700) & (choiceSum2$time<2201),]
rsChoice2=summaryBy(zFR~wasRisky*wire*cluster,data=rsChoice)
head(rsChoice2)

rsWide <- reshape2::dcast(
  rsChoice2,
  wire + cluster ~ wasRisky,
  value.var = "zFR.mean"
)
names(rsWide)=c('wire','cluster','safe','risky')

rsWide$diff=rsWide$risky-rsWide$safe
rsWide$cluster=as.factor(rsWide$cluster)
head(rsWide)

rsWide$Subject=substr(rsWide$wire,1,6)
rsWide$Subject=as.factor(rsWide$Subject)
str(rsWide)
rsWide$cluster=relevel(rsWide$cluster,ref='2')

cDat1=rsWide[rsWide$cluster=='1',]
cDat2=rsWide[rsWide$cluster=='2',]
cDat3=rsWide[rsWide$cluster=='3',]

c1lm1=lmer(diff~1+(1|Subject),data=cDat1)
summary(c1lm1)
confint(c1lm1,method='Wald')

c2lm1=lmer(diff~1+(1|Subject),data=cDat2)
summary(c2lm1)
confint(c2lm1,method='Wald')

c3lm1=lmer(diff~1+(1|Subject),data=cDat3)
summary(c3lm1) 
confint(c3lm1,method='Wald')


# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  Delay - Delay - Delay - Delay
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Load the by-wire data ───────────────────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

datName=sprintf('%s/pccMicro_valueBasedDM_respLocked_smooth_byBlock.csv',baseDir)
delaySumDat <- fread(datName, data.table = FALSE)

delayProbe=summaryBy(smoothDat~time*wire*cluster,data=delaySumDat,FUN=c(mean))
names(delayProbe)=c('time','wire','cluster','fr')

#Subset to the apriori time window
delayPeriod=delayProbe[(delayProbe$time>1500) & (delayProbe$time<2000),]
dpSum1=summaryBy(fr~wire*cluster,data=delayPeriod,FUN=c(mean))

str(dpSum1)
dpSum1$Subject=substr(dpSum1$wire,1,6)
dpSum1$Subject=as.factor(dpSum1$Subject)
dpSum1$cluster=as.factor(dpSum1$cluster)

dpSum11=dpSum1[dpSum1$cluster=='1',]
d1choice=lmer(fr.mean~1+(1|Subject),data=dpSum11)
summary(d1choice)
confint(d1choice,method='Wald')

dpSum12=dpSum1[dpSum1$cluster=='2',]
d2choice=lmer(fr.mean~1+(1|Subject),data=dpSum12)
summary(d2choice)
confint(d2choice,method='Wald')

dpSum13=dpSum1[dpSum1$cluster=='3',]
d3choice=lmer(fr.mean~1+(1|Subject),data=dpSum13)
summary(d3choice)
confint(d3choice,method='Wald')

#Plot
delaySum=summaryBy(smoothDat~time*wire*cluster,data=smoothDat,FUN=c(nanmean))
head(delaySum)
names(delaySum)=c('time','wire','cluster','zFR')

delaySum2=summaryBy(zFR~time*cluster,data=delaySum,FUN=c(mean,se,length))
names(delaySum2)=c('time','cluster','zFR','se','len')
delaySum2$cluster=as.factor(delaySum2$cluster)

delaySum2Trim=delaySum2[delaySum2$time>=500 & delaySum2$time<=2500,]
clusPlot1=ggplot() + theme_bw() + geom_line(data=delaySum2Trim,aes(x=time,y=zFR,color=cluster)) + geom_ribbon(data=delaySum2Trim,aes(x=time,ymin=zFR-se,ymax=zFR+se,fill=cluster),alpha=.3) + geom_vline(xintercept=1500,color='black',linetype='solid') + scale_y_continuous(limits=c(-2,3),breaks=c(-2,-1,0,1,2,3)) + theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank()) + scale_x_continuous(limits=c(500,2500),breaks=c(500,1000,1500,2000,2500))
clusPlot1
plotName=sprintf('%s/pccMicro_valDM_response_byClust.pdf',outDir)
ggsave(plotName,height=7,width=10,clusPlot1)

# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  Clusters x wasRisky
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
#load the data
cname=sprintf('%s/pccMicro_valueBasedDM_byRisk_responseLocked_smooth_byWire.csv',baseDir)
respSum2 <- fread(cname, data.table = FALSE)

rsResp=respSum2[(respSum2$time>1500) & (respSum2$time<2000),]
rsResp2=summaryBy(zFR~wasRisky*wire*cluster,data=rsResp)
head(rsResp2)

rsWide <- reshape2::dcast(
  rsResp2,
  wire + cluster ~ wasRisky,
  value.var = "zFR.mean"
)
names(rsWide)=c('wire','cluster','safe','risky')

rsWide$diff=rsWide$risky-rsWide$safe
rsWide$cluster=as.factor(rsWide$cluster)
head(rsWide)

rsWide$Subject=substr(rsWide$wire,1,6)
rsWide$Subject=as.factor(rsWide$Subject)
str(rsWide)
rsWide$cluster=relevel(rsWide$cluster,ref='2')

rDat1=rsWide[rsWide$cluster=='1',]
rDat2=rsWide[rsWide$cluster=='2',]
rDat3=rsWide[rsWide$cluster=='3',]

r1lm1=lmer(diff~1+(1|Subject),data=rDat1)
summary(r1lm1)
confint(r1lm1,method='Wald')

r2lm1=lmer(diff~1+(1|Subject),data=rDat2)
summary(r2lm1)
confint(r2lm1,method='Wald')

r3lm1=lmer(diff~1+(1|Subject),data=rDat3)
summary(r3lm1) 
confint(r3lm1,method='Wald')

# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  Feedback - Feedback - Feedback - Feedback
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Load and format the data ───────────────────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

datName=sprintf('%s/pccMicro_valueBasedDM_feedbackLocked_smooth_byBlock.csv',baseDir)
feedProbe <- fread(datName, data.table = FALSE)

#Plot firing rate for clusters centered on feedback time
feedProbe2=summaryBy(smoothDat~time*cluster,data=feedProbe,FUN=c(mean,se))
names(feedProbe2)=c('time','cluster','zFR','se')
feedProbe2$cluster=as.factor(feedProbe2$cluster)

feedSumTrim=feedProbe2[feedProbe2$time>=500 & feedProbe2$time<=2500,]
clusPlot1=ggplot() + theme_bw() + geom_line(data=feedSumTrim,aes(x=time,y=zFR,color=cluster)) + geom_ribbon(data=feedSumTrim,aes(x=time,ymin=zFR-se,ymax=zFR+se,fill=cluster),alpha=.3) + geom_vline(xintercept=1500,color='black',linetype='solid') + scale_y_continuous(limits=c(-2,3),breaks=c(-2,-1,0,1,2,3)) + theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank()) + scale_x_continuous(limits=c(500,2500),breaks=c(500,1000,1500,2000,2500))
clusPlot1
plotName=sprintf('%s/pccMicro_valDM_feedback_byClust.pdf',outDir)
ggsave(plotName,height=7,width=10,clusPlot1)


## Stats for apriori window against baseline
feedPeriod2=feedProbe[(feedProbe$time>1500) & (feedProbe$time<2000),]
fpSum1b=summaryBy(fr~wire*cluster,data=feedPeriod2,FUN=c(mean))

str(fpSum1b)
fpSum1b$Subject=substr(fpSum1b$wire,1,6)
fpSum1b$Subject=as.factor(fpSum1b$Subject)
fpSum1b$cluster=as.factor(fpSum1b$cluster)
lmerFeed1b=lmer(fr.mean~cluster+(1|Subject),data=fpSum1b)
summary(lmerFeed1b)

fpSum11b=fpSum1b[fpSum1b$cluster=='1',]
f1choice=lmer(fr.mean~1+(1|Subject),data=fpSum11b)
summary(f1choice)
confint(f1choice,method='Wald')

fpSum12b=fpSum1b[fpSum1b$cluster=='2',]
f2choice=lmer(fr.mean~1+(1|Subject),data=fpSum12b)
summary(f2choice)
confint(f2choice,method='Wald')

fpSum13b=fpSum1b[fpSum1b$cluster=='3',]
f3choice=lmer(fr.mean~1+(1|Subject),data=fpSum13b)
summary(f3choice)
confint(f3choice,method='Wald')

## Do the stats centered on the window of obvious peak
feedPeriod=feedProbe[(feedProbe$time>1250) & (feedProbe$time<=1751),]
fpSum1=summaryBy(fr~wire*cluster,data=feedPeriod,FUN=c(mean))

str(fpSum1)
fpSum1$Subject=substr(fpSum1$wire,1,6)
fpSum1$Subject=as.factor(fpSum1$Subject)
fpSum1$cluster=as.factor(fpSum1$cluster)
lmerFeed1=lmer(fr.mean~cluster+(1|Subject),data=fpSum1)
summary(lmerFeed1)

fpSum11=fpSum1[fpSum1$cluster=='1',]
f1choice=lmer(fr.mean~1+(1|Subject),data=fpSum11)
summary(f1choice)
confint(f1choice,method='Wald')

fpSum12=fpSum1[fpSum1$cluster=='2',]
f2choice=lmer(fr.mean~1+(1|Subject),data=fpSum12)
summary(f2choice)
confint(f2choice,method='Wald')

fpSum13=fpSum1[fpSum1$cluster=='3',]
f3choice=lmer(fr.mean~1+(1|Subject),data=fpSum13)
summary(f3choice)
confint(f3choice,method='Wald')

# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
# ──  Clusters x wasRisky
# ──  ─────────────────────────────────────────────────
# ──  ─────────────────────────────────────────────────
#load the data
cname=sprintf('%s/pccMicro_valueBasedDM_byRisk_feedbackLocked_smooth_byWire.csv',baseDir)
feedSum2 <- fread(cname, data.table = FALSE)

rsResp=feedSum2[(feedSum2$time>1500) & (feedSum2$time<2000),]
rsResp2=summaryBy(zFR~wasRisky*wire*cluster,data=rsResp)
head(rsResp2)

rsWide <- reshape2::dcast(
  rsResp2,
  wire + cluster ~ wasRisky,
  value.var = "zFR.mean"
)
names(rsWide)=c('wire','cluster','safe','risky')

rsWide$diff=rsWide$risky-rsWide$safe
rsWide$cluster=as.factor(rsWide$cluster)
head(rsWide)

rsWide$Subject=substr(rsWide$wire,1,6)
rsWide$Subject=as.factor(rsWide$Subject)
str(rsWide)
rsWide$cluster=relevel(rsWide$cluster,ref='2')

fDat1=rsWide[rsWide$cluster=='1',]
fDat2=rsWide[rsWide$cluster=='2',]
fDat3=rsWide[rsWide$cluster=='3',]

f1lm1=lmer(diff~1+(1|Subject),data=fDat1)
summary(f1lm1)
confint(f1lm1,method='Wald')

f2lm1=lmer(diff~1+(1|Subject),data=fDat2)
summary(f2lm1)
confint(f2lm1,method='Wald')

f3lm1=lmer(diff~1+(1|Subject),data=fDat3)
summary(f3lm1) 
confint(f3lm1,method='Wald')


#Do feedback-centered for cluster 1 like main effect
rsResp=feedSum2[(feedSum2$time>1250) & (feedSum2$time<=1751),]
rsResp2=summaryBy(zFR~wasRisky*wire*cluster,data=rsResp)
head(rsResp2)

rsWide <- reshape2::dcast(
  rsResp2,
  wire + cluster ~ wasRisky,
  value.var = "zFR.mean"
)
names(rsWide)=c('wire','cluster','safe','risky')

rsWide$diff=rsWide$risky-rsWide$safe
rsWide$cluster=as.factor(rsWide$cluster)
head(rsWide)

rsWide$Subject=substr(rsWide$wire,1,6)
rsWide$Subject=as.factor(rsWide$Subject)
str(rsWide)
rsWide$cluster=relevel(rsWide$cluster,ref='2')

fDat1=rsWide[rsWide$cluster=='1',]
fDat2=rsWide[rsWide$cluster=='2',]
fDat3=rsWide[rsWide$cluster=='3',]

f1lm1=lmer(diff~1+(1|Subject),data=fDat1)
summary(f1lm1)
confint(f1lm1,method='Wald')

f2lm1=lmer(diff~1+(1|Subject),data=fDat2)
summary(f2lm1)
confint(f2lm1,method='Wald')

f3lm1=lmer(diff~1+(1|Subject),data=fDat3)
summary(f3lm1) 
confint(f3lm1,method='Wald')

##################
##################
### Plot risk-modulation across timeperiods
##################
##################
#load the data
rm(list=ls())

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/valueDM'

rsName=sprintf('%s/pccMicro_valueBasedDM_byRisk_allTimesSummary_byWire.csv',baseDir)
riskPlotDat <- fread(rsName, data.table = FALSE)

rsWide <- reshape2::dcast(
  riskPlotDat,
  wire + cluster + timeperiod ~ wasRisky,
  value.var = "zFR.mean"
)
names(rsWide)=c('wire','cluster','timeperiod','safe','risky')

rsWide$diff=rsWide$risky-rsWide$safe
rsWide$cluster=as.factor(rsWide$cluster)
head(rsWide)

rsWide$plotX=NaN
rsWide$plotX=ifelse(rsWide$timeperiod=='onset' & rsWide$cluster=='1',2,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='onset' & rsWide$cluster=='2',3,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='onset' & rsWide$cluster=='3',4,rsWide$plotX)

rsWide$plotX=ifelse(rsWide$timeperiod=='response' & rsWide$cluster=='1',6,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='response' & rsWide$cluster=='2',7,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='response' & rsWide$cluster=='3',8,rsWide$plotX)

rsWide$plotX=ifelse(rsWide$timeperiod=='feedback' & rsWide$cluster=='1',10,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='feedback' & rsWide$cluster=='2',11,rsWide$plotX)
rsWide$plotX=ifelse(rsWide$timeperiod=='feedback' & rsWide$cluster=='3',12,rsWide$plotX)

rsWide2=summaryBy(diff~cluster*plotX*timeperiod,data=rsWide,FUN=c(mean,se))
srP1=ggplot() + theme_bw() + geom_bar(data=rsWide2,position='dodge',stat='identity',aes(x=plotX,y=diff.mean,fill=cluster)) + geom_errorbar(data=rsWide2,aes(x=plotX,ymin=diff.mean-diff.se,ymax=diff.mean+diff.se),width=.3) + geom_jitter(data=rsWide,aes(x=plotX,y=diff),color='grey30',alpha=.5,height=0,width=.1,size=3) + scale_y_continuous(limits=c(-10,10),breaks=c(-10,-5,0,5,10))
srP1
plotName=sprintf('%s/pccMicro_valDM_riskMod_byClust_allTimes.pdf',outDir)
ggsave(plotName,height=7,width=10,srP1)

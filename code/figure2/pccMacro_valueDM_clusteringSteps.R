rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(zoo)
library(dendextend)
library(cluster)
library(pheatmap)
library(colorspace)
library(RColorBrewer)

nanmean=function(x)(mean(x,na.rm=T))

se <- function(x) {
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(se=se)
}


#Set the data and output figure paths
outDir='<path2folder>/figures/Figure2'
baseDir='<path2folder>/data/valueDM'

#Remake the entire array to just have that
pccDat=read.csv(sprintf('%s/pccMacro_valueBasedDM_concatTimeseries.csv',baseDir),header=T)
head(pccDat)

#Take just the first block for all electrodes
df_min_block <- pccDat[pccDat$block == ave(pccDat$blockValue, pccDat$elecID, FUN = min),]
str(df_min_block)

# Average across trials for each electrode from that block
choiceProbe=summaryBy(psc~time*elecID,data=pccDat,FUN=c(mean))
names(choiceProbe)=c('time','elecID','psc')

#Smooth the first block data
# Loop through each wire and cluster and plot the choice period time course
smoothDat=data.frame()
wireList=unique(choiceProbe$elecID)
library(zoo)
for(clustI in wireList){
  clustDatC=choiceProbe[choiceProbe$elecID==clustI,]
  
  clustDatC$smoothDat=rollapply(clustDatC$psc, width=50, FUN=mean, align="center", partial=TRUE)
  smoothDat=rbind(smoothDat,clustDatC)
  
}

#This is where the clustering is done. Previous run-throughs helped to identify the best number of clusters for the data, this just runs through now showing those methods.
wideDat <- dcast(smoothDat, elecID ~ time, value.var = "smoothDat")
head(wideDat)
str(wideDat)

pccElecArr=unique(pccDat[,c("elecID",'locName')])
locIDarr=merge(pccElecArr,wideDat,by='elecID')
locIDinfo=locIDarr[,c('elecID','locName')]

#Now create the correlation dissimilarity matrix
corr_mat <- cor(t(as.matrix(wideDat[,-1])), use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - corr_mat)

#First do hierarchical clustering using a dendritic tree approach
hc <- hclust(dist_mat, method = "average")

my_colors <- c("#D19C63","#9B86B7","#C2B458","#6BAE7E","magenta","lightblue","green","grey","yellow","pink")
my_colors_rev=rev(my_colors)

dend <- as.dendrogram(hc)
clusters <- cutree(hc, k = 6)
# Color branches by cluster
dend_colored <- color_branches(dend, k = 6,col=my_colors)
dend_colored <- color_labels(dend_colored, k = 6,col = my_colors)
dend_colored <- set(dend_colored, "branches_lwd", 2)
plot(dend_colored,horiz=T, main = "Hierarchical clustering with colored branches")
dendName=sprintf('%s/pccMacro_valueDM_hcClust.pdf',outDir)
dev.copy(pdf, dendName)
dev.off()

#Plot dissimilarity matrix for visual inspection
base_cols <- brewer.pal(9, "YlOrRd")
base_cols[7:9] <- desaturate(base_cols[7:9], amount = 0.3)
cols <- colorRampPalette(c("#f7f7f7", base_cols))(100)

pheatmap(dist_mat,
         legend_breaks = c(0, 1, 1.8), legend = T,
         main = "", legend_labels = c("0", "1", "2"),
         cluster_rows = T,           
         cluster_cols = T,
         treeheight_row = 0, treeheight_col = 0,
         color = cols)
disMatName1=sprintf('%s/pccMacro_valueDM_dissMatrix.pdf',outDir)
dev.copy(pdf, disMatName1)
dev.off()

wideDat$dendClust=clusters
wideDat$dendClust=as.factor(wideDat$dendClust)


## Clustering Assessments
cmd <- cmdscale(dist_mat, k = 8)
wss <- numeric(15)
for (k in 1:15) {
  km <- kmeans(cmd, centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

plot(1:15, wss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares")


sil_width <- numeric(15)
for (k in 2:15) {
  km <- kmeans(cmd, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(cmd))
  sil_width[k] <- mean(ss[, 3])
}

plot(2:15, sil_width[2:15], type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Average silhouette width")

gap <- clusGap(cmd, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
plot(gap)



clusters_kmeans <- kmeans(cmd, centers = 4)$cluster
wideDat$kClust=clusters_kmeans
wideDat$kClust=as.factor(wideDat$kClust)

#PCA assessment and plotting
idLen=length(wireList)
pca <- prcomp(dist_mat, scale. = TRUE)

pcaDat=data.frame(pc1=numeric(idLen),pc2=numeric(idLen),clust=numeric(idLen))
pcaDat$pc1=as.numeric(pca$x[,1])
pcaDat$pc2=as.numeric(pca$x[,2])
pcaDat$clust=clusters
pcaDat$clust=as.factor(pcaDat$clust)

pcaDat$locName=locIDinfo$locName
pcaDat$locName=as.factor(pcaDat$locName)

str(pcaDat)

kmPlot=ggplot() + theme_bw() + geom_point(data=pcaDat,aes(x=pc1,y=pc2,fill=clust),pch=21,size=6,stroke=.5) + scale_fill_manual(values=my_colors_rev)
kmPlot
plot1Name=sprintf('%s/pccMacro_valueDM_pcaPlot_dim12.pdf',outDir)
ggsave(plot1Name,height=7,width=10,kmPlot)



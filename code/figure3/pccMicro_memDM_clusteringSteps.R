# ── Initiate Script ───────────────────────────────────────────────
rm(list=ls())
library(doBy)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(pheatmap)
library(colorspace)
library(RColorBrewer)
library(cluster)

get_subject <- function(wireVal) {
  sub("_.*", "", wireVal)
}

nanmean=function(x)(mean(x,na.rm=T))
se=function(x)(sd(x,na.rm = T)/sqrt(length(x)))

# ── Step 1: Load data ───────────────────────────────────
outDir='<path2folder>/figures/Figure3'
baseDir='<path2folder>/data/memoryDM'

microDat=read.csv(sprintf('%s/pccMicro_memoryBasedDM_concatTimeseries.csv',baseDir),header=TRUE)

# ── Step 2: Extract training block rows ───────────────────────────────────
smoothTrain=smoothDat[smoothDat$is_train==TRUE,]
wideDat <- dcast(smoothTrain, wire ~ time, value.var = "smoothFR")

# ── Step 3: Dissimilarity matrix ───────────────────────────────────
### Create the dissimilarity matrix
train_mat=as.matrix(wideDat[,-1])
corr_mat <- cor(t(as.matrix(wideDat[,-1])), use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - corr_mat)

hc <- hclust(dist_mat, method = "average")

base_cols <- brewer.pal(9, "YlOrRd")
base_cols[7:9] <- desaturate(base_cols[7:9], amount = 0.3)
cols <- colorRampPalette(c("#f7f7f7", base_cols))(100)

pheatmap(dist_mat,
         legend_breaks = c(0, 1, 1.8), legend = T,
         main = "", legend_labels = c("0", "1", "2"),
         cluster_rows = hc,           # keep kmeans order
         cluster_cols = hc,
         treeheight_row = 0, treeheight_col = 0,
         color = cols)           # cluster features if you wan
#save dissimilarity matrix
disMatName1=sprintf('%s/pccMicro_memoryDM_dissMatrix.pdf',outDir)
dev.copy(pdf, disMatName1)
dev.off()

# ── Step 4: Clustering metrics ───────────────────────────────────
# Silhouette width
K_max <- 12
sil_scores <- sapply(2:K_max, function(k) {
  pam_k   <- pam(dist_mat, k = k, diss = TRUE)
  sil     <- silhouette(pam_k$clustering, dist_mat)
  mean(sil[, "sil_width"])
})

sil_results <- data.frame(
  k         = 2:K_max,
  mean_sil  = sil_scores
)

plot(sil_results$k, sil_results$mean_sil,
     type = "b", pch = 19,
     xlab = "K", ylab = "Mean silhouette width",
     main = "Silhouette width",
     xaxt = "n")
axis(1, at = 2:K_max)

## Elbow method
wss <- sapply(1:10, function(k) {
  cluster_assign <- cutree(hc, k)
  sum(sapply(unique(cluster_assign), function(cl) {
    members <- which(cluster_assign == cl)
    sum(as.matrix(dist_mat)[members, members]) / 2
  }))
})
plot(1:10, wss, type = "b", xlab = "Number of clusters (k)",
     ylab = "Within-cluster dissimilarity")


# Dendrogram
clusters <- cutree(hc, k =4)  # for X clusters
wideDat$clusters=clusters

clustDF <- data.frame(
  wire = wideDat$wire,
  cluster = clusters
)
clustDF$wire=as.factor(clustDF$wire)

my_colors <- c("#cae0cf", "#89ad91", "#506858","#a8a8a8","lightblue","magenta","green","grey","yellow","pink")

library(dendextend)
dend <- as.dendrogram(hc)
clusters <- cutree(hc, k = 4)
# Color branches by cluster
dend_colored <- color_branches(dend, k = 4,col=my_colors)
dend_colored <- color_labels(dend_colored, k = 4,col = my_colors)
dend_colored <- set(dend_colored, "branches_lwd", 1)
plot(dend_colored,horiz=T, main = "Hierarchical clustering with colored branches",leaflab = "none")
dendName=sprintf('%s/pccMicro_memoryDM_hcClust.pdf',outDir)
dev.copy(pdf, dendName)
dev.off()
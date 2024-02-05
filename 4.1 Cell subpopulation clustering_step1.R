rm(list=ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(ConsensusClusterPlus)
library(heatmap3)
library(Seurat)


#MAD提取细胞间高变基因-------------------------------------
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')
BC_Epithelial <- subset(BC,idents ="Epithelial")

BC1_Tumor <- subset(BC_Epithelial,patient=="BC1")
BC1_Tumor <- GetAssayData(object = BC1_Tumor, slot = "data")
BC1_Tumor <- as.matrix(BC1_Tumor)

BC2_Tumor <- subset(BC_Epithelial,patient=="BC2")
BC2_Tumor <- GetAssayData(object = BC2_Tumor, slot = "data")
BC2_Tumor <- as.matrix(BC2_Tumor)

BC3_Tumor <- subset(BC_Epithelial,patient=="BC3")
BC3_Tumor <- GetAssayData(object = BC3_Tumor, slot = "data")
BC3_Tumor <- as.matrix(BC3_Tumor)

BC4_Tumor <- subset(BC_Epithelial,patient=="BC4")
BC4_Tumor <- GetAssayData(object = BC4_Tumor, slot = "data")
BC4_Tumor <- as.matrix(BC4_Tumor)

BC5_Tumor <- subset(BC_Epithelial,patient=="BC5")
BC5_Tumor <- GetAssayData(object = BC5_Tumor, slot = "data")
BC5_Tumor <- as.matrix(BC5_Tumor)

BC6_Tumor <- subset(BC_Epithelial,patient=="BC6")
BC6_Tumor <- GetAssayData(object = BC6_Tumor, slot = "data")
BC6_Tumor <- as.matrix(BC6_Tumor)

BC7_Tumor <- subset(BC_Epithelial,patient=="BC7")
BC7_Tumor <- GetAssayData(object = BC7_Tumor, slot = "data")
BC7_Tumor <- as.matrix(BC7_Tumor)

BC1_Tumor <- BC1_Tumor[names(sort(apply(BC1_Tumor,1,mad),decreasing = T))[1:5000],]
BC2_Tumor <- BC2_Tumor[names(sort(apply(BC2_Tumor,1,mad),decreasing = T))[1:5000],]
BC3_Tumor <- BC3_Tumor[names(sort(apply(BC3_Tumor,1,mad),decreasing = T))[1:5000],]
BC4_Tumor <- BC4_Tumor[names(sort(apply(BC4_Tumor,1,mad),decreasing = T))[1:5000],]
BC5_Tumor <- BC5_Tumor[names(sort(apply(BC5_Tumor,1,mad),decreasing = T))[1:5000],]
BC6_Tumor <- BC6_Tumor[names(sort(apply(BC6_Tumor,1,mad),decreasing = T))[1:5000],]
BC7_Tumor <- BC7_Tumor[names(sort(apply(BC7_Tumor,1,mad),decreasing = T))[1:5000],]

DataList <- list(BC1=BC1_Tumor,BC2=BC2_Tumor,BC3=BC3_Tumor,BC4=BC4_Tumor,
                 BC5=BC5_Tumor,BC6=BC6_Tumor,BC7=BC7_Tumor)

#执行单样本亚型构建
source('4.1 SCC.step1.R')
SCC.step1(NESlist = DataList,data.path = 'Results/4.1 Cell subpopulation clustering_step1/Rdata1',figure.path = 'Results/4.1 Cell subpopulation clustering_step1/Figures1')


# #查看单样本的cluster  
# load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC5_Consensus.rda')
# clusterNum=3
# cluster=results[[clusterNum]][["consensusClass"]]  #从results中提取样本分型信息
# sub <- data.frame(Sample=names(cluster),Cluster=cluster)  
# sub$Cluster <- paste0('C',sub$Cluster)##前面加上C
# table(sub$Cluster)

# #####PAC确定最佳聚类数  Proportion Ambiguous Cluster #####------------------
# library(diceR)
# PACK <- function(conslist){
#   PACbestK <- NULL
#   for (i in names(conslist)) {
#     results <- conslist[[i]]
#     pac <- lapply(results[2:length(results)], function(x){
#       PAC(x$consensusMatrix,lower = 0.05, upper = 0.95)
#     })
#     PACbestK <- c(PACbestK,which.min(as.numeric(pac))+1)
#   }
#   names(PACbestK) <- names(conslist)
#   return(PACbestK)
# }
# 
# load('6.ithCluster/Rdata/BC1_Consensus.rda')
# BC1 <- results
# load('6.ithCluster/Rdata/BC2_Consensus.rda')
# BC2 <- results
# load('6.ithCluster/Rdata/BC3_Consensus.rda')
# BC3 <- results
# load('6.ithCluster/Rdata/BC6_Consensus.rda')
# BC6 <- results
# load('6.ithCluster/Rdata/BC7_Consensus.rda')
# BC7 <- results
# conslist <- list(BC1=BC1,BC2=BC2,BC3=BC3,BC6=BC6,BC7=BC7)
# PACK(conslist)
# 
# 
# #calinsky确定最佳聚类数-------------------------
# calinsky <- function(hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
#   msg <- ""
#   if (is.null(dist)) {
#     require(clue)
#     dist <- sqrt(as.cl_ultrametric(hhc))
#     # message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
#   } else if (attr(dist, "method") != "euclidean") {
#     require(clue)
#     dist <- sqrt(as.cl_ultrametric(hhc))
#     # message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
#   }
#   dist <- as.matrix(dist)^2
#   A <- -dist/2
#   A_bar <- apply(A, 1, mean)
#   totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
#   n <- length(hhc$order)
#   ans <- rep(0, gMax)
#   for (g in 2:gMax) {
#     cclust <- cutree(hhc, k = g)
#     withinSum <- 0
#     for (k in 1:g) {
#       if (sum(cclust == k) == 1) 
#         next
#       A <- as.matrix(-dist/2)[cclust == k, cclust == k]
#       A_bar <- apply(A, 1, mean)
#       withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
#                                      mean(A))
#     }
#     betweenSum <- totalSum - withinSum
#     betweenSum <- betweenSum/(g - 1)
#     withinSum <- withinSum/(n - g)
#     ans[g] <- betweenSum/withinSum
#   }
#   class(ans) <- "calinsky"
#   attr(ans, "message") <- msg
#   return(ans)
# }
# 
# for(i in 1:(length(conslist) - 1)){
#   wData <- t(conslist[[i]])
#   ddist <- dist(wData, method = "euclidean")
#   sHc <- hclust(ddist, method = "ward.D2")
#   aCalinsky <- calinsky(hhc=sHc,dist = NULL, gMax = 10)
#   save(aCalinsky,file = paste0("6.ithCluster/aCalinsky",names(DataList[i],".rda")))
# }
# 
# #查看结果
# load('6.ithCluster/aCalinskyBC1.rda')
# nCls <- which.max(aCalinsky)
# 
# 
# 
# 




# load("3.CNV/BC_T1.rda")
# load("3.CNV/BC_T2.rda")
# load("3.CNV/BC_T3.rda")
# load("3.CNV/BC_T4.rda")
# load("3.CNV/BC_T5.rda")
# load("3.CNV/BC_T6.rda")
# load("3.CNV/BC_T7.rda")
# 
# #构建单样本的肿瘤细胞表达矩阵列表
# T1 <- read.table(file = '3.CNV/sample1_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T1,kmeans_class=='1')
# class3 <- subset(T1,kmeans_class=='3')
# class5 <- subset(T1,kmeans_class=='5')
# Tumor1 <- rbind(class1,class3,class5)
# Tumor1 <- intersect(rownames(Tumor1),colnames(BC_T1))
# BC1_Tumor <- BC_T1[,Tumor1]
# #对数据进行标准化并提取
# Seurat_Tumor1 <- CreateSeuratObject(counts = BC1_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor1 <- NormalizeData(Seurat_Tumor1)
# BC1_Tumor <- GetAssayData(object = Seurat_Tumor1, slot = "data")
# BC1_Tumor <- as.matrix(BC1_Tumor)
# range(BC1_Tumor)
# 
# T2 <- read.table(file = '3.CNV/sample2_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T2,kmeans_class=='1')
# class5 <- subset(T2,kmeans_class=='5')
# class6 <- subset(T2,kmeans_class=='6')
# Tumor2 <- rbind(class1,class5,class6)
# Tumor2 <- intersect(rownames(Tumor2),colnames(BC_T2))
# BC2_Tumor <- BC_T2[,Tumor2]
# #对数据进行标准化并提取
# Seurat_Tumor2 <- CreateSeuratObject(counts = BC2_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor2 <- NormalizeData(Seurat_Tumor2)
# BC2_Tumor <- GetAssayData(object = Seurat_Tumor2, slot = "data")
# BC2_Tumor <- as.matrix(BC2_Tumor)
# range(BC2_Tumor)
# 
# T3 <- read.table(file = '3.CNV/sample3_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T3,kmeans_class=='1')
# class5 <- subset(T3,kmeans_class=='5')
# class6 <- subset(T3,kmeans_class=='6')
# class8 <- subset(T3,kmeans_class=='8')
# Tumor3 <- rbind(class1,class5,class6,class8)
# Tumor3 <- intersect(rownames(Tumor3),colnames(BC_T3))
# BC3_Tumor <- BC_T3[,Tumor3]
# #对数据进行标准化并提取
# Seurat_Tumor3 <- CreateSeuratObject(counts = BC3_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor3 <- NormalizeData(Seurat_Tumor3)
# BC3_Tumor <- GetAssayData(object = Seurat_Tumor3, slot = "data")
# BC3_Tumor <- as.matrix(BC3_Tumor)
# range(BC3_Tumor)
# 
# T4 <- read.table(file = '3.CNV/sample4_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T4,kmeans_class=='1')
# class5 <- subset(T4,kmeans_class=='5')
# class6 <- subset(T4,kmeans_class=='6')
# Tumor4 <- rbind(class1,class5,class6)
# Tumor4 <- intersect(rownames(Tumor4),colnames(BC_T4))
# BC4_Tumor <- BC_T4[,Tumor4]
# #对数据进行标准化并提取
# Seurat_Tumor4 <- CreateSeuratObject(counts = BC4_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor4 <- NormalizeData(Seurat_Tumor4)
# BC4_Tumor <- GetAssayData(object = Seurat_Tumor4, slot = "data")
# BC4_Tumor <- as.matrix(BC4_Tumor)
# range(BC4_Tumor)
# 
# T5 <- read.table(file = '3.CNV/sample5_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T5,kmeans_class=='1')
# class2 <- subset(T5,kmeans_class=='2')
# class5 <- subset(T5,kmeans_class=='5')
# class6 <- subset(T5,kmeans_class=='6')
# class7 <- subset(T5,kmeans_class=='7')
# Tumor5 <- rbind(class1,class2,class5,class6,class7)
# Tumor5 <- intersect(rownames(Tumor5),colnames(BC_T5))
# BC5_Tumor <- BC_T5[,Tumor5]
# #对数据进行标准化并提取
# Seurat_Tumor5 <- CreateSeuratObject(counts = BC5_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor5 <- NormalizeData(Seurat_Tumor5)
# BC5_Tumor <- GetAssayData(object = Seurat_Tumor5, slot = "data")
# BC5_Tumor <- as.matrix(BC5_Tumor)
# range(BC5_Tumor)
# 
# T6 <- read.table(file = '3.CNV/sample6_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T6,kmeans_class=='1')
# class2 <- subset(T6,kmeans_class=='2')
# class6 <- subset(T6,kmeans_class=='6')
# Tumor6 <- rbind(class1,class2,class6)
# Tumor6 <- intersect(rownames(Tumor6),colnames(BC_T6))
# BC6_Tumor <- BC_T6[,Tumor6]
# #对数据进行标准化并提取
# Seurat_Tumor6 <- CreateSeuratObject(counts = BC6_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor6 <- NormalizeData(Seurat_Tumor6)
# BC6_Tumor <- GetAssayData(object = Seurat_Tumor6, slot = "data")
# BC6_Tumor <- as.matrix(BC6_Tumor)
# range(BC6_Tumor)
# 
# T7 <- read.table(file = '3.CNV/sample7_kmeans_df_s.txt', sep = '\t')
# class1 <- subset(T7,kmeans_class=='1')
# class2 <- subset(T7,kmeans_class=='2')
# class6 <- subset(T7,kmeans_class=='6')
# class7 <- subset(T7,kmeans_class=='7')
# class8 <- subset(T7,kmeans_class=='8')
# Tumor7 <- rbind(class1,class2,class6,class7,class8)
# Tumor7 <- intersect(rownames(Tumor7),colnames(BC_T7))
# BC7_Tumor <- BC_T7[,Tumor7]
# #对数据进行标准化并提取
# Seurat_Tumor7 <- CreateSeuratObject(counts = BC7_Tumor, project = "seurat",min.cells = 3)
# Seurat_Tumor7 <- NormalizeData(Seurat_Tumor7)
# BC7_Tumor <- GetAssayData(object = Seurat_Tumor7, slot = "data")
# BC7_Tumor <- as.matrix(BC7_Tumor)
# range(BC7_Tumor)
# 
# save(BC1_Tumor,BC2_Tumor,BC3_Tumor,BC4_Tumor,BC5_Tumor,BC6_Tumor,BC7_Tumor,file = '4.step1/BC_Tumor.rda')



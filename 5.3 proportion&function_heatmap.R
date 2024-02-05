rm(list = ls())
library(data.table)
library(tibble)
library(Biobase)
library(SeuratObject)
library(sp)
library(Seurat)
library(tidyr)

##合并不同样本的细胞亚群-------------------------------------
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC1_Consensus.rda')
results_BC1 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC2_Consensus.rda')
results_BC2 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC3_Consensus.rda')
results_BC3 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC4_Consensus.rda')
results_BC4 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC5_Consensus.rda')
results_BC5 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC6_Consensus.rda')
results_BC6 <- results
load('Results/4.1 Cell subpopulation clustering_step1/Rdata1/BC7_Consensus.rda')
results_BC7 <- results

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

dataList <- list(BC1 = list(consRES = results_BC1, geData = BC1_Tumor, bestK = 4),
                 BC2 = list(consRES = results_BC2, geData = BC2_Tumor, bestK = 4),
                 BC3 = list(consRES = results_BC3, geData = BC3_Tumor, bestK = 4),
                 BC4 = list(consRES = results_BC4, geData = BC4_Tumor, bestK = 3),
                 BC5 = list(consRES = results_BC5, geData = BC5_Tumor, bestK = 3),
                 BC6 = list(consRES = results_BC6, geData = BC6_Tumor, bestK = 5),
                 BC7 = list(consRES = results_BC7, geData = BC7_Tumor, bestK = 3))

load('Results/4.3 Cell subpopulation clustering_step3/Rdata3/aConsensus_allClusters.rda')

G <- dataList[[1]]$bestK
aConsensus <- dataList[[1]]$consRES[[G]]
ss_cluster1 <- as.data.frame(aConsensus$consensusClass)
ss_cluster1$`aConsensus$consensusClass` <- paste0('BC1_c',ss_cluster1$`aConsensus$consensusClass`)
colnames(ss_cluster1) <- "cluster"

G <- dataList[[2]]$bestK
aConsensus <- dataList[[2]]$consRES[[G]]
ss_cluster2 <- as.data.frame(aConsensus$consensusClass)
ss_cluster2$`aConsensus$consensusClass` <- paste0('BC2_c',ss_cluster2$`aConsensus$consensusClass`)
colnames(ss_cluster2) <- "cluster"

G <- dataList[[3]]$bestK
aConsensus <- dataList[[3]]$consRES[[G]]
ss_cluster3 <- as.data.frame(aConsensus$consensusClass)
ss_cluster3$`aConsensus$consensusClass` <- paste0('BC3_c',ss_cluster3$`aConsensus$consensusClass`)
colnames(ss_cluster3) <- "cluster"

G <- dataList[[4]]$bestK
aConsensus <- dataList[[4]]$consRES[[G]]
ss_cluster4 <- as.data.frame(aConsensus$consensusClass)
ss_cluster4$`aConsensus$consensusClass` <- paste0('BC4_c',ss_cluster4$`aConsensus$consensusClass`)
colnames(ss_cluster4) <- "cluster"

G <- dataList[[5]]$bestK
aConsensus <- dataList[[5]]$consRES[[G]]
ss_cluster5 <- as.data.frame(aConsensus$consensusClass)
ss_cluster5$`aConsensus$consensusClass` <- paste0('BC5_c',ss_cluster5$`aConsensus$consensusClass`)
colnames(ss_cluster5) <- "cluster"

G <- dataList[[6]]$bestK
aConsensus <- dataList[[6]]$consRES[[G]]
ss_cluster6 <- as.data.frame(aConsensus$consensusClass)
ss_cluster6$`aConsensus$consensusClass` <- paste0('BC6_c',ss_cluster6$`aConsensus$consensusClass`)
colnames(ss_cluster6) <- "cluster"

G <- dataList[[7]]$bestK
aConsensus <- dataList[[7]]$consRES[[G]]
ss_cluster7 <- as.data.frame(aConsensus$consensusClass)
ss_cluster7$`aConsensus$consensusClass` <- paste0('BC7_c',ss_cluster7$`aConsensus$consensusClass`)
colnames(ss_cluster7) <- "cluster"

cellClust <- rbind(ss_cluster1,ss_cluster2,ss_cluster3,ss_cluster4,ss_cluster5,ss_cluster6,ss_cluster7)

load('Results/4.3 Cell subpopulation clustering_step3/Rdata3/aConsensus_allClusters.rda')
cluster <- data.frame(cluster = aConsensus$consensusClass)
cluster$cluster <- paste0('BCIS',cluster$cluster)
cellClust$cells <- rownames(cellClust)
cellClust <- merge( cluster,cellClust, by.x = 0, by.y = 1) %>%
  column_to_rownames('cells')
colnames(cellClust) <- c('SS_cluster','Cluster')
table(cellClust$Cluster)

save(cellClust,file = "Results/5. BCIS_function/cell_Cluster.rda")

#比例图---------------------------------------
rm(list = ls())
load("Results/5. BCIS_function/cell_Cluster.rda")
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')
BC_Epithelial <- subset(BC,idents ="Epithelial")

ssBC_cells <- subset(BC_Epithelial, cells = rownames(cellClust))
metadata <- ssBC_cells@meta.data
metadata <- merge(metadata,cellClust,by=0)
metadata <- column_to_rownames(metadata,'Row.names')
metadata <- arrange(metadata,Cluster)
ssBC_cells@meta.data <- metadata
ssBC_cells <- SetIdent(ssBC_cells,value = ssBC_cells@meta.data$Cluster) 


metadata <- ssBC_cells@meta.data
metadata <- metadata[,c(4,9)]
metadata <- rownames_to_column(metadata,'ID')
#metadata <- pivot_longer(metadata,cols = 2:3, names_to = 'patient', values_to = 'percent')
library(ggplot2)
library(ggsci)
ggplot( metadata, aes( x = patient,  fill = `Cluster`))+
  geom_bar( position = "stack")+
  ylab(label = 'Percent')+xlab(label = '')+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme_bw(base_size = 1,base_rect_size = 2)+
  theme(axis.title.y = element_text(size = 15,colour = '#b74306', face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13,colour = '#b74306',face = 'bold',angle = 45,hjust = 1.0),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1.5),
        panel.background = element_rect(fill = "#eefdf7", color = NA),
        # panel.grid = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.grid.minor.x = element_line(color = "#cacfd2", linetype = "dashed"),
        legend.text = element_text(size = 12,colour = 'black'),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA, color = NA))+
  scale_fill_manual(values = c('#ED2A32','#82E0BA','#FF733F','#7BC9D9'))
ggsave(filename = 'Results/5. BCIS_function/BCIS_barplot.pdf',width = 4,height = 5)


#单细胞亚型功能热图--------------------------------------------
library(dplyr)
load('Results/5. BCIS_function/NES_allcluster.rda')

allPath <- as.character(unlist(aDEA2))
dup <- sort(unique(allPath[duplicated(allPath)]))
allPath <- allPath[!allPath %in% dup]
aDEA2 <- lapply(aDEA2, function(x) x[x %in% allPath])
outPath <- Reduce(c, aDEA2)
outPath <- data.frame(path = outPath, cluster = rep(names(unlist(lapply(aDEA2, length))), as.numeric(unlist(lapply(aDEA2, length)))))
outPath2 <- outPath %>%
  group_by(cluster) %>%
  slice_head(n=20)

NES_cluster_feature <- NES_allClusters[outPath2$path,]

load("Results/4.3 Cell subpopulation clustering_step3/Rdata3/aConsensus_allClusters.rda", verbose = T)
load("Results/4.3 Cell subpopulation clustering_step3/Rdata3/aMwwGSTs_allClusters.rda", verbose = T)
levels(consClust)[levels(consClust) == "black"] <- "cyan"
levels(consClust)
levels(consClust) <- c("BCIS1", "BCIS2", "BCIS3", "BCIS4")
consClust_df <- as.data.frame(consClust)%>%tibble::rownames_to_column("ID")

cluster <- arrange(consClust_df,consClust)%>%tibble::column_to_rownames("ID")
dd <- NES_cluster_feature[,rownames(cluster)]

#热图
library(pheatmap)
ssgsea_matrix <- dd
#列注释信息
table(cluster$consClust)
annotation_col <- data.frame(
  subtype = c(rep("BCIS1", 8), rep("BCIS2", 6), rep("BCIS3", 7), rep("BCIS4", 5)))
rownames(annotation_col) <- colnames(ssgsea_matrix)

annColors <- list(c("BCIS1"="#436b95", "BCIS2"="#b81d25", "BCIS3"="#2f7e2f","BCIS4"="#800080")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annotation_col)[1] # 如果有多个变量要注释颜色请补充每张list的name

pdf("Results/5. BCIS_function/BCIS_pathway_heatmap20.pdf", width = 12, height = 20)
pheatmap(ssgsea_matrix, 
         scale = "none",
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_rows = FALSE,cluster_cols = FALSE,
         gaps_col = c(8,14,21),
         show_colnames = F,
         annotation_col = annotation_col)
dev.off()
















# 
# #单细胞亚型代谢-----------------------------
# rm(list = ls())
# options(stringsAsFactors = F)
# library(scMetabolism)
# library(Seurat)
# library(ggplot2)
# library(rsvd)
# library(pheatmap)
# library(data.table)
# library(tidyverse)
# library(ggsci)
# library(Seurat)
# library(ComplexHeatmap)
# library(circlize)
# #library(rJava)
# #library(xlsx)
# library(CMScaller)
# library(survival)
# library(survminer)
# library(tibble)
# 
# load("8.Bulk_survival/cell_Cluster.rda")
# load('1.单细胞数据处理/注释后BC.rda')
# BC_Epithelial <- subset(BC,idents ="Epithelial")
# 
# ssBC_cells <- subset(BC_Epithelial, cells = rownames(cellClust))
# metadata <- ssBC_cells@meta.data
# metadata <- merge(metadata,cellClust,by=0)
# metadata <- column_to_rownames(metadata,'Row.names')
# metadata <- arrange(metadata,Cluster)
# ssBC_cells@meta.data <- metadata
# ssBC_cells <- SetIdent(ssBC_cells,value = ssBC_cells@meta.data$Cluster) 
# 
# human_countexp_Seurat<-sc.metabolism.Seurat(obj = ssBC_cells,
#                                             method = "AUCell",  #AUCell
#                                             imputation =F, 
#                                             ncores = 10, 
#                                             metabolism.type = "KEGG")
# 
# save(human_countexp_Seurat,file = 'Results/5. BCIS_function/human_countexp_Seurat.rda')
# 
# load('Results/5. BCIS_function/human_countexp_Seurat.rda')
# input.pathway<-rownames(human_countexp_Seurat@assays[["METABOLISM"]][["score"]])[1:40]
# df$pathway=factor(df$pathway, levels = 一个你想要的顺序通路名称向量, 
#                   labels = 一个你想要的顺序通路名称向量)
# 
# DotPlot.metabolism(obj = human_countexp_Seurat, pathway = input.pathway, 
#                    phenotype = "cluster", norm = "y")
# ggsave(filename = '单细胞亚型代谢/metabolism.pdf',width = 8,height = 10)





rm(list=ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(qusage)

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

##进行功能评分-----------------------
pathways <- readRDS('Results/4.2 Cell subpopulation clustering_step2/Pathways.rds')
source('4.2 SCC.Step2.R')
SCC.step2(dataList,pathways,data.path = 'Results/4.2 Cell subpopulation clustering_step2/Rdata2')






# ###免做 MSigDB genesets ----------------------------------
# hm <- qusage::read.gmt('7.SCC.step2/Reference/h.all.v2023.1.Hs.symbols.gmt')
# kegg <- qusage::read.gmt('7.SCC.step2/Reference/c2.cp.kegg.v2023.1.Hs.symbols.gmt')
# go <- qusage::read.gmt('7.SCC.step2/Reference/c5.go.v2023.1.Hs.symbols.gmt')
# all <- c(hm, kegg, go) # 10768条通路的基因集
# ## 保留基因数量大于等于15条且小于等于200条的通路
# # Ref：2019 - Natrue Protocol - http://www.nature.com/articles/s41596-018-0103-9
# all_filtered <- all[unlist(lapply(all, function(x){length(x)>=15 & length(x)<=200}))]
# saveRDS(all_filtered, file = '7.SCC.step2/Pathways.rds')


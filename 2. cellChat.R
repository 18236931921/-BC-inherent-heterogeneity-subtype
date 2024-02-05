# library(devtools)
#devtools::install_github("sqjin/CellChat")
rm(list = ls())
options(stringsAsFactors = FALSE)
library(patchwork)
library(Seurat)
library(CellChat)
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')

#一、CellChat的数据输入-----------------
#细胞的标准化后的基因表达数据，另一个是用户分配的细胞标签
exp  <- BC@assays$RNA@data  #输入数据必须为matrix
dim(exp)
range(exp)
#exp = normalizeData(exp)
# head(exp)[,1:4]
# meta1 = BC@meta.data
meta2 =  data.frame(Cell = rownames(BC@meta.data), 
                    cell_type = BC@active.ident ) 
table(meta2$cell_type)

#创建CellChat 对象
cellchat <- createCellChat(object = exp, meta = meta2, group.by = "cell_type")

#设置受体配体数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB  # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

#预处理用于细胞间通讯的表达数据分析
cellchat <- subsetData(cellchat)
library(future)
future::plan("multisession", workers = 8)   #multiprocess
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#二、细胞间通信网络的推理---------------------------
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#提取推断的cellchat网络作为数据框架
df.net <- subsetCommunication(cellchat)

#在信号通路级别推断细胞-细胞通信
cellchat <- computeCommunProbPathway(cellchat)

#计算整合的细胞通信网络
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

save(cellchat,file = 'Results/2. cellChat/cellchat.rda')

#数量与强度差异热图
par(mfrow=c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat,measure = 'weight')
h1+h2

#三、细胞间通信网络的可视化,使用层次结构图、圆图可视化特定信号通路或弦图---------
cellchat@netP[["pathways"]]

pathways.show <- cellchat@netP[["pathways"]] 
vertex.receiver = seq(1,5)  
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# group.cellType <- c(rep("", 4), rep("", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#计算每个配体-受体对整体的贡献 信号通路和可视化由 单配体-受体对
netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# #自动保存所有推断网络的绘图以便快速勘探
# pathways.show.all <- cellchat@netP$pathways
# levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "chord")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0('3.5 细胞间通讯/',pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

#可视化由多个介导的细胞间通信配体受体或信号通路,气泡图,和弦图
netVisual_bubble(cellchat,sources.use = 1, targets.use = c(1:4),remove.isolate = F)   #sources.use指定目标亚群，targets.use指定与哪些亚群关系

netVisual_chord_gene(cellchat,sources.use = 1, targets.use = c(1:4),lab.cex = 0.5,legend.pos.y = 30)


#使用小提琴/点绘制信号基因表达分布 
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)


#四、细胞间通信网络的系统分析-------------------

#识别细胞的信号传导角色（例如，显性发送者、接收者）组以及主要贡献信号
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#可视化主要发送方（源）和接收方（目标） 2D 空间
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = NULL)
gg1 + gg2

#识别对传出或传入信令贡献最大的信号 某些细胞组
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#识别全局通信模式以探索多个细胞群如何和信号通路协调在一起----------

#识别并可视化目标的传出通信模式
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

#识别并可视化目标的传入通信模式
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

#根据功能相似性识别信号基团
library(CellChat)
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")  #Error: No such strategy for futures: ‘multiprocess’
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

#根据结构相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)  #Error: No such strategy for futures: ‘multiprocess’
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

# saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")









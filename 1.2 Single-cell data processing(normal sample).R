rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)

BCN <- fread(file = "Datasets/GSM5329919_BCN_gene_cell_exprs_table.txt")
BCN <- BCN[,-1]
BCN=aggregate(.~Symbol,max,data=BCN)
BCN <- BCN%>%
  tibble::column_to_rownames('Symbol')
patientN <- data.frame(rep('BCN',6035))
rownames(patientN) <- colnames(BCN)
colnames(patientN) <- "patient"
BCN <- CreateSeuratObject(counts = BCN, meta.data=patientN, project = "seurat",min.cells = 3)

#质控
BCN[["percent.mt"]] <- PercentageFeatureSet(BCN, pattern = "^MT-")
BCN <- subset(BCN, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

#SCTransform
BCN <- SCTransform(BCN)
#运行 PCA
BCN <- RunPCA(BCN)

#弯道图（Elbow plot）法,选择合适的PC数，一般选择拐点作为降维的维度---------------
ElbowPlot(BCN) 
ggsave(filename = 'Results/1.2 Single-cell data processing(normal sample)/BCN_ElbowPlot.pdf',width = 6,height = 4)

#选择最佳聚类数
BCN <- FindNeighbors(BCN, reduction = "pca", dims = 1:10)
BCN <- FindClusters(BCN)
head(Idents(BCN), 5) #Idents获取对象的标识类

#UMAP聚类
BCN <- RunUMAP(BCN, dims = 0:16)
DimPlot(BCN, reduction = "umap",label = T)+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:17]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:17]))+  
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = 'BCN')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold'))
ggsave(filename = 'Results/1.2 Single-cell data processing(normal sample)/umap.pdf',width = 8,height = 6)
#TSNE聚类
BCN <- RunTSNE(BCN, dims = 0:16)
DimPlot(BCN, reduction = "tsne",label = T)+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:17]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:17]))+  
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = 'BCN')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold'))
ggsave(filename = 'Results/1.2 Single-cell data processing(normal sample)/tsne.pdf',width = 6,height = 4)


#查看细胞亚群的差异基因--------------------------------------------------------
BCN_markers <- FindAllMarkers(BCN, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
save(BCN_markers,file = 'Results/1.2 Single-cell data processing(normal sample)/BCN_markers.rda')

load('1、单细胞数据处理/BC_markers.rda')
top <- BCN_markers %>% group_by(cluster) %>% slice_max(n = 20,order_by=avg_log2FC)

#细胞注释
new.cluster.ids <- c("Fibroblast","Epithelial","Endothelial","Fibroblast","Epithelial","Epithelial","Epithelial","Epithelial","T cell","Epithelial","Epithelial","Epithelial","Myeloid/Macrophage","Epithelial",
                     "Epithelial","Myeloid/Macrophage","Myeloid/Macrophage")
names(new.cluster.ids) <- levels(BCN)
BCN <- RenameIdents(BCN, new.cluster.ids)
DimPlot(BCN, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:5]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:5]))+  #对边框色
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = '')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
ggsave(filename = 'Results/1.2 Single-cell data processing(normal sample)/post-annotated_UMAP.pdf',width = 9,height = 6)

DimPlot(BCN, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:5]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:5]))+  #对边框色
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = '')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
ggsave(filename = 'Results/1.2 Single-cell data processing(normal sample)/post-annotated_TSNE.pdf',width = 6,height = 6)


#提取上皮细胞类群
BCN_Epithelial <- subset(BCN,idents ="Epithelial")
save(BCN_Epithelial,file = 'Results/1.2 Single-cell data processing(normal sample)/BCN_Epithelial.rda')







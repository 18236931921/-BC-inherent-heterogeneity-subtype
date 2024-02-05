rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)

#输入来自样本的肿瘤组织的单细胞数据，创建Seurat后合并-------------------
BC1 <- fread(file = "Datasets/GSM4006644_BC1_gene_cell_exprs_table.txt")
BC1 <- BC1[,-1]
BC1=aggregate(.~Symbol,max,data=BC1)
BC1 <- BC1%>%
  tibble::column_to_rownames('Symbol')
patient1 <- data.frame(rep('BC1',3218))
rownames(patient1) <- colnames(BC1)
colnames(patient1) <- "patient"
BC1 <- CreateSeuratObject(counts = BC1, meta.data=patient1, project = "seurat",min.cells = 3)

BC2 <- fread(file = "Datasets/GSM4006645_BC2_gene_cell_exprs_table.txt")
BC2 <- BC2[,-1]
BC2=aggregate(.~Symbol,max,data=BC2)
BC2 <- BC2%>%
  tibble::column_to_rownames('Symbol')
patient2 <- data.frame(rep('BC2',5006))
rownames(patient2) <- colnames(BC2)
colnames(patient2) <- "patient"
BC2 <- CreateSeuratObject(counts = BC2, meta.data=patient2, project = "seurat",min.cells = 3)

BC3 <- fread(file = "Datasets/GSM4006646_BC3_gene_cell_exprs_table.txt")
BC3 <- BC3[,-1]
BC3=aggregate(.~Symbol,max,data=BC3)
BC3 <- BC3%>%
  tibble::column_to_rownames('Symbol')
patient3 <- data.frame(rep('BC3',3100))
rownames(patient3) <- colnames(BC3)
colnames(patient3) <- "patient"
BC3 <- CreateSeuratObject(counts = BC3, meta.data=patient3, project = "seurat",min.cells = 3)

BC4 <- fread(file = "Datasets/GSM4006647_BC4_gene_cell_exprs_table.txt")
BC4 <- BC4[,-1]
BC4=aggregate(.~Symbol,max,data=BC4)
BC4 <- BC4%>%
  tibble::column_to_rownames('Symbol')
patient4 <- data.frame(rep('BC4',5131))
rownames(patient4) <- colnames(BC4)
colnames(patient4) <- "patient"
BC4 <- CreateSeuratObject(counts = BC4, meta.data=patient4, project = "seurat",min.cells = 3)

BC5 <- fread(file = "Datasets/GSM4006648_BC5_gene_cell_exprs_table.txt")
BC5 <- BC5[,-1]
BC5=aggregate(.~Symbol,max,data=BC5)
BC5 <- BC5%>%
  tibble::column_to_rownames('Symbol')
patient5 <- data.frame(rep('BC5',8643))
rownames(patient5) <- colnames(BC5)
colnames(patient5) <- "patient"
BC5 <- CreateSeuratObject(counts = BC5, meta.data=patient5, project = "seurat",min.cells = 3)

BC6 <- fread(file = 'Datasets/GSM4751267_BC6_gene_cell_exprs_table.txt')
BC6 <- BC6[,-1]
BC6=aggregate(.~Symbol,max,data=BC6)
BC6 <- BC6%>%
  tibble::column_to_rownames('Symbol')
patient6 <- data.frame(rep('BC6',6167))
rownames(patient6) <- colnames(BC6)
colnames(patient6) <- "patient"
BC6 <- CreateSeuratObject(counts = BC6, meta.data=patient6, project = "seurat",min.cells = 3)

BC7 <- fread(file = 'Datasets/GSM4751268_BC7_gene_cell_exprs_table.txt')
BC7 <- BC7[,-1]
BC7=aggregate(.~Symbol,max,data=BC7)
BC7 <- BC7%>%
  tibble::column_to_rownames('Symbol')
patient7 <- data.frame(rep('BC7',5523))
rownames(patient7) <- colnames(BC7)
colnames(patient7) <- "patient"
BC7 <- CreateSeuratObject(counts = BC7, meta.data=patient7, project = "seurat",min.cells = 3)

BC <- merge(BC1,y = c(BC2,BC3,BC4,BC5,BC6,BC7), merge.data = TRUE)
save(BC,file = 'Results/1.1 Single-cell data processing/BC.rda')

#细胞聚类--------------------------------------------------------
rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(SingleCellExperiment)

load('Results/1.1 Single-cell data processing/BC.rda')
# mat <- as.matrix(BC@assays$RNA@data)
# df <- as.data.frame(mat)
# write.csv(df, file = "Results/1.1 Single-cell data processing/Single-cell RNA-seq data.csv", row.names = FALSE)

#质控
BC[["percent.mt"]] <- PercentageFeatureSet(BC, pattern = "^MT-")
## 使用小提琴图可视化QC指标
VlnPlot(BC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##细胞过滤
BC <- subset(BC, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(BC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#归一化/标准化，作用是让测序数据量不同的细胞的基因表达量具有可比性
BC <- NormalizeData(BC, normalization.method = "LogNormalize", scale.factor = 10000)

#查找高变基因
BC <- FindVariableFeatures(BC, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(BC), 10)
plot1 <- VariableFeaturePlot(BC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#PCA 
#Scaling the data,对数据做了标准化和中心化，使数据成高斯分布，为PCA做准备
all.genes <- rownames(BC)
BC <- ScaleData(BC, features = all.genes)
BC <- RunPCA(BC, features = VariableFeatures(object = BC))
# 画图 
PCAPlot(BC,split.by = "patient") 

immune.combined <- RunUMAP(BC, reduction = "pca", dims = 1:30)
DimPlot(immune.combined, reduction = "umap", group.by = "patient")
print(BC[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(BC, dims = 1:2, reduction = "pca")

###JackStraw为聚类做伏笔，确定设置多少个dim来聚类合适
BC <- JackStraw(BC)
BC <- ScoreJackStraw(BC, dims = 1:20)
JackStrawPlot(BC, dims = 1:20)
#弯道图（Elbow plot）法,选择合适的PC数，一般选择拐点作为降维的维度
ElbowPlot(BC,ndims = 50) 

#选择最佳聚类数
BC <- FindNeighbors(BC, reduction = "pca", dims = 1:30)
BC <- FindClusters(BC,resolution = 1) #resolution可调细化聚类
head(Idents(BC), 5) #Idents获取对象的标识类

#UMAP聚类
BC <- RunUMAP(BC, dims = 1:30)
DimPlot(BC, reduction = "umap",label = T) + NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:31]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:31]))+  #对边框色
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
ggsave(filename = 'Results/1.1 Single-cell data processing/True_umap.pdf',width = 8,height = 6)

DimPlot(BC,group.by = "patient",reduction = "umap",label = T) + NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:7]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:7]))+  #对边框色
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = 'Patient')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
ggsave(filename = 'Results/1.1 Single-cell data processing/umap_patient.pdf',width = 8,height = 6)

#TSNE聚类
BC <- RunTSNE(BC, dims = 1:30)
DimPlot(BC, reduction = "tsne",label = T)+ NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:31]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:31]))+  #对边框色
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
ggsave(filename = 'Results/1.1 Single-cell data processing/True_tsne.pdf',width = 8,height = 6)

DimPlot(BC,group.by = "patient",reduction = "tsne",label = T) + NoLegend()+
  scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:7]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(18)[1:7]))+  #对边框色
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = 'Patient')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
ggsave(filename = 'Results/1.1 Single-cell data processing/tsne_patient.pdf',width = 8,height = 6)

save(BC,file = 'Results/1.1 Single-cell data processing/post-clustering_BC.rda')

#查看细胞亚群的差异基因--------------------------------------------------------
load('Results/1.1 Single-cell data processing/聚类后BC.rda')
BC_markers <- FindAllMarkers(BC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
save(BC_markers,file = 'Results/1.1 Single-cell data processing/AllMarkers.rda')

#Marker_Vlnplot-----
load('Results/1.1 Single-cell data processing/AllMarkers.rda')
top <- BC_markers %>% group_by(cluster) %>% slice_max(n = 3,order_by=avg_log2FC)
VlnPlot(BC, features = top$gene,  
        stack=T,pt.size=0,  
        #cols = my36colors,
        combine = "horizontal")+   
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())+
  guides(fill = FALSE)  # 隐藏图例
ggsave(filename = 'Results/1.1 Single-cell data processing/Marker_Vlnplot_top3.pdf', width = 15, height = 10)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置
genes <- c("EPCAM", "KRT8",'KRT18','CD14','CSF1R','AIF1','CD2',"CD3D",'CD3E',
           'DCN','PDPN','TAGLN','PECAM1','VWF','CLDN5')
VlnPlot(BC, features = genes,  
        stack=T,pt.size=0,  
        cols = my36colors,
        combine = T)+   #"horizontal"
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())+
  guides(fill = FALSE)  # 隐藏图例
ggsave(filename = 'Results/1.1 Single-cell data processing/Marker_Vlnplot.pdf', width = 6, height = 10)


#Marker气泡图
DotPlot(BC, features = unique(top$gene))+coord_flip()+theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave(filename = 'Results/1.1 Single-cell data processing/Marker_DotPlot.pdf', width = 15, height = 30)

FeaturePlot(BC, features = c("EPCAM", "KRT8",'KRT18'))
ggsave(filename = 'Results/1.1 Single-cell data processing/Marker_FeaturePlot.pdf', width = 10, height = 8)

#Marker_DoHeatmap
top5 <- BC_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(BC, features = top5$gene) + NoLegend()
ggsave(filename = 'Results/1.1 Single-cell data processing/Marker_DoHeatmap.pdf', width = 30, height = 15)




#细胞注释--------------------------------------------------
new.cluster.ids <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Epithelial","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Epithelial","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Epithelial","Epithelial",
                     "Myeloid/Macrophage","Epithelial","Epithelial","T cells","Fibroblast",
                     "Fibroblast","Epithelial","Endothelial","Epithelial","Epithelial",
                     "T cells")
names(new.cluster.ids) <- levels(BC)
BC <- RenameIdents(BC, new.cluster.ids)
save(BC,file = "Results/1.1 Single-cell data processing/post-annotated_BC.rda")

#注释后绘图
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')
DimPlot(BC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+
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
ggsave(filename = 'Results/1.1 Single-cell data processing/post-annotated_UMAP.pdf',width = 9,height = 6)

DimPlot(BC, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()+
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
ggsave(filename = 'Results/1.1 Single-cell data processing/post-annotated_TSNE.pdf',width = 8,height = 6)



#绘制各个样本中的细胞类型相对含量-------------------------------------
library(data.table)
library(tibble)
library(Biobase)
library(SeuratObject)
library(sp)
library(Seurat)
library(tidyr)
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')
active_ident <- BC@active.ident
# 将active.ident信息添加到meta.data中
BC <- AddMetaData(BC, metadata = active_ident, col.name = "active.ident")
metadata <- BC@meta.data
metadata <- metadata[,c(4,8)]
metadata <- rownames_to_column(metadata,'ID')
colnames(metadata)[3] <- "Celltype"

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(agricolae)
p1 <- ggplot(metadata,aes(patient,fill=Celltype))+geom_bar(position="fill")+#coord_flip()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + #scale_fill_npg()
  scale_fill_manual(values = c('#E41A1C','#377EB8','#4DAF4A','#984EA3',"#F29403")) 
p1 + theme_bw() + labs(y="Proportion") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = 'Results/1.1 Single-cell data processing/Celltype_percent_barplot.pdf',width = 4,height = 5)


#Cell cycle score---------------------------
library(data.table)
library(tibble)
library(Biobase)
library(SeuratObject)
library(sp)
library(Seurat)
library(tidyr)
library(ggplot2)
library(ggsci)
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')

cc.genes
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(BC)) #提pbmc拒阵中的G2M期marker基医
s_genes <- cc.genes$s.genes #获取S期marker基因
s_genes <- CaseMatch(search=s_genes,match=rownames(BC)) #提取pbmc拒阵中的S期marker基因
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对pbmc进行细胞周期评分
BC <- CellCycleScoring(BC,g2m.features = g2m_genes,s.features = s_genes)
colnames(BC@meta.data)
table(BC$Phase)

DimPlot(BC,group.by = "Phase",reduction = "umap",label = T) + NoLegend()+
  scale_fill_manual(values = c('#377EB8','#000000','#E41A1C'))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
  scale_color_manual(values = c('#377EB8','#000000','#E41A1C'))+  #对边框色
  theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
  labs(title = 'Cell Cycle Score')+  #labs设置图形标签
  theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
        axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
        axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
        axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
        panel.grid = element_blank(),  #网格线。element_blank()清空
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
        panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
        #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
ggsave(filename = 'Results/1.1 Single-cell data processing/umap_Phase.pdf',width = 8,height = 6)



#提取上皮细胞类群,再聚类---------------------------
#rm(list = ls())
# library(data.table)
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(Seurat)
# library(patchwork)
# library(ggplot2)
# library(ggsci)
# BC_Epithelial <- subset(BC,idents ="Epithelial")
# 
# 
# sce=BC_Epithelial
# sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000) 
# sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 3000)
# all.genes <- rownames(sce)
# sce <- ScaleData(sce, features = all.genes)
# sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
# 
# elbow_plot <- ElbowPlot(sce,ndims = 50)
# ggsave("1.单细胞数据处理/Epithelial_elbow_plot.pdf", width = 10, height = 10,plot = elbow_plot)
# 
# sce <- FindNeighbors(sce, dims = 1:22)
# sce <- FindClusters(sce, resolution = 1)  #, resolution = 1 
# # Look at cluster IDs of the first 5 cells
# head(Idents(sce), 5)
# table(sce$seurat_clusters) 
# 
# sce <- RunUMAP(sce, dims = 1:22)
# DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',
#         label = TRUE, pt.size = 0.5) + NoLegend()+
#   scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:25]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
#   scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:25]))+  #对边框色
#   theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
#   labs(title = '')+  #labs设置图形标签
#   theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
#         axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
#         axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
#         axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
#         panel.grid = element_blank(),  #网格线。element_blank()清空
#         panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
#         panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
#         #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
#         plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题。 seurat_clusters 
# ggsave(filename = '1.单细胞数据处理/Epithelial_UMAP.pdf',width = 8,height = 6)
# 
# DimPlot(sce, reduction = 'umap', group.by = 'patient',
#         label = TRUE, pt.size = 0.5) + NoLegend()+
#   scale_fill_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:25]))+  #在ggsci包中，对图形的填充色进行美化。pal_npg选择调色板类型，alpha选择透明度
#   scale_color_manual(values = c(pal_d3('category20',alpha = 0.9)(20)[1:25]))+  #对边框色
#   theme_bw(base_rect_size = 1)+  #base_rect_size设置边框线条的粗细
#   labs(title = '')+  #labs设置图形标签
#   theme(axis.text.y = element_text(size = 12, colour = 'black'),   #theme自定义图形中的非数据组件。对y轴上的文本进行编辑。
#         axis.text.x = element_text(size = 12, colour = 'black'),   #对X轴上的文本进行编辑（字体，颜色，字体面，角度，横向调整）
#         axis.title.y = element_text(size = 15, colour = 'darkred',face = 'bold'),  #坐标轴y的名称。element_text对Y轴的文本进行编辑（字体，颜色，字体面）
#         axis.title.x = element_text(size = 15, colour = 'darkred',face = 'bold'),   #坐标轴x的名称。element_blank()不绘制任何内容，也不分配任何空间
#         panel.grid = element_blank(),  #网格线。element_blank()清空
#         panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),  #主网格线。定义网格线颜色、线的类型
#         panel.background = element_rect(fill='#f3f6f6'),  #背景。fill填充颜色
#         #legend.position = 'none',   #图例的位置("none"， "left"， "right"， "bottom"， "top")
#         plot.title = element_text(hjust = 0.5,size = 18,colour = 'darkred',face='bold')) #图标题 
# ggsave(filename = '1.单细胞数据处理/Epithelial_Sample_UMAP.pdf',width = 8,height = 6)
# 

# #基因表达气泡图
# genes_to_check = c("EPCAM", "KRT8",'KRT18')
# DotPlot(sce, group.by = 'seurat_clusters',
#         features = unique(genes_to_check)) + RotatedAxis()
# 
# # p1=DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',
# #            label = TRUE, pt.size = 0.5) + NoLegend()
# # p2=DotPlot(sce, group.by = 'seurat_clusters',
# #            features = unique(genes_to_check)) + RotatedAxis()
# # library(patchwork)
# # p1+p2
# save(sce,file = '1.单细胞数据处理/Epithelial_cluster.Rdata')





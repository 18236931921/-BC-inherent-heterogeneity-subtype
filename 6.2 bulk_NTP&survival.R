rm(list=ls())
library(data.table)
library(tidyverse)
library(ggsci)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
#library(rJava)
#library(xlsx)
library(CMScaller)
library(survival)
library(survminer)

#单细胞中注释细胞分群信息---------------------------------------
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

#单细胞数据集与bulk数据集的基因取交集------------------------------
load('Datasets/Bladder cancer.rda')
library(clusterProfiler)
library(org.Hs.eg.db)

new_total_expr_list <- total_expr_list[!names(total_expr_list) %in% c("GSE48075", "GSE52219")]
new_list <- list()
for (i in names(new_total_expr_list)) {
  expr <- total_expr_list[[i]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
new_list[[i]] <- expr
}
expr <- total_expr_list[["GSE48075"]]
new_list[["GSE48075"]] <- expr

interID <- Reduce(intersect,list(rownames(new_list[["GSE13507"]]),
                                 rownames(new_list[["GSE154261"]]),
                                 rownames(new_list[["GSE19423"]]),
                                 rownames(new_list[["GSE31684"]]),
                                 rownames(new_list[["GSE37815"]]),
                                 rownames(new_list[["GSE39281"]]),
                                 rownames(new_list[["GSE48276"]]),
                                 rownames(new_list[["GSE69795"]]),
                                 rownames(new_list[["GSE70691"]]),
                                 rownames(new_list[["GSE48075"]]),
                                 rownames(new_list[["IMvigor210"]]),
                                 rownames(new_list[["TCGA_BLCA"]])))
interID <- intersect(interID,rownames(ssBC_cells))
ssBC_cells <- subset(ssBC_cells, features = interID)

#提取细胞亚群的差异基因------------------------------------------------
#方法一FindAllMarkers（不合适）--------------
ssBC_markers <- FindAllMarkers(ssBC_cells,logfc.threshold = 0,min.pct = 0.1,only.pos = T,return.thresh = 0.01)
table(ssBC_markers$cluster)

markers <- subset(ssBC_markers,p_val_adj<0.05&avg_log2FC>0.7)  #p_val_adj<0.01&&avg_log2FC>0.3
table(markers$cluster)
template <- markers[,c(6,7)]
colnames(template) <- c('class','probe')

#方法二FindMarkers（不合适）------------------------------
# ssBC_cells <- SetIdent(ssBC_cells,value = ssBC_cells@meta.data$cluster)
cluster_levels <- as.data.frame(table(ssBC_cells@meta.data$Cluster))[,1]
all_markers <- lapply(cluster_levels, function(Cluster) {
  FindMarkers(ssBC_cells, ident.1 = Cluster)#, min.pct = 0.1 默认0.1
})
# 添加聚类身份到每个数据框中
all_markers_df <- lapply(seq_along(all_markers), function(i) {
  df <- all_markers[[i]]
  df$cluster <- cluster_levels[i]
  return(df)
})
# 合并所有的数据框
all_markers_df <- do.call(rbind, all_markers_df)
table(all_markers_df$cluster)

markers <- subset(all_markers_df,p_val_adj<0.05&avg_log2FC>0.5)  #&avg_log2FC>0.3
markers <- rownames_to_column(markers,'gene')
table(markers$cluster)
template <- markers[,c(1,7)]
template <- template[,c(2,1)]
colnames(template) <- c('class','probe')

#方法三wilcoxauc（不合适）-----------------------------
library(presto) 
genes <- wilcoxauc(ssBC_cells, 'Cluster')
dplyr::count(genes, group)#查看每个cluster中有多少基因（与矩阵的基因数一致)#
genes <- as.data.frame(genes)
range(genes$logFC)
range(genes$padj)
range(genes$auc)

scRNA_markers2 <- subset(genes,padj<0.05&logFC>0&auc>0.6)  #&logFC>0.35
# scRNA_markers2 <- scRNA_markers2 %>% group_by(group) %>% slice_max(n = 200,order_by=logFC)
table(scRNA_markers2$group)
template <- scRNA_markers2[,c(2,1)]
colnames(template) <- c('class','probe')

#方法四COSG（选用此方法）-----------------------------
library(COSG)
dim(ssBC_cells@assays$RNA@data)  #标准化矩阵
table(Idents(ssBC_cells)) #聚类分群结果

marker_cosg <- cosg(
  ssBC_cells,
  groups='all', #考虑全部分组
  assay='RNA',
  slot='data',
  mu=1,         #惩罚项参数，值越大
  remove_lowly_expressed=TRUE,   #是否过滤低表达基因
  expressed_pct=0.1,             #设置低表达的阈值
  n_genes_user=1000      #每个cluster定义Top-N个marker gene
)

top_list<-data.frame()
for (group in colnames(marker_cosg$names)){
  top_i<-marker_cosg$names[group][1:550,1]
  top_x <- marker_cosg$scores[group][1:550,1]
  #top_list<-c(top_list,top_i)
  top <- data.frame(gene=top_i,cluster=rep(group,550),scores=top_x)
  top_list <- rbind(top_list,top)
}
template <- top_list[,c(2,1,3)]
colnames(template) <- c('class','probe',"scores")
table(template$class)

#为bisque做准备
markers <- template
markers$class <- gsub("BCIS", "C", markers$class)
save(markers,file = "Results/6. bulk_NTP&survival/markers.rda")

#bulk中执行NTP-------------------------------------------------------------
#TCGA_BCLA--------------------------
load("Results/6. bulk_NTP&survival/markers.rda")
template <- markers
load('Datasets/Bladder cancer.rda')
emat <- new_list[["TCGA_BLCA"]][,substr(colnames(new_list[["TCGA_BLCA"]]),14,16) == '01A']
TCGA_clin <- total_clin_list[["TCGA_BLCA"]]%>%
  tibble::column_to_rownames("ID")
clin <- TCGA_clin[colnames(emat),]
range(emat)
emat <- as.data.frame(emat)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by.x = 1,by.y = 0)
table(dd$prediction)

dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

#保存结果进行后续分析
data <- emat[,dd$ID]
save(dd,data,file = "Results/6. bulk_NTP&survival/TCGA_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('TCGA_BLCA (n=289)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/TCGA_OS_survival.pdf',width = 8,height = 8)


#PFS
dd$PFS <- as.numeric(dd$PFS)
dd$PFS.time <- as.numeric(dd$PFS.time)
fit <- survfit(Surv(PFS.time,PFS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)
sdf <- survdiff(Surv(PFS.time,PFS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = PFS,time = PFS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Progression-free survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('TCGA_BLCA (n=289)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/TCGA_PFS_survival.pdf',width = 8,height = 8)

#DSS
dd$DSS <- as.numeric(dd$DSS)
dd$DSS.time <- as.numeric(dd$DSS.time)
fit <- survfit(Surv(DSS.time,DSS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)
sdf <- survdiff(Surv(DSS.time,DSS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = DSS,time = DSS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Disease-special survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('TCGA_BLCA (n=289)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/TCGA_DSS_survival.pdf',width = 8,height = 8)


#IMvigor210-----------------------------------------
expr <- new_list[["IMvigor210"]]
clin <- total_clin_list[["IMvigor210"]]
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/IMvigor210_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)

dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

#保存结果进行后续分析
IMvigor210_data <- emat[,dd$ID]
IMvigor210_dd <- dd
save(IMvigor210_dd,IMvigor210_data,file = "Results/6. bulk_NTP&surviva/IMvigor210_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+ 
  #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('IMvigor210 (n=253)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/IMvigor210_survival.pdf',width = 8,height = 8)



#Atezolizumab免疫治疗反应-箱线图
load('Results/6. bulk_NTP&surviva/IMvigor210_cluster.rda')
mypalette <- c(alpha("#003153",0.8),alpha('#008080',0.8),"#c0c0c0",alpha("darkred",0.8))
plotdata <- IMvigor210_dd[,c(1,2,7)]
x <- aggregate(plotdata$prediction,by = list(Cluster = plotdata$prediction,
                                           Response = plotdata$Atezolizumab),length)
data2 <- x
p <- fisher.test(matrix(as.numeric(data2$x),nrow = 2))$p.value
ct12 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C2'),]$x),nrow = 2))$p.value
ct23 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C3'),]$x),nrow = 2))$p.value
ct34 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C3','C4'),]$x),nrow = 2))$p.value
ct14 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C4'),]$x),nrow = 2))$p.value
ct13 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C3'),]$x),nrow = 2))$p.value
ct24 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C4'),]$x),nrow = 2))$p.value
data2$percent <- c(42/70,39/55,21/24,60/70,28/70,16/55,3/24,10/70)
data2$label <- round(data2$percent*100,1)

ggplot(data2, aes( x = Cluster, y = percent, fill = Response))+
  geom_bar( stat = 'identity' ,position = position_stack())+
  geom_text(aes(label = label),position = position_stack(vjust = 0.5),size=5,color="black")+
  labs(title = 'Response')+
  ylab(label = 'Percent')+xlab(label = NULL)+
  #scale_y_continuous(expand = c(0,0.05))+
  scale_y_continuous(limits = c(-0.05, 1.15), breaks = c(0,0.25,0.5,0.75,1.0))+
  # ylim(0.0,1.05)+
  theme_bw(base_rect_size = 2)+
  theme(axis.title.y = element_text(size = 18,colour = 'black', face = 'bold'),
        axis.title.x = element_text(size = 18,colour = 'black', face = 'bold',vjust = -0.5),
        axis.text.x = element_text(size = 16,colour = 'black',face = 'bold',hjust =0.5),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1.5),
        panel.background = element_rect(fill = '#f3f6f6', color = NA),
        # panel.grid = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.grid.minor.x = element_line(color = '#f3f6f6', linetype = "dashed"),
        legend.text = element_text(size = 12,colour = 'black'),
        legend.title = element_text(size = 14,colour = 'black', face = 'bold'),
        legend.position = 'bottom',
        legend.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(size = 18,colour = 'black', face = 'bold',hjust = 0.5,vjust = 1.5))+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  annotate("text",x = c(1.5),y = -0.05,label = c("ns"),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(2.5),y = -0.05,label = c('ns'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3.5),y = -0.05,label = c('ns'),size =6,fontface = "bold")+
  annotate("text",x = c(2),y = 1.03,label = c('**'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3),y = 1.09,label = c('*'),size =6,fontface = "bold")+
  annotate("text",x = c(2.5),y = 1.15,label = c('***'),size =6,fontface = "bold")+
  annotate("segment",x = 1.2,xend = 1.8, y = -0.03,yend = -0.03,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.06,'cm')))+
  annotate("segment",x = 2.2,xend = 2.8, y = -0.03,yend = -0.03,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.06,'cm')))+
  annotate("segment",x = 3.2,xend = 3.8, y =-0.03,yend = -0.03,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.06,'cm')))+
  annotate("segment",x = 1,xend = 3, y = 1.02,yend = 1.02,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.08,'cm')))+
  annotate("segment",x = 2,xend = 4, y = 1.08,yend = 1.08,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.08,'cm')))+
  annotate("segment",x = 1,xend = 4, y = 1.14,yend = 1.14,size = 0.8,
           arrow = arrow(ends = 'both',angle = 90,length = unit(0.08,'cm')))
ggsave(filename = 'Results/6. bulk_NTP&surviva/Response.pdf',width = 4.7,height = 9)


#旧图画法
data <- dd[,c(1,2,7)]
Atezo <- column_to_rownames(data,var = "ID")
colnames(Atezo)[1] <- 'Cluster'
Atezo1 <- Atezo[Atezo$Cluster == 'C1',]
Atezo2 <- Atezo[Atezo$Cluster == 'C2',]
Atezo3 <- Atezo[Atezo$Cluster == 'C3',]
Atezo4 <- Atezo[Atezo$Cluster == 'C4',]

table(Atezo1$Atezo)
table(Atezo2$Atezo)
table(Atezo3$Atezo)
table(Atezo4$Atezo)

d<- data.frame(Var1 =c(rep('C1',2),rep('C2',2),rep('C3',2),rep('C4',2)),kind=c('NR','R'),amunt=c(42,28,39,16,21,3,60,10))
data <- pivot_wider(d,values_from = 'amunt',names_from = 'Var1') %>% 
  column_to_rownames('kind')
ct <- chisq.test(data)
pvalue <- ct$p.value

ggplot(d, aes( x = Var1, y = amunt, fill = kind))+
  geom_bar( stat = 'identity' ,position = "fill")+
  ylab(label = 'Percent')+xlab(label = '')+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme_bw(base_rect_size = 2)+
  theme(axis.title.y = element_text(size = 18,colour = 'black', face = 'bold'),
        axis.title.x = element_text(size = 18,colour = 'black', face = 'bold',vjust = -0.5),
        axis.text.x = element_text(size = 14,colour = 'black',vjust = -0.5),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1.5),
        panel.background = element_rect(fill = '#f3f6f6', color = NA),
        # panel.grid = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.grid.minor.x = element_line(color = '#f3f6f6', linetype = "dashed"),
        legend.text = element_text(size = 12,colour = 'black'),
        # legend.title = element_blank(),
        legend.position = 'right',
        legend.background = element_rect(fill = NA, color = NA))+ 
  labs(fill = "Atezo")+ # 设置图例标题
  scale_fill_manual(values = c('#e07a5f','#118ab2')) #,'#81b29a','#f2cc8f','#8d85c1','#c1a180','#ff0000','#00ced1','#ffa07a
ggsave(filename = "Results/6. bulk_NTP&surviva/Atezo_barplot.pdf", width =5.2, height = 7)



#E-MTAB-4321--------------------------
load('Results/6. bulk_NTP&surviva/E-MTAB-4321.rda')
emat <- expr
emat <- t(scale(t(emat)))
emat <- na.omit(emat)
range(emat)
emat <- as.data.frame(emat)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/E-MTAB-4321_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by.x = 1,by.y = 0)
table(dd$prediction)
dd$PFS <- as.numeric(dd$PFS)
dd$PFS.time <- as.numeric(dd$PFS.time)
fit <- survfit(Surv(PFS.time,PFS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

#保存结果进行后续分析
EMTAB_data <- emat[,dd$ID]
EMTAB_dd <- dd
save(EMTAB_dd,EMTAB_data,file = "Results/6. bulk_NTP&surviva/EMTAB_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(PFS.time,PFS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = PFS,time = PFS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) + geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Progression-Free survival',x='Time (years)')+  
  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.9),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.8,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('E-MTAB-4321 (n=340)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/E-MTAB-4321_survival.pdf',width = 8,height = 8)




#GSE13507-----------------------
expr <- total_expr_list[["GSE13507"]]
clin <- total_clin_list[["GSE13507"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE13507_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

#保存结果进行后续分析
GSE13507_data <- emat[,dd$ID]
GSE13507_dd <- dd
save(GSE13507_dd,GSE13507_data,file = "Results/6. bulk_NTP&surviva/GSE13507_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE13507 (n=124)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE13507_survival.pdf',width = 8,height = 8)


#GSE70691-----------------------------------------
expr <- total_expr_list[["GSE70691"]]
clin <- total_clin_list[["GSE70691"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE70691_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

#保存结果进行后续分析
GSE70691_data <- emat[,dd$ID]
GSE70691_dd <- dd
save(GSE70691_dd,GSE70691_data,file = "Results/6. bulk_NTP&surviva/GSE70691_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE70691 (n=31)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE70691_survival.pdf',width = 8,height = 8)



#GSE154261-----------------------
expr <- total_expr_list[["GSE154261"]]
clin <- total_clin_list[["GSE154261"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$PFS <- as.numeric(dd$PFS)
dd$PFS.time <- as.numeric(dd$PFS.time)
fit <- survfit(Surv(PFS.time,PFS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE154261_survival.pdf',width = 8,height = 8)

#GSE19423 -----------------------
expr <- total_expr_list[["GSE19423"]]
clin <- total_clin_list[["GSE19423"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE19423_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE19423 (n=36)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE19423_survival.pdf',width = 8,height = 8)

#GSE31684可试-----------------------------------------
expr <- total_expr_list[["GSE31684"]]
clin <- total_clin_list[["GSE31684"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE31684_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE31684 (n=61)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE31684_survival.pdf',width = 8,height = 8)

#GSE48075 可试-----------------------------------------
expr <- total_expr_list[["GSE48075"]]
clin <- total_clin_list[["GSE48075"]]
# ID <- rownames(expr)
# ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
# expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
# expr <- expr%>%
#   group_by(SYMBOL)%>%
#   summarise_all(max)%>%
#   filter(SYMBOL!='')%>%
#   column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE48075_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE48075 (n=47)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE48075_survival.pdf',width = 8,height = 8)

#GSE48276 可试-----------------------------------------
expr <- total_expr_list[["GSE48276"]]
clin <- total_clin_list[["GSE48276"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/GSE48276_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('GSE48276 (n=45)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE48276_survival.pdf',width = 8,height = 8)

#GSE37815-----------------------------------------
expr <- total_expr_list[["GSE37815"]]
clin <- total_clin_list[["GSE37815"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE37815_survival.pdf',width = 8,height = 8)

#GSE39281 无预测出的样本-----------------------------------------
expr <- total_expr_list[["GSE39281"]]
clin <- total_clin_list[["GSE39281"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE39281_survival.pdf',width = 8,height = 8)

#GSE69795 -----------------------------------------
expr <- total_expr_list[["GSE69795"]]
clin <- total_clin_list[["GSE69795"]]
ID <- rownames(expr)
ID <- bitr(ID,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
expr <- merge(ID,expr,by.x=1,by.y=0)[,-1]
expr <- expr%>%
  group_by(SYMBOL)%>%
  summarise_all(max)%>%
  filter(SYMBOL!='')%>%
  column_to_rownames('SYMBOL')
emat <- as.data.frame(expr)

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
res <- ntp_res[ntp_res$FDR<0.2,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by= 1)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)

ggsave(filename = 'Results/6. bulk_NTP&surviva/GSE69795_survival.pdf',width = 8,height = 8)


#Meat队列--------------------------
load('Results/6. bulk_NTP&surviva/markers.rda')
load('Results/6. bulk_NTP&surviva/4Metadata.Rda')

emat <- dd
clin <- cc
range(emat)
emat <- as.data.frame(emat)
template <- markers

ntp_res <- ntp(emat      = emat,
               templates = template,
               doPlot    = T,
               nPerm     = 1000,
               distance  = "cosine",
               nCores    = 10,
               seed      = 2020104,
               verbose   = T)
#绘制ntp热图
pdf(file="Results/6. bulk_NTP&surviva/Meat_subheatmap.pdf",width = 10,height = 9)
subHeatmap(emat = emat, res = ntp_res, templates = template,
           classCol = c("#436b95","#b81d25",'#2f7e2f','#800080'),heatCol = c('seagreen','orange'))
dev.off()

res <- ntp_res[ntp_res$FDR<0.01,]%>%tibble::rownames_to_column('ID')

dd <- merge(res[,1:2],clin,by.x = 1,by.y = 0)
table(dd$prediction)
dd$OS <- as.numeric(dd$OS)
dd$OS.time <- as.numeric(dd$OS.time)
fit <- survfit(Surv(OS.time,OS)~prediction,dd)
ggsurvplot(fit,pval = T,risk.table = T)


#保存结果进行后续分析
Meta_data <- emat[,dd$ID]
Meta_dd <- dd
save(Meta_dd,Meta_data,file = "Results/6. bulk_NTP&surviva/Metadata_cluster.rda")

#优化绘图
#测试生存曲线差异并定义P值
sdf <- survdiff(Surv(OS.time,OS)~prediction,data = dd) #生存曲线差异
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)   #卡方检验中的pchisq给出分布函数
if(p.val<0.0001){
  p <- 'P < 0.001'
}else{
  p <- paste0('P = ',round(p.val,3))
}
p
library(ggplot2)
library(ggkm)
ggplot(dd, aes( status = OS,time = OS.time,color = prediction)) +      #aes构建美学绘图
  geom_km(size=1.2) +geom_kmticks()+    #geom_km生存率曲线制法，geom_kmticks在Kaplan Meier曲线上显示标记
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+  #scale_color_manual给图形调色
  #设置图形在x,Y轴上的标签
  labs(y='Overall survival',x='Time (years)')+  
  #修改主题的组件(标题、标签、字体、背景、网格线和图例)
  theme_bw(base_rect_size = 1.5)+    #base_rect_size设置边框线条的粗细
  theme(legend.position = c(0.82,0.87),   #图例的位置c(0.82,0.84)，"none"
        legend.text = element_text(size = 12,colour = 'black'),   #图例上的文本
        legend.title = element_blank(),   #图例的标题
        legend.key = element_rect(fill = "#f5f3eb", color = NA),  #图例下面的基调色
        legend.background = element_rect(fill = "#f5f3eb", color = NA),  #图例背景的颜色
        axis.text = element_text(size = 12,colour = 'black'),  #x,y轴上的文本字体设置
        axis.title = element_text(face = "bold",colour = "darkred",size = 18),  #x,y轴上的标题字体设置
        panel.background = element_rect(fill = "#f5f3eb", color = NA),  #背景。fill填充颜色
        panel.grid = element_blank(),  #设置网格线
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "#cacfd2", linetype = "dashed"),
        axis.line = element_line(color = "#606F7B"),
        plot.title = element_text(face = "bold",colour = "darkred",size = 18,hjust = 0.5))+   #设置图的标题
  annotate('text',x = 0.25,y=0.25,label=paste0('Log-rank\n',p),hjust=0,size=5,fontface = 'italic')+  #创建注释层
  ggtitle('Meta Chort (n=154)')
ggsave(filename = 'Results/6. bulk_NTP&surviva/4MetaChort_survival.pdf',width = 8,height = 8)



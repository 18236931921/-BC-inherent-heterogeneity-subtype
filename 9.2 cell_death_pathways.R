rm(list = ls())
library(tidyverse)
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
library(pheatmap)

#表达数据整理
gene <- read.table("Results/9. Function/cell_death_pathways/input_genes.txt",sep="\t",header=F,check.names=F)
load("Results/6. bulk_NTP&survival/TCGA_cluster.rda")
expr <- data[rownames(data)%in%gene$V1,]#采取细胞死亡基因

#自定义函数---------------------------------------------------------------------
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
} 

# 计算组间统计差异
cross_subtype_compr <- function(expr = NULL,
                                subt = NULL,
                                subt.label = "Subtype",
                                two_sam_compr_method = "wilcox",
                                multi_sam_compr_method = "kruskal",
                                res.path = NULL) {
  
  if (!is.element(two_sam_compr_method, c("t.test", "wilcox"))) {stop("Two samples comparison should be t.test or wilcox!\n") }
  if (!is.element(multi_sam_compr_method, c("anova", "kruskal"))) {stop("multiple samples comparison should be kruskal or anova!\n") }
  
  subt.name <- unique(subt[,subt.label])
  n.subt <- length(subt.name)
  if(n.subt < 2) {stop("The number of subtype should be greater than 2!\n")}
  
  comprTab <- NULL
  
  # 两个亚型且为非参数检验
  if(n.subt == 2 & two_sam_compr_method == "wilcox") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      wt <- wilcox.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = wt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 两个亚型且为参数检验
  if(n.subt == 2 & two_sam_compr_method == "t.test") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp1 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[1]),,drop = F])])
      tmp2 <- as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[2]),,drop = F])])
      tt <- t.test(tmp1,tmp2)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = tt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为非参数检验
  if(n.subt > 2 & multi_sam_compr_method == "kruskal") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      kt <- kruskal.test(value ~ subt,data = tmp)
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = kt$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 多个亚型且为参数检验
  if(n.subt > 2 & multi_sam_compr_method == "anova") {
    for (i in 1:nrow(expr)) {
      display.progress(index = i,totalN = nrow(expr))
      tmp.list <- list()
      for (n in 1:n.subt) {
        tmp.list[[n]] <- data.frame(value = as.numeric(expr[i,rownames(subt[which(subt[,subt.label] == subt.name[n]),,drop = F])]),
                                    subt = subt.name[n],
                                    stringsAsFactors = F)
      }
      tmp <- do.call(rbind,tmp.list)
      at <- summary(aov(value ~ subt,data = tmp))
      comprTab <- rbind.data.frame(comprTab,
                                   data.frame(gene = rownames(expr)[i],
                                              nominal.p.value = at[[1]][1,5],
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
    }
  }
  
  # 调整p值
  comprTab$adjusted.p.value = p.adjust(comprTab$nominal.p.value,method = "BH")
  # 按p值排序
  #comprTab <- comprTab[order(comprTab$adjusted.p.value, decreasing = F),] 
  
  write.table(comprTab,file.path(res.path,"comprTab.txt"),sep = "\t",row.names = F,quote = F)
  return(comprTab)
}


## 读取感兴趣的基因表达矩阵
mygene_data <- expr

## 读取分组信息
colnames(dd)[2] <- "Subtype"
Subtype <- dd[,c(1,2)]%>%tibble::column_to_rownames('ID')

# 查看各组名字和sample数量，画图时要用
table(Subtype$Subtype)




#绘制各种细胞死亡通路热图-------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(gplots) #用于丰富热图颜色
library(stringr)
library(pheatmap)
library(grid)
# library(iCluster)
library(iClusterPlus)
library(dplyr)
library(circlize)
library(ggsci)

# GSVA--------------
#表达数据整理
load("Results/6. bulk_NTP&survival/TCGA_cluster.rda")
exp <- data
sub <- dd[order(dd$prediction),c(1,2)]
exp <- exp[,sub$ID]
biocarta <- qusage::read.gmt("Results/9. Function/cell_death_pathways/c2.cp.biocarta.v7.5.1.symbols.gmt")
go <- qusage::read.gmt("Results/9. Function/cell_death_pathways/c5.go.bp.v7.5.1.symbols.gmt")
kegg <- qusage::read.gmt("Results/9. Function/cell_death_pathways/c2.cp.kegg.v7.5.1.symbols.gmt")
pid <- qusage::read.gmt("Results/9. Function/cell_death_pathways/c2.cp.pid.v7.5.1.symbols.gmt")
reactome <- qusage::read.gmt("Results/9. Function/cell_death_pathways/c2.cp.reactome.v7.5.1.symbols.gmt")
hall <- qusage::read.gmt("Results/9. Function/cell_death_pathways/h.all.v7.5.1.symbols.gmt")

#gsva
gsva_es1 <- gsva(as.matrix(exp), reactome)
gsva_es2 <- gsva(as.matrix(exp), kegg)
gsva_es3 <- gsva(as.matrix(exp), pid)
gsva_es4 <- gsva(as.matrix(exp), hall)
gsva_es5 <- gsva(as.matrix(exp), biocarta)
gsva_es6 <- gsva(as.matrix(exp), go)
gsva_cell_death <- rbind(gsva_es1,gsva_es2,gsva_es3,gsva_es4,gsva_es5,gsva_es6)
save(gsva_cell_death,file="Results/9. Function/cell_death_pathways/gsva_cell_death.rda")
# 把通路的表达量保存到文件
write.csv(gsva_cell_death, "Results/9. Function/cell_death_pathways/gsva_output.csv", quote = F)

# # #铁死亡
# ferroptosis <- qusage::read.gmt("Results/9. Function/cell_death_pathways/WP_FERROPTOSIS.v7.5.1.gmt")
# Ferroptosis <- gsva(as.matrix(exp), ferroptosis)


#细胞死亡热图---------------------------------------------------------------------
load('Results/9. Function/cell_death_pathways/gsva_cell_death.rda')
gsva_cell_death <- gsva_cell_death[,sub$ID]
table(sub$prediction)

Apoptosis <- gsva_cell_death[c(2,4,98,1072,1597,1729,2007,8481,8803,8804,8805,9886),]
Necroptosis <- as.data.frame(t(gsva_cell_death[8528,]))
rownames(Necroptosis) <- "GOBP_NECROPTOTIC_SIGNALING_PATHWAY"
Autophagy <- gsva_cell_death[c(1485,1490,1725,2291,3745,3746,3747,4152,4153,
                               4154,4156,4533,5304,6089,9189,9190,9277,9430,9460),]
pyroptosis <- gsva_cell_death[c(919,7752),]
Necrosis <- gsva_cell_death[c(853,854,1564,3856,3857,5093,5281,8047,9204,9264),]
lysosome <- gsva_cell_death[c(2071,2419,4900,7621,7666,8325,8400,8620,9122,
                              9545,763),]
Top = HeatmapAnnotation(Cluster=sub$prediction,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                     title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        gap = unit(5.5, "mm"),
                        col=list(Cluster=c("C1"="#003C67CC","C2"="#b81d25","C3"='#2f7e2f',"C4"='#800080')),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
table(sub$prediction)
# range(Ferroptosis)
range(Apoptosis)
range(Autophagy)
range(pyroptosis)
range(Necrosis)
range(lysosome)
plotdata <- rbind(Apoptosis,Autophagy,Necrosis,lysosome,pyroptosis,Necroptosis)

range(plotdata)
pdf('Results/9. Function/cell_death_pathways/celldeath_heatmap2.pdf',height = 10,width = 16)
Heatmap(as.matrix(plotdata),##对基因表达量进行限定
        name='Z-score',#热图的名字
        top_annotation = Top,#顶部注释  
        cluster_rows = F,#不对行进行聚类
        col=colorRamp2(c(-0.3,0,0.3),
                       c(pal_npg("nrc", alpha =0.6)(2)[2],
                         'white',
                         pal_npg("nrc", alpha =0.6)(2)[1])),#颜色
        color_space = "RGB",
        #rect_gp = gpar(col = "white", lwd = 2),#热图内部矩形的参数，如设置举行边框为白色
        
        width = unit(8, "cm"),
        cluster_columns = F,
        
        border = T,
        border_gp = gpar(col = "black",lwd=2),
        
        row_order=NULL,
        column_order=NULL,
        column_names_rot = 45,
        row_names_side = 'left',
        show_column_names = F,
        show_row_names = T,  
        row_names_gp = gpar(#col = c(rep("red", 10), rep("blue", 8))
          fontsize = 6),
        
        
        column_split = c(rep("C1",119),rep("C2",74),rep("C3",40),rep("C4",56)),
        row_split = factor(rep(c('Apoptosis','Necroptosis','Autophagy','pyroptosis','Necrosis','lysosome'),
                               times=c(12,1,19,2,10,11)),
                           levels = c('Apoptosis','Necroptosis','Autophagy','pyroptosis','Necrosis','lysosome')),
        row_gap = unit(c(1.5), "mm"),
        column_gap = unit(c(1.5), "mm"),
        
        column_title = NULL,
        column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 5),
        row_title_gp = gpar(fontsize=12),
        
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), 
                                  title_gp = gpar(fontsize = 10, fontface = "bold"))) 

dev.off()







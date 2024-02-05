#GSE13507+GSE70691------------------------------------------
rm(list = ls())
library(data.table)
library(tibble)
library(dplyr)

load('Datasets/Bladder cancer.rda')
#对多个GEO表达数据的合并
GSE13507_expr <- total_expr_list[["GSE13507"]]
GSE13507_clin <- total_clin_list[["GSE13507"]]
GSE70691_expr <- total_expr_list[["GSE70691"]]
GSE70691_clin <- total_clin_list[["GSE70691"]]

gene <- intersect(rownames(GSE13507_expr),rownames(GSE70691_expr))
expr_GSE70691 <- GSE70691_expr[gene,]
expr_GSE13507 <- GSE13507_expr[gene,]

ee <- Reduce(cbind,list(expr_GSE13507,expr_GSE70691))
#ee2 <- na.omit(ee)  #na.omit处理对象中缺失的值

#对临床信息的合并
GSE13507_cc2 <- GSE13507_clin[,c('ID','OS','OS.time')]%>% mutate(Cohort='GSE13507')  #mutate函数创建、修改和删除列
GSE70691_cc2 <- GSE70691_clin[,c('ID','OS','OS.time')]%>% mutate(Cohort='GSE70691')
cc <- Reduce(rbind,list(GSE13507_cc2,GSE70691_cc2))
ee <- ee[,cc$ID]

#批次校正前绘制PCA预测队列差异
library(sva)
library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(ee), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",
             pointshape=21,
             fill.ind = cc$Cohort,
             palette = "npg", 
             col.ind = cc$Cohort,
             alpha.ind=0.7,
             addEllipses = TRUE, 
             legend.title = "Cohort",
             title = 'Before Removing Batch'
)+theme_bw(base_rect_size = 1.5)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12,colour = 'black'),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'))
# ggsave(filename = 'Results/6. bulk_NTP&survival/Before-batch-PCA.pdf',width = 5.8,height = 4.4)

### 校正批次效应,model可以有也可以没有，
### 如果有，也就是告诉combat，有些分组本来就有差别，不要给我矫枉过正
#model <- model.matrix(~as.factor(cc$OS))
batdata <- ComBat(dat = ee, batch = cc$Cohort)

#批次校正后绘制PCA预测队列差异
ddb.pca <- PCA(t(batdata), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",
             pointshape=21,
             fill.ind = cc$Cohort,
             palette = "npg", 
             col.ind = cc$Cohort,
             alpha.ind=0.7,
             addEllipses = TRUE, 
             legend.title = "Cohort",
             title = 'After Removing Batch'
)+theme_bw(base_rect_size = 1.5)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12,colour = 'black'),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'black'))
# ggsave(filename = 'Results/6. bulk_NTP&survival/After-batch-PCA.pdf',width = 5.8,height = 4.4)

dd <- batdata%>%as.data.frame()

#ID转换
load('Results/6. bulk_NTP&survival/GRCh38_v39.gtf.rda')
aa <- gtf_data[,c("gene_id","gene_name")]
aa <- aa[!duplicated(aa$gene_id),]  #基因名去重
aa$gene_id <- substr(aa$gene_id,1,15)  #substr函数截取字符串
dd <- merge(aa,dd,by.x = 1,by.y = 0)
dd <-  dd%>%
  dplyr::select(-gene_id)
dd <- dd[!duplicated(dd$gene_name),]  #基因名去重
rownames(dd) <- dd$gene_name
dd <-  dd%>%
  dplyr::select(-gene_name)
range(dd)

save(dd,cc,file = 'Results/6. bulk_NTP&survival/Metadata.Rda')



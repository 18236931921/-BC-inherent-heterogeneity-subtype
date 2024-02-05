rm(list = ls())
library(tidyr)
library(dplyr)
library(rjags)
library(infercnv)
library(rtracklayer)
library(data.table)
###上皮细胞提取---------------------------------------------------
load('Results/1.1 Single-cell data processing/post-annotated_BC.rda')

load('Results/1.2 Single-cell data processing(normal sample)/BCN_Epithelial.rda')
#正常上皮细胞的原始数据矩阵
BCN <- fread(file = "Datasets/GSM5329919_BCN_gene_cell_exprs_table.txt")
BCN <- BCN[,-1]
BCN=aggregate(.~Symbol,max,data=BCN)
BCN <- BCN%>%
  tibble::column_to_rownames('Symbol')
BC_Normal <- BCN_Epithelial@meta.data
Con <- intersect(rownames(BC_Normal),colnames(BCN))
BC_N <- BCN[,Con]
save(BC_N,file = 'Results/3. inferCNV/BC_N.rda')

load('Results/3. inferCNV/BC_N.rda')
#肿瘤上皮的原始数据
BC_Epithelial <- subset(BC,idents ="Epithelial")
BC_Tumor <- BC_Epithelial@meta.data

BC1 <- fread(file = "Datasets/GSM4006644_BC1_gene_cell_exprs_table.txt")
BC1 <- BC1[,-1]
BC1 <- aggregate(.~Symbol,max,data=BC1)
BC1 <- BC1%>%
  tibble::column_to_rownames('Symbol')
old_cellname <- colnames(BC1)
new_cellname <- gsub('\\.','.1_',old_cellname)
colnames(BC1) <- new_cellname
BC_T1 <- subset(BC_Tumor,patient=="BC1")
BC_T1 <- BC1[,rownames(BC_T1)]
save(BC_T1,file = 'Results/3. inferCNV/BC_T1.rda')

BC2 <- fread(file = "Datasets/GSM4006645_BC2_gene_cell_exprs_table.txt")
BC2 <- BC2[,-1]
BC2=aggregate(.~Symbol,max,data=BC2)
BC2 <- BC2%>%
  tibble::column_to_rownames('Symbol')
BC_T2 <- subset(BC_Tumor,patient=="BC2")
Con <- intersect(substr(rownames(BC_T2),1,18),colnames(BC2))
BC_T2 <- BC2[,Con]
save(BC_T2,file = 'Results/3. inferCNV/BC_T2.rda')

BC3 <- fread(file = "Datasets/GSM4006646_BC3_gene_cell_exprs_table.txt")
BC3 <- BC3[,-1]
BC3=aggregate(.~Symbol,max,data=BC3)
BC3 <- BC3%>%
  tibble::column_to_rownames('Symbol')
BC_T3 <- subset(BC_Tumor,patient=="BC3")
Con <- intersect(substr(rownames(BC_T3),1,18),colnames(BC3))
BC_T3 <- BC3[,Con]
save(BC_T3,file = 'Results/3. inferCNV/BC_T3.rda')

BC4 <- fread(file = 'Datasets/GSM4006647_BC4_gene_cell_exprs_table.txt')
BC4 <- BC4[,-1]
BC4=aggregate(.~Symbol,max,data=BC4)
BC4 <- BC4%>%
  tibble::column_to_rownames('Symbol')
BC_T4 <- subset(BC_Tumor,patient=="BC4")
Con <- intersect(substr(rownames(BC_T4),1,18),colnames(BC4))
BC_T4 <- BC4[,Con]
save(BC_T4,file = 'Results/3. inferCNV/BC_T4.rda')

BC5 <- fread(file = 'Datasets/GSM4006648_BC5_gene_cell_exprs_table.txt')
BC5 <- BC5[,-1]
BC5=aggregate(.~Symbol,max,data=BC5)
BC5 <- BC5%>%
  tibble::column_to_rownames('Symbol')
BC_T5 <- subset(BC_Tumor,patient=="BC5")
Con <- intersect(substr(rownames(BC_T5),1,18),colnames(BC5))
BC_T5 <- BC5[,Con]
save(BC_T5,file = 'Results/3. inferCNV/BC_T5.rda')

BC6 <- fread(file = 'Datasets/GSM4751267_BC6_gene_cell_exprs_table.txt')
BC6 <- BC6[,-1]
BC6=aggregate(.~Symbol,max,data=BC6)
BC6 <- BC6%>%
  tibble::column_to_rownames('Symbol')
BC_T6 <- subset(BC_Tumor,patient=="BC6")
Con <- intersect(substr(rownames(BC_T6),1,18),colnames(BC6))
BC_T6 <- BC6[,Con]
save(BC_T6,file = 'Results/3. inferCNV/BC_T6.rda')

BC7 <- fread(file = 'Datasets/GSM4751268_BC7_gene_cell_exprs_table.txt')
BC7 <- BC7[,-1]
BC7=aggregate(.~Symbol,max,data=BC7)
BC7 <- BC7%>%
  tibble::column_to_rownames('Symbol')
BC_T7 <- subset(BC_Tumor,patient=="BC7")
Con <- intersect(substr(rownames(BC_T7),1,18),colnames(BC7))
BC_T7 <- BC7[,Con]
save(BC_T7,file = 'Results/3. inferCNV/BC_T7.rda')


#执行inferCNV------------------------------------------------------------
####tumor1--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
library(png)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T1.rda')

CNV1 <- intersect(rownames(BC_T1),rownames(BC_N))
BC_N <- BC_N[CNV1,]
BC_T1 <- BC_T1[CNV1,]
CNV1 <- cbind(BC_T1,BC_N)
CNV1 <- CNV1 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann1 <- data.frame(V1=c(rep("tumor",2913),
                       rep("normal",3243)))
rownames(ann1) <- colnames(CNV1)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering1 <- gtf[,c(12,1,2,3)]%>%
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')
inter <- intersect(rownames(geneOrdering1),rownames(CNV1))
geneOrdering1 <- geneOrdering1[inter,]
geneOrdering1 <- arrange(geneOrdering1,seqnames) 

CNV1 <- CNV1[inter,]
CNV1 <- CNV1[match(rownames(geneOrdering1),rownames(CNV1)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV1, #可以直接提供矩阵对象
                                    annotations_file=ann1,
                                    gene_order_file=geneOrdering1,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV1",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV
                             

####tumor2--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T2.rda')

CNV2 <- intersect(rownames(BC_T2),rownames(BC_N))
BC_N <- BC_N[CNV2,]
BC_T2 <- BC_T2[CNV2,]
CNV2 <- cbind(BC_T2,BC_N)
CNV2 <- CNV2 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann2 <- data.frame(V1=c(rep("tumor",4783),
                        rep("normal",3243)))
rownames(ann2) <- colnames(CNV2)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering2 <- gtf[,c(12,1,2,3)]
geneOrdering2 <- geneOrdering2 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')
inter <- intersect(rownames(geneOrdering2),rownames(CNV2))
geneOrdering2 <- geneOrdering2[inter,]
geneOrdering2 <- arrange(geneOrdering2,seqnames) 

CNV2 <- CNV2[inter,]
CNV2 <- CNV2[match(rownames(geneOrdering2),rownames(CNV2)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV2, #可以直接提供矩阵对象
                                    annotations_file=ann2,
                                    gene_order_file=geneOrdering2,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV2dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV


####tumor3--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T3.rda')

CNV3 <- intersect(rownames(BC_T3),rownames(BC_N))
BC_N <- BC_N[CNV3,]
BC_T3 <- BC_T3[CNV3,]
CNV3 <- cbind(BC_T3,BC_N)
CNV3 <- CNV3 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann3 <- data.frame(V1=c(rep("tumor",2907),
                        rep("normal",3243)))
rownames(ann3) <- colnames(CNV3)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering3 <- gtf[,c(12,1,2,3)]

geneOrdering3 <- geneOrdering3 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')

inter <- intersect(rownames(geneOrdering3),rownames(CNV3))
geneOrdering3 <- geneOrdering3[inter,]
geneOrdering3 <- arrange(geneOrdering3,seqnames) 
CNV3 <- CNV3[inter,]
CNV3 <- CNV3[match(rownames(geneOrdering3),rownames(CNV3)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV3, #可以直接提供矩阵对象
                                    annotations_file=ann3,
                                    gene_order_file=geneOrdering3,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV3dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV


####tumor4--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T4.rda')

CNV4 <- intersect(rownames(BC_T4),rownames(BC_N))
BC_N <- BC_N[CNV4,]
BC_T4 <- BC_T4[CNV4,]
CNV4 <- cbind(BC_T4,BC_N)
CNV4 <- CNV4 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann4 <- data.frame(V1=c(rep("tumor",4871),
                        rep("normal",3243)))
rownames(ann4) <- colnames(CNV4)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering4 <- gtf[,c(12,1,2,3)]

geneOrdering4 <- geneOrdering4 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')

inter <- intersect(rownames(geneOrdering4),rownames(CNV4))
geneOrdering4 <- geneOrdering4[inter,]
geneOrdering4 <- arrange(geneOrdering4,seqnames) 
CNV4 <- CNV4[inter,]
CNV4 <- CNV4[match(rownames(geneOrdering4),rownames(CNV4)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV4, #可以直接提供矩阵对象
                                    annotations_file=ann4,
                                    gene_order_file=geneOrdering4,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV4dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV

####tumor5--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T5.rda')

CNV5 <- intersect(rownames(BC_T5),rownames(BC_N))
BC_N <- BC_N[CNV5,]
BC_T5 <- BC_T5[CNV5,]
CNV5 <- cbind(BC_T5,BC_N)
CNV5 <- CNV5 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann5 <- data.frame(V1=c(rep("tumor",8116),
                        rep("normal",3243)))
rownames(ann5) <- colnames(CNV5)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering5 <- gtf[,c(12,1,2,3)]

geneOrdering5 <- geneOrdering5 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')

inter <- intersect(rownames(geneOrdering5),rownames(CNV5))
geneOrdering5 <- geneOrdering5[inter,]
geneOrdering5 <- arrange(geneOrdering5,seqnames) 
CNV5 <- CNV5[inter,]
CNV5 <- CNV5[match(rownames(geneOrdering5),rownames(CNV5)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV5, #可以直接提供矩阵对象
                                    annotations_file=ann5,
                                    gene_order_file=geneOrdering5,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV5dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV


####tumor6--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T6.rda')

CNV6 <- intersect(rownames(BC_T6),rownames(BC_N))
BC_N <- BC_N[CNV6,]
BC_T6 <- BC_T6[CNV6,]
CNV6 <- cbind(BC_T6,BC_N)
CNV6 <- CNV6 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann6 <- data.frame(V1=c(rep("tumor",5969),
                        rep("normal",3243)))
rownames(ann6) <- colnames(CNV6)


#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering6 <- gtf[,c(12,1,2,3)]

geneOrdering6 <- geneOrdering6 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')

inter <- intersect(rownames(geneOrdering6),rownames(CNV6))
geneOrdering6 <- geneOrdering6[inter,]
geneOrdering6 <- arrange(geneOrdering6,seqnames) 
CNV6 <- CNV6[inter,]
CNV6 <- CNV6[match(rownames(geneOrdering6),rownames(CNV6)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV6, #可以直接提供矩阵对象
                                    annotations_file=ann6,
                                    gene_order_file=geneOrdering6,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV6dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) #是否基于HMM预测CNV

####tumor7--------------------------------------------
rm(list = ls())
library(dplyr)
library(infercnv)
library(rtracklayer)
load('Results/3. inferCNV/BC_N.rda')
load('Results/3. inferCNV/BC_T7.rda')

CNV7 <- intersect(rownames(BC_T7),rownames(BC_N))
BC_N <- BC_N[CNV7,]
BC_T7 <- BC_T7[CNV7,]
CNV7 <- cbind(BC_T7,BC_N)
CNV7 <- CNV7 %>% 
  tibble::rownames_to_column('gene_name')%>%
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)%>%
  tibble::column_to_rownames("gene_name")

#样本注释信息文件
ann7 <- data.frame(V1=c(rep("tumor",5286),
                        rep("normal",3243)))
rownames(ann7) <- colnames(CNV7)

#基因位置信息文件
gtf <- import('Results/3. inferCNV/gencode.v38.annotation.gtf.gz')
gtf <- as.data.frame(gtf)
gtf <- gtf %>% 
  filter(type=="gene")
geneOrdering7 <- gtf[,c(12,1,2,3)]

geneOrdering7 <- geneOrdering7 %>% 
  mutate(rowMean = (.$end)-(.$start)) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>%
  tibble::column_to_rownames('gene_name')

inter <- intersect(rownames(geneOrdering7),rownames(CNV7))
geneOrdering7 <- geneOrdering7[inter,]
geneOrdering7 <- arrange(geneOrdering7,seqnames) 
CNV7 <- CNV7[inter,]
CNV7 <- CNV7[match(rownames(geneOrdering7),rownames(CNV7)),]

#执行CNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=CNV7, #可以直接提供矩阵对象
                                    annotations_file=ann7,
                                    gene_order_file=geneOrdering7,
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="Results/3. inferCNV/output_CNV7dir",  # 输出文件夹
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) #是否基于HMM预测CNV

















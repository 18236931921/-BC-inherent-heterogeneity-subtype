rm(list = ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(data.table)
library(tibble)
library(dplyr)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(maftools)
library(circlize)
library(ggsci)

# data processing ---------------------------
maf <-  read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T) # 读取突变数据
maf1 <- as.data.frame(maf@data)
maf1 <- maf1[,c(1,2,4:13,16)]

### 计算96种三核苷酸变化 ###
rmSilence = F #是否移除沉默突变（根据实际情况还可能移除其他SNP类型）;移除则改为T
if (rmSilence) {
  maf1 <- as.data.frame(maf1[which(maf$Variant_Type == "SNP" & maf1$Variant_Classification != "Silent"),]) #仅考虑SNP并移除沉默突变
} else {
  maf1 <- as.data.frame(maf1[which(maf1$Variant_Type == "SNP"),]) #仅使用SNP (突变签名研究本身只考虑SNP)
}

table(maf1$Chromosome)
maf1$Tumor_Sample_Barcode <- substr(maf1$Tumor_Sample_Barcode,1,12)

# maf$Chromosome <- paste0("chr",maf$Chromosome) 
# 这句代码是为保证染色体格式为chrx，符合deconstructsigs输入需求,这里已经符合要要求，无需运行

snp.count <- mut.to.sigs.input(mut.ref = maf1, 
                               sample.id = "Tumor_Sample_Barcode", 
                               chr = "Chromosome", 
                               pos = "Start_Position", 
                               ref = "Reference_Allele", # 参考位点碱基
                               alt = "Tumor_Seq_Allele2", # 突变位点碱基
                               bsg = BSgenome.Hsapiens.UCSC.hg38) # hg19参考基因组

### 计算单样本突变签名的权重 ###
cut.off <- 0.06 # 默认cutoff为6%，如例文所示
mut.wt <- data.frame()
sigs.out.list <- list()
index <- 1
for (sample in rownames(snp.count)) {
  cat(paste0(sample," starts and ",139-index," samples remain to be analyzed!\n"))
  
  tmp <- whichSignatures(tumor.ref = snp.count, 
                         signatures.ref = signatures.nature2013, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome', # 标准化方式如例文所示
                         signature.cutoff = cut.off)
  
  index <- index + 1
  
  # 生成突变签名矩阵
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt,tmp)
}

write.table(mut.wt,"Results/10.2 Mutation/mutsig.weightMatrix.txt",sep = "\t",row.names = T,col.names = NA)

# Mutant landscape -------------------------------------------------------------------
rm(list = ls())
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)

## 分组信息
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
cluster <- dd[,c(1,2)]
cluster <- cluster[order(cluster$prediction),]
colnames(cluster) <- c('ID','Cluster')
cluster$ID <- substr(cluster$ID,1,12)

maf <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T)%>%subsetMaf(tsb = cluster$ID)
ID1 <- unique(maf@data$Tumor_Sample_Barcode)
ID2 <- unique(maf@clinical.data$Tumor_Sample_Barcode)
ID3 <- intersect(ID1,ID2)
ID4 <- intersect(ID1,cluster$ID)
ID <- ID4

# maf1 <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T)%>%subsetMaf(tsb = ID)
#oncoplot(maf1,writeMatrix = T,top = 30,) #用来生成onco_matrix.txt
dd <- read.table('Results/10.2 Mutation/onco_matrix.txt',h=T,sep = '\t',check.names = F)
mut <- ifelse(dd=="",0,1)
mut <- as.data.frame(mut)  ##准备完毕

## 突变频数
ID2 <- colnames(mut)
maf2 <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T)%>%subsetMaf(tsb = ID2)
tmb <- maf2@variants.per.sample
colnames(tmb) <- c("ID","log10TMB")
tmb$log10TMB <- log10(tmb$log10TMB)

tmb <- as.data.frame(tmb)
rownames(tmb) <- tmb[,1]  ##准备完毕

## 突变签名
mutsig <- read.table("Results/10.2 Mutation/mutsig.weightMatrix.txt", row.names = 1, header = T) 

## 拷贝数GISTIC2.0结果（利用SNP segment文件在GenePattern上获得）
cna.region <- read.table("Results/10.2 Mutation/gdac.broadinstitute.org_BLCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_lesions.conf_99.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cna.gene <- read.table("Results/10.2 Mutation/gdac.broadinstitute.org_BLCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

## 亚型数据
ID <- intersect(ID,colnames(mut))
cna.gene <- as.data.frame(t(cna.gene))
cna.gene <- cna.gene%>%rownames_to_column('ID')
cna.gene$ID <- substr(cna.gene$ID,1,12)
ID1 <- intersect(ID,cna.gene$ID)  

cna.region <- as.data.frame(t(cna.region))
cna.region <- cna.region%>%rownames_to_column('ID')
cna.region[9:416,]$ID <- substr(cna.region[9:416,]$ID,1,12)
ID2 <- intersect(ID1,cna.region$ID) 

subt <- cluster[cluster$ID%in%ID2,]
table(subt$Cluster)

##分组信息不再处理
a <- subt[,c(1,2)]
rownames(a) <- NULL
subt <- a%>%column_to_rownames('ID')
table(subt$Cluster)  
# C1 C2 C3 C4 
# 117    72    35    49 

cna.gene2 <- cna.gene[cna.gene$ID%in%ID2,]
cna.gene3 <- rbind.data.frame(cna.gene[1:2,],cna.gene2)
rownames(cna.gene3) <- cna.gene3[,1]
cna.gene <- cna.gene3[,-1]
cna.gene <- as.data.frame(t(cna.gene)) 

cna.region2 <- cna.region[cna.region$ID%in%ID2,]
cna.region3 <- rbind.data.frame(cna.region[1:8,],cna.region2)
rownames(cna.region3) <- cna.region3[,1]
cna.region <- cna.region3[,-1]
cna.region <- as.data.frame(t(cna.region)) 

mutsig <- mutsig%>%rownames_to_column('ID')
mutsig <- mutsig[mutsig$ID%in%ID2,]
rownames(mutsig) <- mutsig[,1]
mutsig <- mutsig[,-1] 

# 设置亚型颜色
clust.col <- c('#436b95','#b81d25','#2f7e2f','#800080')
blue   <- '#ee9269'
  red    <- '#7dbea8'
    
  # 处理突变签名数据
  mutsig <- mutsig[,c("Signature.1B","Signature.2","Signature.5","Signature.13")] # 参考文献PMID：23945592，找到该肿瘤对应的signature
  
  mutsig$Cluster <- subt[rownames(mutsig),"Cluster"] # 添加亚型结果
  mutsig <- mutsig[order(mutsig$Cluster,decreasing = F),] # 确定整个热图的排序，High-Low Group
  
  # 挑选要展示的基因
  mutgene <- rownames(dd)[1:15]
  # 制作oncoprint的输入数据
  onco.input <- mut[mutgene,rownames(mutsig),]
  onco.input[onco.input == 1] <- "Mutated" # 二值矩阵中1记为突变
  onco.input[onco.input != "Mutated"] <- "" # 非“突变”给予空值
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
    },
    Mutated = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A60000", col = "#A60000")) 
    }
  )
  col = c("Mutated" ="#A60000") # 突变颜色，注意这里只给了主图像的图例
  
  my_ann <- subt[rownames(mutsig),,drop = F]
  my_annotation = HeatmapAnnotation(df = my_ann, 
                                    col = list(Cluster = c("C1" = clust.col[1],
                                                           "C2" = clust.col[2],
                                                           "C3" = clust.col[3],
                                                           "C4" = clust.col[4])))
  
  # 突变主区域的上部注释（突变负荷柱状图）
  top_anno <- anno_barplot(as.numeric(tmb[rownames(mutsig),"log10TMB"]),
                           border = FALSE,
                           gp = gpar(fill ='#7dbea8',#"#3379B4",
                                     border =NA,lty="blank"), 
                           height = unit(2.5, "cm"))
  
  # 突变主区域的上部注释（突变签名柱状图）
  tmp <- mutsig[,c("Signature.1B","Signature.2","Signature.5","Signature.13")] # 只取感兴趣Signature
  tmp$Others <- 1 - rowSums(tmp) # 计算其他签名的比例
  top_anno2 <- anno_barplot(as.matrix(tmp),
                            border = FALSE,
                            gp = gpar(fill = c('#b7795e','#5eb7a5','#b75e9d','#800080'), 
                                      border = NA, # 无边框
                                      lty = "blank"),
                            height = unit(2, "cm")) # 高度
  
  tmp <- as.data.frame(t(mut[mutgene,rownames(mutsig),]))
  mut.order <- names(sort(colSums(tmp),decreasing = T)) # 根据突变频数高低排序展示突变的顺序
  tmp$Cluster <- subt[rownames(tmp),"Cluster"]
  pct <- NULL # 计算各个基因突变的百分比
  for (i in mut.order) {
    tmp1 <- tmp[,c(i,"Cluster")]
    tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$Cluster))[2,]/sum(tmp1[,1])
    pct <- rbind.data.frame(pct,tmp1)
  }
  rownames(pct) <- mut.order
  
  # 添加右侧百分比堆积柱状图
  right_anno <- anno_barplot(as.matrix(pct),
                             which = "row",
                             border = FALSE,
                             gp = gpar(fill = clust.col,border=NA,lty="blank"), 
                             bar_width = 0.6,
                             width = unit(1.8, "cm"),
                             height = unit(1, "cm"))
  
  #rownames(my_ann) <- my_ann[,1]
  dd <- onco.input[mut.order,rownames(my_ann)]
  op1 <- oncoPrint(dd, # 排序的突变矩阵
                   alter_fun = alter_fun,  # 主区域的函数，包括各单元格大小、背景颜色等等
                   col = col, # 突变颜色
                   bottom_annotation = NULL, # 无底部注释
                   top_annotation = c(HeatmapAnnotation(TMB = top_anno), # 顶部第一个注释：TMB
                                      my_annotation, # 顶部第二个注释：亚型
                                      HeatmapAnnotation(MutSig = top_anno2)), # 顶部第三个注释：突变签名
                   column_order = rownames(my_ann), # 样本的排序，根据突变签名的顺序
                   right_annotation = rowAnnotation(PCT = right_anno), # 右侧堆叠柱状图注释
                   show_pct = T, # 展示左侧的百分比
                   column_title = "", # 不显示主题
                   show_heatmap_legend = T, # 展示图例
                   column_split = my_ann$Cluster, # 根据亚型切分热图
                   column_title_gp = gpar(fontsize = 8),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8))
  op1
  
  
  # 拷贝数变异 ---
  # 选择要展示的拷贝数区域
  cna <- cna.region[1:(nrow(cna.region)/2),c(1,8,9:(ncol(cna.region)))] # 选取有效列
  cna <- distinct(cna,Descriptor,.keep_all = T)
  rownames(cna) <- paste0(gsub(" ","",cna$Descriptor),"-", substr(rownames(cna),1,3)) # 重命名行以确定扩增和缺失的位点
  cna.modified <- cna[1:nrow(cna),3:ncol(cna)]
  cna.count <- sapply(cna.modified, as.numeric) %>% 
    as.data.frame()
  rownames(cna.count) <- rownames(cna.modified)
  cna.count[cna.count == 2] <- 1
  cna.count$sum <- rowSums(cna.count)
  cna.order <- cna.count[order(cna.count$sum,decreasing = T),]
  amp <- rownames(cna.order)[1:15][which(rownames(cna.order)%like%'Amp')] %>% na.omit()
  del <- rownames(cna.order)[1:15][which(rownames(cna.order)%like%'Del')] %>% na.omit()
  lesion.sig <- c(amp,del)
  
  onco.input2 <- cna.modified[lesion.sig,rownames(mutsig)] # 选取要展示的拷贝数变异
  tmp11 <- onco.input2[1:6,] # 前6个为扩增
  
  tmp1 <- as.data.frame(lapply(tmp11,as.numeric))
  rownames(tmp1) <- rownames(tmp11)
  colnames(tmp1) <- colnames(tmp11) 
  
  tmp1[tmp1 == 1] <- "Gain" # 数值大于0即为Gain
  tmp1[tmp1 == 2] <- "Gain"
  tmp1[tmp1 == 0] <- ""
  
  tmp22 <- onco.input2[7:15,] # 后9个为缺失
  tmp2 <- as.data.frame(lapply(tmp22,as.numeric))
  rownames(tmp2) <- rownames(tmp22)
  colnames(tmp2) <- colnames(tmp22) 
  
  tmp2[tmp2 == 1] <- "Loss"
  tmp2[tmp2 == 2] <- "Loss"
  tmp2[tmp2 == 0] <- ""
  onco.input2 <- rbind.data.frame(tmp1,tmp2)
  
  alter_fun2 = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
    },
    Gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = red, col = red)) 
    },
    Loss = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = blue, col = blue)) 
    }
  )
  col2 = c("Gain" = red,
           "Loss" = blue)
  # 确定展示的顺序（看自己喜好，我这里是按照臂的顺序来的）
  
  lesion.order <- lesion.sig
  tmp <- as.data.frame(t(cna.modified[lesion.order,rownames(mutsig),]))
  
  tmp3 <- as.data.frame(lapply(tmp,as.numeric))
  rownames(tmp3) <- rownames(tmp)
  colnames(tmp3) <- colnames(tmp) 
  
  tmp3[tmp3 > 0] <- 1 # 所有大于1的均改为1以便计算变异频数
  tmp <- tmp3
  
  tmp$Cluster <- as.character(subt[rownames(tmp),"Cluster"])
  pct <- NULL
  for (i in lesion.order) {
    tmp111 <- tmp[,c(i,"Cluster")]
    tmp111 <- as.data.frame.array(table(tmp111[,1],tmp111$Cluster))[2,]/sum(tmp111[,1])
    pct <- rbind.data.frame(pct,tmp111)
  }
  rownames(pct) <- lesion.order
  
  # 右侧堆叠百分比柱状图
  right_anno2 <- anno_barplot(as.matrix(pct),
                              which = "row",
                              border = FALSE,
                              gp = gpar(fill = clust.col,
                                        border = NA,
                                        lty = "blank"), 
                              bar_width = 0.6,
                              width = unit(1.8, "cm"),
                              height = unit(1, "cm"))
  
  # 同样的方式绘制热图
  op2 <- oncoPrint(onco.input2[lesion.order,rownames(my_ann)], 
                   alter_fun = alter_fun2, 
                   col = col2, 
                   bottom_annotation = NULL, 
                   top_annotation = NULL,
                   column_order = rownames(my_ann),
                   right_annotation = rowAnnotation(PCT = right_anno2),
                   row_order = lesion.order, 
                   show_pct = T,
                   column_title = "", 
                   show_heatmap_legend = F, 
                   column_split = my_ann$Cluster,
                   column_title_gp = gpar(fontsize = 8),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8))
  
  op2
  
  
  # 构建额外图例
  lgd.mutsig = Legend(labels = c("Sig.1B","Sig.4","Sig.6"), 
                      title = "MutSig", 
                      legend_gp = gpar(fill = c('#b7795e','#5eb7a5','#b75e9d','#800080')))
  
  lgd.cna.region = Legend(labels = c("Gain","Loss"), 
                          title = "CNV",  #(arm-level)
                          legend_gp = gpar(fill = c(red,blue)))
  
  
  
  lgd_list <- list(lgd.mutsig, lgd.cna.region)#,lgd.cna.gene
  
  # 合并热图
  pdf("Results/10.2 Mutation/Mutational landscape in TCGA_BLCA.pdf", width = 10,height = 12)
  draw(op1 %v% op2, # 垂直叠加热图%v% op3
       annotation_legend_list = lgd_list) # 添加自定义的图例
  invisible(dev.off())
  
  
  
# gene mutation frequency -----------------------------------------------
  rm(list = ls())
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(ggsci)
  library(circlize)
  library(maftools)
  ## 分组信息
  load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
  cluster <- dd[,c(1,2)]
  cluster <- cluster[order(cluster$prediction),,drop=F]
  colnames(cluster) <- c('ID','Cluster')
  cluster$ID <- substr(cluster$ID,1,12)
  
  maf <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T)%>%subsetMaf(tsb = cluster$ID)
  ID1 <- unique(maf@data$Tumor_Sample_Barcode)
  ID2 <- unique(maf@clinical.data$Tumor_Sample_Barcode)
  ID3 <- intersect(ID1,ID2)
  ID4 <- intersect(ID1,cluster$ID)
  ID <- ID4
  
  # maf1 <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz',isTCGA = T)%>%subsetMaf(tsb = ID)
  # oncoplot(maf = maf1,top=15,writeMatrix = T,removeNonMutated = F)
  pp <- t(read.table('Results/10.2 Mutation/onco_matrix2.txt',h=T,row.names = 1,sep = '\t',check.names = F))%>%as.data.frame()
  pp[pp=='0'] <- 'x'
  pp[pp==''] <- 'x'
  ID <- cluster[cluster$ID %in% ID,]
  C1 <- C2 <- C3 <- C4 <-c()
  for (i in colnames(pp)) {
    C1[i] <- sum(pp[ID$ID[ID$Cluster=='C1'],i]!='x')/length(ID$ID[ID$Cluster=='C1'])
    C2[i] <- sum(pp[ID$ID[ID$Cluster=='C2'],i]!='x')/length(ID$ID[ID$Cluster=='C2'])
    C3[i] <- sum(pp[ID$ID[ID$Cluster=='C3'],i]!='x')/length(ID$ID[ID$Cluster=='C3'])
    C4[i] <- sum(pp[ID$ID[ID$Cluster=='C4'],i]!='x')/length(ID$ID[ID$Cluster=='C4'])
    
  }
  my <- data.frame(C1=C1,C2=C2,C3=C3,C4=C4)
  my$p <- NA
  my$sig <- NA
  for (i in 1:nrow(my)) {
    p <- chisq.test(data.frame(x=c(my[i,1]*144,144*(1-my[i,1])),y=c(my[i,2]*59,59*(1-my[i,2])),z=c(my[i,3]*152,152*(1-my[i,3]))))$p.value
    my$p[i] <- p
    my$sig[i] <- ifelse(p<0.0001,'****',ifelse(p<0.001,'***',ifelse(p<0.01,'**',ifelse(p<0.05,'*',''))))
  }
  my$ID <- rownames(my)
  #my2 <- pivot_longer(my,1:2,names_to = 'Group',values_to = 'Rate')
  
  Cluster <- c('#436b95','#b81d25', '#2f7e2f','#800080')
  names(Cluster) <- c('C1','C2','C3','C4')
  
  my$C1 <- sprintf('%0.2f',my$C1)%>%as.numeric()
  my$C2 <- sprintf('%0.2f',my$C2)%>%as.numeric()
  my$C3 <- sprintf('%0.2f',my$C3)%>%as.numeric()
  my$C4 <- sprintf('%0.2f',my$C4)%>%as.numeric()
  hh <- as.matrix(my[,1:4])
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", hh[i, j]), x, y, gp = gpar(fontsize = 15))
  }
  right <- HeatmapAnnotation(dd = anno_text(x = my$sig,which = 'row',just = -0.1),which = 'row')
  top <- HeatmapAnnotation(Cluster = anno_block(gp = gpar(fill = c('#436b95','#b81d25', '#2f7e2f','#800080'),lwd = 2),
                                                labels = c('C1','C2','C3','C4'),
                                                labels_gp = gpar(col = "white", fontsize = 11,just = "center")),
                           border = T,
                           annotation_name_side = 'right',
                           annotation_name_gp = gpar(fontsize=8),
                           show_annotation_name = F,
                           na_col = NA,
                           show_legend = F,
                           gap = unit(1, "mm"),#'#436b95','#b81d25','#2f7e2f'
                           col = list(Cluster = c("C1" = '#436b95', "C2" = '#b81d25', "C3" = '#2f7e2f', "C4" = '#800080')))
  plo <- Heatmap(hh,cell_fun = cell_fun,
                 border = T,
                 top_annotation = top,
                 right_annotation = right,
                 col = colorRamp2(c(0,0.2,0.5),colors = c('white',pal_npg("nrc", alpha =0.6)(9)[5],pal_npg("nrc", alpha =0.8)(9)[1])),
                 column_split = 1:4,
                 rect_gp = gpar(color='balck'),
                 column_title = NULL,show_column_names = F,
                 column_title_rot = 90,
                 cluster_rows = F,cluster_columns = F,
                 row_names_side = 'left',
                 row_names_gp = gpar(size = 10,color='balck'),
                 heatmap_legend_param = list(title= "Frequency",
                                             at=c(0,0.3,0.6),
                                             title_position = "topleft")
  )
  plo
  draw(plo,heatmap_legend_side = "right")
  
  library(export)
  graph2pdf(file='Results/10.2 Mutation/Mutational-Frequency.pdf',height=6,width=8)
  
  
# FGG ----------------------------
  rm(list = ls())
  library(magrittr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(Rmisc)
  library(dplyr)
  library(tidyverse)
  library(ggsci)
  Segment <- read.delim("Results/10.2 Mutation/TCGA-BLCA.cnv.tsv.gz")
  Segment$bases=Segment$End-Segment$Start
  
  load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
  cluster <- dd[,c(1,2)]
  colnames(cluster) <- c("ID","group")
  
  data=data.frame()
  for (i in 1:length(table(Segment$sample))) {
    tmp=Segment[Segment$sample==names(table(Segment$sample))[i],]
    FGA=sum(tmp[abs(tmp$value)>0.2,"bases"])/ sum(tmp[,"bases"])
    FGG=sum(tmp[tmp$value>0.2,"bases"])/sum(tmp[,"bases"])
    FGL=sum(tmp[tmp$value< -0.2,"bases"])/sum(tmp[,"bases"])
    tmp=data.frame(Patient=names(table(Segment$sample))[i],FGA=FGA,FGG=FGG,FGL=FGL)
    data=rbind(data,tmp)
  }
  data$Patient <-substr(data$Patient,1,12) 
  data <- data%>%
    distinct(Patient, .keep_all = TRUE)
  
  tcga <- cluster
  sub <- tcga
  rownames(sub) <- NULL
  sub$ID <- substr(sub$ID,1,12)
  data2 <- merge(sub,data,by.x=1,by=1)
  data2 <- data2[,c("ID","group","FGA","FGG","FGL")]
  save(data2,file = "Results/10.2 Mutation/fraction of genomic alteration.rda")
  
  
  load("Results/10.2 Mutation/fraction of genomic alteration.rda")
  ggdata <- pivot_longer(data2,3:5,names_to = 'FF',values_to = 'VV')
  my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
  p1 <- ggplot(ggdata,aes(group,VV,color=group))+
    geom_boxplot(outlier.colour = NA,aes(fill=group),color='black',size=0.6,alpha=0.65)+
    geom_violin(alpha=0.5,aes(fill=group),color=NA,trim = T)+
    geom_signif(comparisons = my_comparisons,step_increase = 0.12,#两组直接对比框的高度
                map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
                test =t.test,size=0.9,#线条粗细
                textsize = 2.5,#P值字号大小
                position ='identity', #位置
                show.legend=FALSE,
                aes(col=group))+
    facet_wrap(~FF)+
    scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
    scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
    theme_bw(base_rect_size = 1.8)+
    labs(y='Percent Genome Altered')+
    theme(axis.text.y = element_text(size = 12,colour = 'black'),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 18,colour = 'black',face = 'bold'),
          axis.title.x = element_blank(),
          strip.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size=12,colour = 'black'),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          panel.background = element_rect(fill='#f3f6f6'))
  p1
  ggsave(filename = 'Results/10.2 Mutation/FGA_FGG_FGL.pdf',width = 7,height = 5)
  
  
# TMB+SNP+INDEL ------------------------------------------
  BLCA <- read.maf('Results/10.2 Mutation/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf.gz', isTCGA = T)
  dd <- BLCA@data
  write.table(dd,'Results/10.2 Mutation/All_maf.txt',row.names = F,sep = '\t',quote = F)
  maf <- read.maf('Results/10.2 Mutation/All_maf.txt',isTCGA = T)
  
  # 准备group
  load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
  cluster <- dd[,c(1,2)]
  cluster <- cluster[order(cluster$prediction),]
  colnames(cluster) <- c('ID','Cluster')
  cluster$ID <- substr(cluster$ID,1,12)
  
  dd <- maf@variant.type.summary%>%as.data.frame()
  dd <- merge(dd,cluster,by=1)
  ID <- dd$Tumor_Sample_Barcode%>%as.character()
  dd$INDEL <- dd$INS+dd$DEL
  dd <- dd[,c(6,5,4,7)]
  colnames(dd)[2] <- 'TMB'
  dd2 <- pivot_longer(dd,2:4,names_to = 'tt',values_to = 'count')
  dd2$tt <- factor(dd2$tt,levels = c('TMB','SNP','INDEL'))

  my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
  library(ggsci)
  ggplot(dd2,aes(Cluster,count+0.1,color=Cluster))+
    geom_boxplot(outlier.colour = NA,aes(fill=Cluster),color='black',size=0.6,alpha=0.65)+
    geom_violin(alpha=0.5,aes(fill=Cluster),color=NA,trim = T)+
    scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+ 
    scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
    scale_y_log10()+
    geom_signif(comparisons = my_comparisons,step_increase = 0.12,#两组直接对比框的高度
                map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
                test = t.test,size=0.85,#线条粗细
                textsize = 2.5,#P值字号大小
                position ='identity', #位置
                show.legend=FALSE,
                aes(col=Cluster)
    )+
    facet_wrap(~tt)+
    expand_limits(y=max(dd2[,'count'])*1.1)+
    theme_bw(base_rect_size = 2)+
    labs(y='Mutation Counts')+
    theme(axis.text.y = element_text(size = 12,colour = 'black'),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 18,colour = 'black',face = 'bold'),
          axis.title.x = element_blank(),
          strip.text = element_text(size=14),
          legend.title = element_blank(),
          legend.text = element_text(size=12,colour = 'black'),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          panel.background = element_rect(fill='#f3f6f6'))
  ggsave(filename = 'Results/10.2 Mutation/TMB.pdf',width =7,height = 5)
  
  #TMB
  load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
  maf <- read.maf('Results/10.2 Mutation/All_maf.txt',isTCGA = T)
  BLCA_tmb = tmb(maf = maf)
  BLCA_tmb$Tumor_Sample_Barcode <- substr(BLCA_tmb$Tumor_Sample_Barcode,1,12)
  exp <- data
  colnames(exp) <- str_sub(colnames(exp),1,12)
  dd$ID <- str_sub(dd$ID,1,12)
  id <- intersect(colnames(exp),dd$ID)
  
  rr2 <- dd %>% 
    .[dd$ID %in% id,1:2,drop=F] 
  exp <- exp[,id] %>% 
    .[,match(rr2$ID,colnames(.))] 
  
  corinput <- merge(rr2,BLCA_tmb,by = 1)
  colnames(corinput)[1:2] <- c('sample','Cluster')
  colnames(corinput)
  
  corinput$Cluster <- factor(corinput$Cluster)
  colnames(corinput)[5] <- 'TMB_log'
  
  my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
  corinput %>%
    ggplot(aes(Cluster,TMB_log)) +
    geom_boxplot(aes(color=Cluster), outlier.colour = NA, size=1, fill=NA) +
    geom_jitter(aes(fill=Cluster,color=Cluster), width = 0.25, shape=21, size=2, alpha=0.7) +
    geom_signif(comparisons = my_comparisons,step_increase = 0.12,#两组直接对比框的高度
                map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
                test = t.test,size=0.85,#线条粗细
                textsize = 2.5,#P值字号大小
                position ='identity', #位置
                show.legend=FALSE,
                aes(col=Cluster)
    )+
    theme_bw(base_rect_size = 1.5) +
    labs(x=NULL,y='Relative Expression', title = 'BLCA-TMB') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5,size = 12),
          axis.title.y = element_text(size = 12,colour = 'black'),
          axis.text.x = element_text(size = 12,colour = 'black'),
          axis.text.y = element_text(size = 12,colour = 'black'),
          axis.ticks = element_line(size = 1),
          panel.grid = element_blank()) +
    scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
    scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))
  
  ggsave(filename = 'Results/10.2 Mutation/TMB_log-boxplot.pdf', width = 2.5, height = 4)
  
  
# Neoantigen -----------------------------
  rm(list = ls())
  gc()
  library(data.table)
  library(tibble)
  library(ggpubr)
  load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
  expr <- data
  colnames(expr) <- substr(colnames(expr),1,12)
  
  nal <- fread('Results/10.2 Mutation/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv',data.table = F)
  nal <- data.frame(ID=nal$sample,NAL=nal$neoantigen_num)
  nal <- nal[nal$ID %in% colnames(expr),]
  
  Cluster <- dd[,c(1,2)]
  Cluster$ID <- substr(Cluster$ID,1,12)
  
  nal_df <- merge(nal,Cluster,by = 1)
  nal_df$NAL <- log(nal_df$NAL+1)
  colnames(nal_df)[3] <- 'Cluster'
  nal_df$Cluster <- factor(nal_df$Cluster )
  colnames(nal_df) <- c('ID','NAL_log','Cluster')
  
  #箱型图
  my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
  nal_df %>%
    ggplot(aes(Cluster,NAL_log)) +
    geom_boxplot(aes(color=Cluster), outlier.colour = NA, size=1, fill=NA) +
    geom_jitter(aes(fill=Cluster,color=Cluster), width = 0.25, shape=21, size=2, alpha=0.7) +
    geom_signif(comparisons = my_comparisons,step_increase = 0.12,#两组直接对比框的高度
                map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
                test = t.test,size=0.85,#线条粗细
                textsize = 2.5,#P值字号大小
                position ='identity', #位置
                show.legend=FALSE,
                aes(col=Cluster)
    )+
    theme_bw(base_rect_size = 1.5) +
    labs(x=NULL,y='Relative Expression', title = 'Neoantigen') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5,size = 12),
          axis.title.y = element_text(size = 12,colour = 'black'),
          axis.text.x = element_text(size = 12,colour = 'black'),
          axis.text.y = element_text(size = 12,colour = 'black'),
          axis.ticks = element_line(size = 1),
          panel.grid = element_blank()) +
    scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
    scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))
  
  ggsave(filename = 'Results/10.2 Mutation/Neoantigen-boxplot.pdf', width = 2.3, height = 4)
  
  
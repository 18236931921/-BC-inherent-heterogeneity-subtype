rm(list = ls())
library(GSVA)
library(tibble)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)
library(circlize)
library(ggplot2)
library(grid)
library(gridExtra)
library(magrittr)
library(tidyverse)
library(ggpubr)
library(data.table)
library(clusterProfiler)


# Calculation of Immune cell infiltration -----------------------------------------------------------------
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
load('Results/10.1 TIME/cellMarker_ssGSEA.Rdata')
TCGA_expr <- data
ssGSEA_output <- gsva(as.matrix(TCGA_expr),cellMarker, method = "ssgsea")
ssGSEA_output <- as.data.frame(t(ssGSEA_output))
save(ssGSEA_output,file = 'Results/10.1 TIME/ssGSEA_output.Rda')

# Immune cell infiltration ---------------------------------------------------------------------
load('Results/10.1 TIME/ssGSEA_output.Rda')
TCGA_clin <- dd[,c(1,2)]
ssGSEA_output <- rownames_to_column(ssGSEA_output,var = "ID")
ssGSEA_output <- merge(ssGSEA_output,TCGA_clin,by = 1)
colnames(ssGSEA_output)[30] <- c('cluster')
heatmapinput <- ssGSEA_output
heatmapinput <- arrange(heatmapinput,cluster)
df <- column_to_rownames(heatmapinput,var = "ID")

dd2 <- pivot_longer(df,cols = 1:28,names_to = 'Cells',values_to = 'Value')
colnames(dd2)[1] <- c('Subtype')
ggplot(dd2,aes(Cells,Value,fill = Subtype))+
  geom_boxplot(outlier.colour = 'black',
               outlier.size = 0.1,
               outlier.alpha = 0.6,
               size=0.4)+
  stat_compare_means(label = 'p.signif',method = "anova",vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")))+
  theme_bw(base_rect_size = 1.5)+
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 11),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust = 1,colour = 'black'),
        axis.title.y = element_text(size=16))+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  #guides(fill=F)+
  labs(fill='Subtype',y='Relative Expression',x=NULL)
ggsave(filename = 'Results/10.1 TIME/Immune_cell_box.pdf',height = 5,width = 12)

#热图
ml <- df[, c(1:28)] #画图的ssGSEA数据
ml <- as.data.frame(t(apply(ml, 2, scale))) #scale标化
colnames(ml) <- heatmapinput$ID
library(RColorBrewer)
mycolor1<-"#436b95"
mycolor2<-"#b81d25"
mycolor3<-"#2f7e2f"
mycolor4<-"#800080"

ml[ml>2] <- 2
ml[ml< -2] <-  -2
col_fun <- colorRamp2(c(-2, 0, 2), c("#376EB3",'white','#C95550'))#蓝色——白——红色
#col_fun <- colorRamp2(c(-2, 0, 2), c("#21b6af",'white',"#eeba4d"))#设置绘图的颜色，可以根据数据实际情况调整
#大于2的值都映射为eeba4d色黄色，小于-2的值都映射为21b6af色绿色    绿色 ——白——黄色
#c('#0080FF','white','#EA0000'))#蓝色——白——红色
#临床信息注释
Cluster <- heatmapinput[, "cluster"]

ha = HeatmapAnnotation(Cluster=Cluster,
                       col = list(Cluster = c("C1" = mycolor1 ,#黄色
                                              "C2" = mycolor2, #绿色
                                              "C3" = mycolor3,#红色
                                              "C4" = mycolor4)#蓝
                       ),na_col = "white",   #NA值的颜色由na_col参数控制
                       show_legend = rep(TRUE, 1),#是否要显示annotation legend 图例
                       annotation_height = unit(rep(3, 4), "mm"),#临床annotation的高度
                       annotation_legend_param = list(Cluster = list(title = "Cluster")))

ml <- as.matrix(ml)
ht <- Heatmap(ml, col = col_fun,
              name = "Z-score",
              cluster_rows = TRUE, cluster_columns = F,
              show_row_names = TRUE, show_column_names = FALSE,
              top_annotation = ha, column_title = 'TCGA BLCA',  #"TCGA BLCA",#膀胱尿路上⽪癌
              clustering_method_columns = "ward.D2",#可用于指定进行层次聚类的方法。允许的值是hclust()函数支持的值，包括 "ward.D2"，“single”, “complete”， “average”
              clustering_distance_columns = "euclidean",
              clustering_distance_rows = "euclidean",
              clustering_method_rows  = "ward.D2", column_dend_height = unit(20, "mm"),
              column_split = heatmapinput$Cluster,
              column_gap = unit(2, "mm"),
              border = T
              
              
)
pdf("Results/10.1 TIME/immun_infiltration_ssGSEA_heatmap_TCGA.pdf", 7, 5)
draw(ht, annotation_legend_side = "left", heatmap_legend_side = "left")
dev.off()


# Immune check point(ICP) -------------------------------------------------------------------
ic <- ssGSEA_output[,c('ID','cluster',"Macrophage","Neutrophil",
                       "Natural killer cell","Activated CD4 T cell",
                       "Activated CD8 T cell","Gamma delta T cell")]
IC <- rbind(data.frame(Type='B7-CD28',
                       Gene=c('CD274','CD276','CTLA4','HHLA2','ICOS',
                              'ICOSLG','PDCD1','PDCD1LG2','TMIGD2','VTCN1')),
            data.frame(Type='TNF superfamily',
                       Gene=c('BTLA','CD27','CD40','CD40LG','CD70',
                              'TNFRSF18','TNFRSF4','TNFRSF9','TNFSF14','TNSF4','TNSF9')),
            data.frame(Type='Others',
                       Gene=c('C10ORF54','ENTPD1','FGL1','HAVCR2','IDO1',
                              'LAG3','NCR3','NT5E','SIGLEC15')))
IC$Type <- factor(IC$Type,levels = unique(IC$Type))
x <- intersect(IC$Gene,rownames(TCGA_expr))
IC <- IC[IC$Gene%in%x,]
IC <- IC[order(IC$Type,IC$Gene),]
ee <- TCGA_expr[IC$Gene,]
ee <- ee[,ic$ID]

dd2 <- merge(ic,t(ee),by.x=1,by.y = 0)
# 对数据再次检查，排序
dd2[,3:35] <- scale(dd2[,3:35])
rownames(dd2) <- dd2$ID
colnames(dd2)[2] <- c('cluster')
dd2 <- dd2[order(dd2$cluster),]
dd2 <- dd2[,-1]

# 为临床信息定颜色 注释信息的构建
Group <- c('#436b95','#b81d25','#2f7e2f','#800080')
names(Group) <- c('C1','C2','C3','C4')

Top = HeatmapAnnotation(Group=anno_block(gp = gpar(fill = c('#436b95','#b81d25','#2f7e2f','#800080'),lwd = 2),
                                         labels = c('C1','C2','C3','C4'),
                                         labels_gp = gpar(col = "white", fontsize = 12,just = "center",fontface = 'bold')),
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 12),
                                                     title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                     ncol=1),
                        #gap=unit(c(1, rep(0, 4), 1, rep(0, 5)), "mm"),
                        border = T,
                        col=list(Group=Group))

table(dd2$cluster)

col = c(rep(pal_npg("nrc", alpha =0.9)(9)[4],10),
        rep(pal_npg("nrc", alpha =0.9)(10)[9],9),rep(pal_d3("category20",alpha = 0.9)(20)[2],8))
col_title <- c(pal_npg("nrc", alpha =0.9)(9)[4],
               pal_npg("nrc", alpha =0.9)(10)[9],pal_d3("category20",alpha = 0.9)(20)[2])

pdf(file = 'Results/10.1 TIME/ICP_heatmap.pdf',width=7,height = 10)
Heatmap(t(dd2[,c(8:34)]),##对基因表达量进行限定，将>4-->4  <-4-->-4,也可以进行转置t后scale进行过滤
        name='Z-score',#热图的名字
        top_annotation = Top,#顶部注释  
        cluster_rows = FALSE,#不对行进行聚类
        col=colorRamp2(c(-2,0,2),c(pal_npg("nrc", alpha =0.7)(9)[2],'white',pal_npg("nrc", alpha =0.9)(9)[1])),#颜色
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        #show_row_names = FALSE,  
        row_names_gp = gpar(fontsize = 12,col = col),
        column_split = c(rep(1,118),rep(2,74),rep(3,40),rep(4,56)),
        row_split = factor(rep(c('B7-CD28','TNF superfamily','Others'),
                               times=c(10,9,8)),levels = c('B7-CD28','TNF superfamily','Others')),
        gap = unit(1.5, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 9),
        row_title_gp = gpar(fontsize=14,col = col_title),
        #width=unit(8, "cm"),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 12,colour = col), 
                                  title_gp = gpar(fontsize = 12 ,fontface = "bold"))) 
dev.off()

#箱线图
dd <- dd2[,c(1,8:34)]
dd2 <- pivot_longer(dd,cols = 2:28,names_to = 'ICPs',values_to = 'Value')

sort(colnames(dd)[2:11])
sort(colnames(dd)[12:20])
sort(colnames(dd)[21:28])

dd2$ICPs <- factor(dd2$ICPs,levels = c(sort(colnames(dd)[2:11]),sort(colnames(dd)[12:20]),sort(colnames(dd)[21:28])))
colnames(dd2)[1] <-c('Subtype') 
ggplot(dd2,aes(ICPs,Value,fill = Subtype))+
  geom_boxplot(outlier.colour = 'black',
               outlier.size = 0.1,
               outlier.alpha = 0.6,
               size=0.4)+
  stat_compare_means(label = 'p.signif',method = "anova",vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")))+
  theme_bw(base_rect_size = 1.5)+
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 10),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust = 1,colour = 'black'),
        axis.title.y = element_text(size=16))+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  #guides(fill=F)+
  labs(fill='Subtype',y='Relative Expression',x=NULL)
ggsave(filename = 'Results/10.1 TIME/ICP_box.pdf',height = 4,width = 10.2)


# ICTs-----------------------------
rm(list = ls())
library(data.table)
library(xlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(compareGroups)
library(reshape2)

mm <- read.xlsx2("Results/10.1 TIME/immune_related_molecular.xlsx",sheetIndex = 1)
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

TCGA_expr <- data
gene <- intersect(rownames(TCGA_expr),mm$GENE.SYMBLE)
mm <- mm[mm$GENE.SYMBLE %in% gene,]
sp <- intersect(colnames(TCGA_expr),dd$ID)
mm_df <- TCGA_expr[gene,sp] %>% 
  t()%>%
  as.data.frame() %>% 
  merge(dd,.,by.x = 1,by.y = 0) %>% 
  column_to_rownames("ID")
mm_df <- mm_df[,-c(2:18)]
colnames(mm_df)[1] <- 'Cluster'

#--比较两组的差异
mm_df$Cluster <- factor(mm_df$Cluster)
dd <- descrTable(Cluster~.,data = mm_df,method = 2)
print(dd)
pvals <- getResults(dd, "p.overall")
p.adjust(pvals, method = "BH")


#---绘图boxplot
gene <-  colnames(mm_df)[which(pvals<0.1)+1]
mm1 <- mm[mm$GENE.SYMBLE %in% gene,]
mm_df1 <- cbind(Cluster=mm_df$Cluster,mm_df[,gene])
ggdata <- melt(mm_df1,id.vars = "Cluster",
               variable.name = "Gene",value.name = "expr") %>% 
  merge(.,mm1,by.x = 2,by.y = 1,all = T) 
ggdata <- ggdata[order(ggdata$Type,ggdata$Gene),]
ggdata$Gene <- factor(ggdata$Gene)

ggdata_ca <- ggdata[ggdata$Type == "Co-stimulatory immune checkpoint targets",]
ggdata_ca <- ggdata_ca[-c(177:352,1409:1584,2289:2640),]
ggdata_ci <- ggdata[ggdata$Type == "Co-inhibitory immune checkpoint targets",]
ggdata_ci <- ggdata_ci[-c(1:176,881:1056,1409:1584,1585:1760,2289:2464,3345:3520),]
ggdata_mhc <- ggdata[ggdata$Type %in% c("MHC Class II" ,"MHC Class I"),]



ggplot(ggdata_ca,aes(Gene,expr,fill=Cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  theme_bw()+
  labs(title="Co-stimulatory ICTs",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
ggsave(filename = "Results/10.1 TIME/Co-stimulatory ICTs.pdf",width = 8,height = 10)

gg_ci <- ggplot(ggdata_ci,aes(Gene,expr,fill=Cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  theme_bw()+
  labs(title="Co-inhibitory ICTs",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
ggsave(gg_ci,filename = "Results/10.1 TIME/Co-inhibitory ICTs.pdf",width = 8,height = 10)

gg_mhc <- ggplot(ggdata_mhc,aes(Gene,expr,fill=Cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  theme_bw()+
  labs(title="MHC molecules",x=NULL,y="Expression")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30,hjust = 0.9),
        plot.title = element_text(size = 12,face = "bold",hjust = 0.3,vjust=0.8),
        legend.title = element_text(size = 12,face = "bold"))+
  stat_compare_means(aes(label=..p.signif..))
ggsave(gg_mhc,filename = "Results/10.1 TIME/MHC molecules.pdf",width = 8,height = 10)




# Immunophenoscore(IPS) --------------------------------------------------------------------------
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

## Assign colors 
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
gene_expression <- TCGA_expr
sample_names<-names(gene_expression)

## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
# For different 
IPSG<-read.table("Results/10.1 TIME/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL
# Gene names in expression file
GVEC<-row.names(gene_expression)
# Genes names in IPS genes file
VEC<-as.vector(IPSG$GENE)
# Match IPS genes with genes in expression file
ind<-which(is.na(match(VEC,GVEC)))
# List genes missing or differently named
MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
  cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {	
  GE<-gene_expression[[i]]
  mGE<-mean(GE)
  sGE<-sd(GE)
  Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])}

DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
write.table(DF,file="Results/10.1 TIME/IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")

colnames(TCGA_clin)[1] <- "SAMPLE"
cluster <- TCGA_clin
ips <- read.table("Results/10.1 TIME/IPS.txt",h=T)
ips <- ips %>% 
  dplyr::select(SAMPLE,IPS) %>% 
  merge(cluster, by.x=1, by.y=1)
ips[,2] <- scale(ips[,2])

colnames(ips)[2:3] <- c('IPS','cluster')
ips %>%
  ggplot(aes(cluster,IPS)) +
  geom_boxplot(aes(color=cluster), outlier.colour = NA, size=1, fill=NA) +
  # geom_jitter(aes(fill=cluster,color=cluster), width = 0.2, shape=21, size=3, alpha=0.7) +
  stat_compare_means(method = "kruskal.test",label = 'p.signif', size=9,vjust = 0,hjust = -2.5) +
  geom_violin(alpha=0.5,aes(fill=cluster),color=NA,trim = T)+
  theme_bw(base_rect_size = 2) +
  labs(x=NULL,y='Relative IPS', title = ' ') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.y = element_text(size = 16,colour = 'black',face="bold"),
        axis.text.x = element_text(size = 12,colour = 'black',face = "bold"),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.ticks = element_line(size = 1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.background = element_rect(fill='#f3f6f6')) +
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080')) +
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))
ggsave("Results/10.1 TIME/IPS.pdf", width = 5, height = 8)

# Antigen presentation score(APS) ---------------------------------
rm(list = ls())
library(xlsx)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)
#提取参考文献中APS基因
imsig <- read.xlsx("Results/10.1 TIME/NIHMS958212-supplement-7.xlsx", sheetIndex = 1,header = T)
imsig <- imsig[1:78,1:8]
ap <- imsig[imsig$Super.Category == "Antigen presentation",]
ap_gene <- ap$Gene

#16个APS基因的表达
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
apExpr <- data[rownames(data) %in% ap_gene,]

#合并分组信息
metadata <-  as.data.frame(t(apExpr)) %>% 
  merge(.,dd, by.x=0,by.y=1) %>% 
  column_to_rownames("Row.names")
metadata <- metadata[,-c(16:32)]
colnames(metadata)[15] <- 'Cluster'
metadata$Cluster <- factor(metadata$Cluster)
plotinput <- pivot_longer(metadata,cols = 1:14, names_to = "signature", values_to = "value")

mycolor3<-"#436b95"
mycolor4<-"#b81d25"
mycolor5 <- "#2f7e2f"
mycolor6<-"#800080"
cols <- c(mycolor3,mycolor4,mycolor5,mycolor6)
ggplot(plotinput,aes(signature,value,fill=Cluster))+
  geom_boxplot(width=0.9,position=position_dodge(0.9),
               outlier.colour = NA,color='grey30')+
  scale_fill_manual(values = cols)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
        legend.position ="top" ,#'none'
        axis.ticks.x = element_blank())+
  rotate_x_text(60)+
  ylab("Relative expression")+xlab("")+
  stat_compare_means(label = 'p.signif',label.y = 5,size=6)

#gsva计算APS评分
library(GSVA)
expr <- data
APS <- list('APS'= ap_gene)
GSVA <- gsva(expr = as.matrix(expr),
             gset.idx.list = APS,
             method="gsva")
APS <- data.frame(APS=t(GSVA))%>%rownames_to_column('ID')
save(APS,file = 'Results/10.1 TIME/APS.Rda')


#作图
rm(list = ls())
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
load('Results/10.1 TIME/APS.Rda')
dat <-  APS %>% 
  merge(.,dd, by.x=1,by.y=1) %>% 
  column_to_rownames("ID")
dat <- dat[,c(1,2)]
colnames(dat)[2] <- "Cluster"
dat$Cluster <- factor(dat$Cluster)
mycolor3<-"#436b95"
mycolor4<-"#b81d25"
mycolor5 <- "#2f7e2f"
mycolor6<-"#800080"
cols <- c(mycolor3,mycolor4,mycolor5,mycolor6)
ggplot(dat,aes_string("Cluster","APS"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(dat[,1])*1.05,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(dat[,1])*1.1)+
  theme_bw(base_rect_size = 2)+
  labs(y='Antigen presentation score',x=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
        legend.position ="top" ,#'none'
        axis.ticks.x = element_blank())
ggsave(filename = "Results/10.1 TIME/APS.pdf",width = 3.5,height = 4.5)

# immunogenicity associated indicators ----------------------------------
rm(list = ls())
library(ggplot2)
library(data.table)
library(tidyverse)
options(digits = 3)
#-读取immune landscape
load("Results/10.1 TIME/immune_landscape.RData")
ABSOLUTE_scores$ID <- substr(rownames(ABSOLUTE_scores),1,12)
ABSOLUTE_scores <- ABSOLUTE_scores[!duplicated(ABSOLUTE_scores$ID),-5]
rownames(ABSOLUTE_scores) <- substr(rownames(ABSOLUTE_scores),1,12)

BLCA_im <- immune_landscape[immune_landscape$TCGA.Study == "BLCA",]
BLCA_ab <- ABSOLUTE_scores[rownames(BLCA_im),]
BLCA_il<- merge(BLCA_im,BLCA_ab,by=0) %>%  column_to_rownames("Row.names")
rm(immune_landscape,ABSOLUTE_scores,BLCA_im,BLCA_ab)

load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
dd$ID <- substr(dd$ID,1,12)

#-合并分组信息
metadata <- merge(dd,BLCA_il,by.x = 1,by.y=0) %>% 
  column_to_rownames("ID")
ggdata <- metadata[,-c(2:21)] 
colnames(ggdata)[1] <- 'Cluster'
ggdata[,-c(1)] <- apply(ggdata[,-c(1)],2,as.numeric) %>% 
  as.data.frame()
ggdata$Cluster <- factor(ggdata$Cluster)

#-初步比较两组
library(compareGroups)
dd <- descrTable(Cluster~.,data = ggdata)
summary(dd)
print(dd)
pvals <- getResults(dd, "p.overall")
p.adjust(pvals, method = "BH")
export2csv(dd,file = "Results/10.1 TIME/immune_landscape.csv")
gg <- read.csv( "Results/10.1 TIME/immune_landscape.csv",header = T,skip = 1,check.names = F,row.names = 1)
colnames(gg) <- c("C1","C2","c3","c4",'p')

colnames(ggdata)
img <- c("Cluster",
         "SNV.Neoantigens","Indel.Neoantigens","CTA.Score",
         
         "Intratumor.Heterogeneity",  
         "Number.of.Segments","Fraction.Altered",
         "LOH_n_seg","LOH_frac_altered",
         
         "Homologous.Recombination.Defects","Aneuploidy.Score",
         
         "TCR.Richness","TCR.Shannon")

dat <- ggdata[,img] 
rownames(dat) <- rownames(metadata)
write.csv(dat,file = "13.Immune/immunogenicity associated indicators.csv")


df <- data.frame()
for (i in 2:13) {
  m1= mean(dat[dat$Cluster=="C1",i][!is.na(dat[dat$Cluster=="C1",i])])/mean(dat[,i][!is.na(dat[,i])]);m1=log2(m1)
  m2= mean(dat[dat$Cluster=="C2",i][!is.na(dat[dat$Cluster=="C2",i])])/mean(dat[,i][!is.na(dat[,i])]);m2=log2(m2)
  m3= mean(dat[dat$Cluster=="C3",i][!is.na(dat[dat$Cluster=="C3",i])])/mean(dat[,i][!is.na(dat[,i])]);m3=log2(m3)
  m4= mean(dat[dat$Cluster=="C4",i][!is.na(dat[dat$Cluster=="C4",i])])/mean(dat[,i][!is.na(dat[,i])]);m4=log2(m4)
  df <- rbind(df,data.frame("C1"=m1,"C2"=m2,"C3"=m3,"C4"=m4))
}
df <- df %>% as.matrix()

rownames(df) <- c("SNV Neoantigens","Indel Neoantigens","CTA Score",
                  
                  "Intratumor Heterogeneity",  
                  "Number of Segments","Fraction Altered",
                  "LOH_n_seg","LOH_frac_altered",
                  
                  "Homologous Recombination Defects","Aneuploidy Score",
                  
                  "TCR Richness","TCR Shannon")

colnames(df) <- c("C1","C2","C3","C4")

library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(export)
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", df[i, j]),
            x, y, gp = gpar(fontsize = 10.5))
}
top <- HeatmapAnnotation(Cluster=c("C1","C2","C3","C4"),
                         border = T,
                         annotation_name_side = 'right',
                         annotation_name_gp = gpar(fontsize=10),
                         show_annotation_name = F,
                         na_col = NA,
                         show_legend = F,
                         gap = unit(1, "mm"),
                         col = list(Cluster=c('C1'="#436b95",'C2'="#b81d25","C3"="#2f7e2f","C4"="#800080")))
plo <- Heatmap(df,
               #cell_fun = cell_fun,
               border = T,
               heatmap_width = unit(10, "cm"),
               heatmap_height = unit(1, "cm")*nrow(df),
               #top_annotation = top,
               #right_annotation = right,
               col = colorRamp2(c(-1,0,1),colors = c("#009698","white","#cd5c5c")),
               #column_split = 1:2,
               row_split = rep(c(1,2,3,4),time= c(3,5,2,2)),
               gap = unit(2, "mm"),
               rect_gp = gpar(color='grey',lwd=1.2,size=0.4),
               
               show_column_names = F,
               cluster_rows = F,cluster_columns = F,
               
               row_names_side = 'right',
               row_names_gp = gpar(fontface='bold',color='balck',fontsize = 12),
               
               heatmap_legend_param = list(title= "log2(Ratio)",
                                           at=c(-1,0,1),
                                           title_position = "topcenter",
                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                           legend_direction="horizontal"))
plo
library(ComplexHeatmap)
draw(plo,heatmap_legend_side = "bottom")
graph2pdf(file='Results/10.1 TIME/immunogenicity associated indicators.pdf',height=7,width=5)



# Radar map -----------------------------------------------------------------------
rm(list = ls())
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(tibble)
library(org.Hs.eg.db)
library(reshape2)
library(plyr)
library(GSVA)
library(circlize)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(data.table)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
immunogram <- fread('Results/10.1 TIME/immunogram geneset.txt', data.table = F)
immunogram <- immunogram[order(immunogram$Axis),3:4]%>%column_to_rownames('Axis')
immunogram$`Required marker genes` <- gsub(' ','',immunogram$`Required marker genes`)
a <- list()
for (i in 1:nrow(immunogram)) {
  a[[i]] <- unlist(strsplit(immunogram[i,1],','))
}
names(a) <- rownames(immunogram)
immunogram <- a
rm(a)

load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
tcga_gsva <- as.data.frame(t(gsva(as.matrix(data), immunogram, method = "ssgsea")))
str(tcga_gsva)
TCGA_cluster <- dd[,c(1,2)]
data <- merge(tcga_gsva,TCGA_cluster,by.x=0,by.y=1)
colnames(data)[ncol(data)] <- "Cluster"
radarinput <- arrange(data,Cluster)
radarinput <- radarinput[order(radarinput$Cluster),]
radarinput <-   column_to_rownames(radarinput, 'Row.names')
ml <- radarinput[, 1:10] #画图的ssGSEA数据
#minmax标准化
min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}
ml <- as.data.frame(apply(ml, 2, min.max.norm)) #scale标化

library(ggradar)
library(palmerpenguins)
library(tidyverse)
library(scales)
library(showtext)

font_add_google("Lobster Two", "lobstertwo")
font_add_google("Roboto", "roboto")
# Showtext will be automatically invoked when needed
showtext_auto()
data <- merge(TCGA_cluster,ml,by.x=1,by.y=0)%>%column_to_rownames('ID')
colnames(data) <- gsub(' ','_',colnames(data))
colnames(data)[c(1,3,7,8)] <- c('Cluster','Priming_and_activation','Inhibitory_cells_Tregs','Inhibitory_cells_MDSCs')
data$Cluster <- as.factor(data$Cluster)
data <- data[order(data$Cluster),]
table(data$Cluster)
data <- data[,-1]
radar <- data.frame(group = c('C1','C2','C3',"C4"))
for (i in colnames(data)) {
  radar <- cbind(radar, i = c(mean(data[1:119,i]),mean(data[120:193,i]),mean(data[194:233,i]),mean(data[234:289,i])))
}
colnames(radar) <- c('group',colnames(data))
#整体矫正而非每列矫正
radar[,-1] <- min.max.norm(radar[,-1])
range(radar[,-1])
###画图
plt <- radar %>%
  ggradar(
    font.radar = "roboto",
    grid.label.size = 4,  # Affects the grid annotations (0%, 50%, etc.)
    axis.label.size = 4, # Afftects the names of the variables
    group.point.size = 4,   # Simply the size of the point
    group.colours = c('#436b95','#b81d25','#2f7e2f','#800080')
  )
plt
# 1. Set the position legend to bottom-right
# 2. Bottom-right justification
# 3. Customize text size and family
# 4. Remove background and border color for the keys
# 5. Remove legend background
plt <- plt +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.text = element_text(size = 14, family = "roboto"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.background = element_blank()
  )
plt
# * The panel is the drawing region, contained within the plot region.
#   panel.background refers to the plotting area
#   plot.background refers to the entire plot
plt <- plt +
  labs(title = "Radar plot of penguins species") +
  theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.title.position = "plot", # slightly different from default
    plot.title = element_text(
      family = "Arial",
      size = 14,
      face = "bold",
      color = "#2a475e"
    )
  )
plt
ggsave(filename = "Results/10.1 TIME/radar.pdf",width = 14,height = 8)


# Stromal cell infiltration -------------------------------------------------
rm(list = ls())
library(ggcor)
library(parallel)
library(IOBR)
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

expr <- data[,dd$ID]
mcp <- deconvo_tme(eset = expr, method = "mcpcounter")

dd <- merge(dd[,c(1,2)],mcp,by=1)
colnames(dd)[1:2] <- c('sample','Subtype')
dd <- tibble::column_to_rownames(dd,"sample")
dd1 <- dd[,c(1,10,11)]
colnames(dd1) <- c('Subtype','Endothelial cells','Fibroblasts')
dd2 <- pivot_longer(dd1,cols = 2:3,names_to = "Scores",values_to = 'Value')

my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
library(ggsci)
ggplot(dd2,aes(Subtype,Value,color=Subtype))+
  geom_boxplot(outlier.colour = NA,aes(fill=Subtype),color='black',size=0.6,alpha=0.65)+
  #geom_jitter(width = 0.2,size=1.5,shape=21,aes(fill=cluster))+
  geom_violin(alpha=0.5,aes(fill=Subtype),color=NA,trim = T)+
  #stat_compare_means(label = 'p.signif',method = "anova",label.x = 1.95)+
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.12,#两组直接对比框的高度
              map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
              test = wilcox.test,size=0.85,#线条粗细
              textsize = 3,#P值字号大小
              tip_length = 0.03,
              y_position =1.5, #P值位置
              show.legend= FALSE
  )+
  #scale_y_log10()+
  ylim(-2,4)+
  facet_wrap(~Scores,scale = "free")+
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+ 
  theme_bw(base_rect_size = 2)+
  labs(y='Scores')+
  guides(fill=F)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        #axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=15,colour = 'black'),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.background = element_rect(fill='#f3f6f6'))
ggsave(filename = 'Results/10.1 TIME/Stroma_cell.pdf',width = 4.2,height = 3.5)


# TME -----------------------------------------------
rm(list = ls())
library(ggcor)
library(parallel)
library(IOBR)
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')
load('Results/6. bulk_NTP&survival/logtpm3.rda')
expr <- data
cluster <- dd[,c(1,2)]
expr <- expr[,cluster$ID]
estimate <- function(dat,pro){
  input.f=paste0('Results/10.1 TIME/',pro,'_estimate_input.txt')
  output.f=paste0('Results/10.1 TIME/',pro,'_estimate_gene.gct')
  output.ds=paste0('Results/10.1 TIME/',pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- 'BLCA'
scores <- estimate(expr,pro)

TumorPurity <- cos(0.6049872018+0.0001467884 * scores[,3])
names(TumorPurity) <- colnames(expr)
scores <- as.data.frame(scores)
scores$TumorPurity <- TumorPurity
score <- scores[,-3]
score <- scale(score)
cluster$ID <- gsub('-','.',cluster$ID)
score <- merge(cluster,score,by.x=1,by.y=0)
rownames(score) <- score[,1]
score <- score[,-1]
dd2 <- pivot_longer(score,cols = 2:4,names_to = 'Scores',values_to = 'Value')
colnames(dd2)[1] <- c('Subtype')

my_comparisons <- list(c("C1","C2"), c("C1","C3"),c('C1','C4'),c('C2','C3'),c('C2','C4'),c('C3','C4'))
library(ggsci)
ggplot(dd2,aes(Subtype,Value,color=Subtype))+
  geom_boxplot(outlier.colour = NA,aes(fill=Subtype),color='black',size=0.6,alpha=0.65)+
  #geom_jitter(width = 0.2,size=1.5,shape=21,aes(fill=cluster))+
  geom_violin(alpha=0.5,aes(fill=Subtype),color=NA,trim = T)+
  #stat_compare_means(label = 'p.signif',method = "anova",label.x = 2, vjust=0.4,size=5)+
  #scale_y_log10()+
  geom_signif(comparisons = my_comparisons,step_increase = 0.1,#两组直接对比框的高度
              map_signif_level = T,#这里可以设置两组之间差异显示为数值或者***，F为数值，T为***
              test = wilcox.test,size=0.9,#线条粗细
              textsize = 2.5,#P值字号大小
              y_position = 2, #位置
              show.legend= FALSE
  )+
  ylim(-3,4.2)+
  facet_wrap(~Scores,scales = 'free')+
  scale_color_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+
  scale_fill_manual(values = c('#436b95','#b81d25','#2f7e2f','#800080'))+ 
  theme_bw(base_rect_size = 2)+
  guides(fill=F)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=15,colour = 'black'),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.background = element_rect(fill='#f3f6f6'))
ggsave(filename = 'Results/10.1 TIME/TME.pdf',width = 8,height = 3.5)


# TIS --------------------------------
rm(list = ls())
options(stringsAsFactors = F)
library(MCPcounter)
library(GSVA)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggsci)
load("Results/6. bulk_NTP&survival/TCGA_cluster.rda")
expr <- data
rm(TCGA_expr)
table(substr(colnames(expr),14,16))


TISG <- c("CCL5", "CD27", "CD274", "CD276",
          "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1",
          "HLA-DRB1", "HLA-E", "IDO1",
          "LAG3", "NKG7", "PDCD1LG2",
          "PSMB10", "STAT1", "TIGIT")
GSVA <- gsva(expr=as.matrix(expr),
             gset.idx.list=list(TISG), method="gsva")
TIS <- data.frame(TIS=t(GSVA))%>%rownames_to_column('ID')

#作图
dat <-  TIS %>%
  merge(.,dd[,c(1,2)], by.x=1,by.y=1) %>%
  column_to_rownames("ID")
colnames(dat)[2] <- 'Cluster'
dat$Cluster <- factor(dat$Cluster)
mycolor1<-"#436b95"
mycolor2<-"#b81d25"
mycolor3<-"#2f7e2f"
mycolor4<-"#800080"

#作图
cols <- c(mycolor1,mycolor2,mycolor3,mycolor4)
ggplot(dat,aes_string("Cluster","TIS"))+
  geom_jitter(shape = 21,size=2,width = 0.2,aes_string(fill="Cluster",color="Cluster"))+
  geom_boxplot(outlier.colour = NA,aes_string(fill="Cluster"),color='black',size=0.6,alpha=0.65)+
  geom_violin(alpha=0.5,aes_string(fill="Cluster"),color=NA,trim = T)+
  stat_compare_means(label.y = max(dat[,1])*1.05,color='black',size=5,hjust = 0.1)+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  expand_limits(y=max(dat[,1])*1.1)+
  theme_bw(base_rect_size = 2)+
  labs(y='Tumor Inflammation Signature',x=NULL)+
  theme(axis.text.y = element_text(size = 11, colour = 'black'),
        axis.text.x = element_text(size = 14,colour = 'black',angle = 0),
        axis.title.x = element_text(size = 15,colour = 'darkred',face='bold'),
        axis.title.y = element_text(size = 15,colour = 'darkred',face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        plot.title = element_text(hjust = 0.5,size = 15,colour = 'darkred',face='bold'),
        legend.position = 'top',#'none'
        axis.ticks.x = element_blank())
ggsave(filename = "Results/10.1 TIME/TIS.pdf",width = 3.5,height = 4.5)




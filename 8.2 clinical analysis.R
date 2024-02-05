
#bulk亚型的临床相关性热图-------------------------------------
rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(cowplot)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

TCGA_expr <- data
TCGA_clin <- dd[,c(1,2,11:18)]
clin <- TCGA_clin[,c("ID","prediction","Age","Gender","Stage","T","N","M",'Grade','Subtype')]
colnames(clin) <- c("ID","Cluster","Age","Gender","Stage","Ts","Ns","Ms",'Grade','Subtype')
clin <- clin[order(clin$Cluster),]


Cluster <- c('#436b95','#b81d25','#2f7e2f','#800080')
names(Cluster) <- c('C1','C2','C3','C4')
Age <- c(pal_npg("nrc", alpha =0.7)(9)[3:4])
names(Age) <- c("≤65",">65")
Gender <-c(pal_d3("category20",alpha = 0.7)(20)[16:17])
names(Gender) <- c("Male","Female")
Stage <- c(pal_d3("category20c",alpha = 0.7)(20)[c(17,12,7,2)])
names(Stage) <- c('I','II','III','IV')
Grade <- c('#ffe6e6','#ff9999')  
names(Grade) <- c('High','Low')
Ts <-  c(pal_d3("category20c",alpha = 0.7)(20)[c(18,13,8,3)])
names(Ts) <- c('T1','T2','T3','T4')
Ns <- c(pal_d3("category20c",alpha = 0.7)(20)[c(6,1,9,13)])
names(Ns) <- c('N0','N1','N2','N3')
Ms <- c(pal_d3("category20b",alpha = 0.7)(20)[c(10,5)])
names(Ms) <- c('M0','M1')
Subtype <- c(pal_d3("category20",alpha = 0.7)(20)[c(15,5)])
names(Subtype)<- c('Non-Papillary','Papillary')


e <- NA
names(e) <- "Cluster"
dd <- c(e)
for (i in colnames(clin)[2:c(ncol(clin)-1)]) {
  a <- xtabs(~clin[,i]+Cluster,data=clin)
  b <- chisq.test(a)
  p.value <- b$p.value
  names(p.value) <- i
  dd <- c(dd,p.value)
}
DD <- as.data.frame(dd)
DD$dd <- round(DD$dd,4)
DD$dd <- ifelse(DD$dd<0.05,'<0.05',DD$dd)
DD$dd[as.numeric(DD$dd) > 0.05] <- "ns"
DD$dd[DD$dd < 0.05 & DD$dd > 0.01] <- "*"
DD$dd[DD$dd < 0.01 & DD$dd > 0.001] <- "**"
DD$dd[DD$dd < 0.001] <- "***"
DD$dd[DD$dd %in% 0.001] <- "***"

DD[1,1] <- ""
dd1 <- data.frame(row.names = colnames(clin)[2:10],value=DD$dd)


dd1 <- dd1[2:9,,drop=F]


ha <-  HeatmapAnnotation(Cluster = anno_block(gp = gpar(fill = c('#436b95','#b81d25','#2f7e2f','#800080'),lwd = 2),
                                              labels = c('C1','C2','C3','C4'),
                                              labels_gp = gpar(col = "white", fontsize = 11,just = "center")),
                         Age = clin$Age, Gender = clin$Gender,Stage = clin$Stage,Grade = clin$Grade,Ts = clin$Ts,Ns = clin$Ns,Ms= clin$Ms,Subtype=clin$Subtype,na_col = 'grey90',
                         annotation_legend_param = list(
                           Cluster = list(title = "Cluster"),
                           Age = list(title = "Age"),
                           Gender = list(title = "Gender"),
                           Stage = list(title = "Stage"),
                           Grade = list(title = 'Grade'),
                           Ts = list(title = 'Ts'), 
                           Ns = list(title = 'Ns'), 
                           Ms = list(title = 'Ms'), 
                           Subtype = list(title = 'Subtype'),
                           labels_gp = gpar(fontsize = 10),
                           title_gp = gpar(fontsize = 12),
                           ncol=1),#annotation legend的标签
                         border = T,
                         gp = gpar(lwd = 1),
                         gap = unit(c(1,rep(0,9),1),'mm'),
                         col = list(Cluster = Cluster, Age = Age, Gender = Gender, Stage = Stage,Grade = Grade,Ts = Ts, Ns = Ns, Ms = Ms,Subtype = Subtype),
                         show_annotation_name = TRUE,
                         annotation_name_gp = gpar(fontsize = 12),
                         annotation_name_side="left",
                         annotation_height = unit(rep(5,9), "mm")
                         #各个临床信息的颜色，Anatomic_location颜色输入也同其他，更好颜色搭配信息请参考FigureYa28color
)
col_fun <- colorRamp2(c(-2, 0, 2),c(pal_npg("nrc", alpha =0.7)(9)[2], "white", pal_npg("nrc", alpha =0.9)(9)[1]))
hm <- t(data.frame(row.names = clin$ID,  d_e_l_e_t = rep(0, nrow(clin))))
table(clin$Cluster)
ht <- Heatmap(hm,name = "Z-score", col = col_fun,
              height = 0,
              cluster_rows = TRUE , cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              top_annotation = ha, 
              column_title = " ",
              border_gp = gpar(col = "black", lwd = 4),
              border = F,
              column_split = c(rep(1,119),rep(2,74),rep(3,40),rep(4,56)),
              show_heatmap_legend = F
)


annotation_titles <- data.frame(row.names  = c(Age = 'Age', Gender = 'Gender',  #Cluster =' Cluster', 
                                              Stage = 'Stage',Grade= 'Grade',Ts = 'Ts', Ns = 'Ns', Ms = 'Ms',Subtype = 'Subtype'),
                                P = c('***','***','ns','ns','ns','ns','ns','***'))   #'P value',

par(mar = c(5, 10, 5, 10)) 
pdf('Results/8. COX&clinical analysis/clin_Heatmap.pdf',6,9)
draw(ht,annotation_legend_side = "bottom",heatmap_legend_side = "left")
for(an in rownames(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an,'P'], unit(15, "cm"), just = "left",gp = gpar(fontsize = 12))#对齐方向是右边
  })
}
dev.off()

rownames(dd1) <- c('Age','Gender','Stage','Grade','Ts','Ns', 'Ms','Subtype')
save(dd1,file = 'Results/8. COX&clinical analysis/clin.rda')

#临床性状百分比图---------------------------------------------------------------
#输入文件准备
mypalette <- c(alpha("#003153",0.8),alpha('#008080',0.8),"#c0c0c0",alpha("darkred",0.8))
plotdata <- clin
x0 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Age = plotdata$Age),length)
x1 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Gender = plotdata$Gender),length)
x2 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Stage = plotdata$Stage),length)
x3 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Grade = plotdata$Grade),length)
x4 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           T = plotdata$T),length)
x5 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           M = plotdata$M),length)
x6 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           N = plotdata$N),length)
x7 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Growth = plotdata$Subtype),length)
#Growth Pattern-----------------------------
data2 <- x7
p <- fisher.test(matrix(as.numeric(data2$x),nrow = 2))$p.value
ct12 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C2'),]$x),nrow = 2))$p.value
ct23 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C3'),]$x),nrow = 2))$p.value
ct34 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C3','C4'),]$x),nrow = 2))$p.value
ct14 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C4'),]$x),nrow = 2))$p.value
ct13 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C3'),]$x),nrow = 2))$p.value
ct24 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C4'),]$x),nrow = 2))$p.value
data2$percent <- c(77/117,61/73,26/40,28/56,40/117,12/73,14/40,28/56)
data2$label <- round(data2$percent*100,1)

ggplot(data2, aes( x = Cluster, y = percent, fill = Growth))+
  geom_bar( stat = 'identity' ,position = position_stack())+
  geom_text(aes(label = label),position = position_stack(vjust = 0.5),size=5,color="black")+
  labs(title = 'Growth Pattern')+
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
  annotate("text",x = c(1.5),y = -0.05,label = c("**"),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(2.5),y = -0.05,label = c('*'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3.5),y = -0.05,label = c('ns'),size =6,fontface = "bold")+
  annotate("text",x = c(2),y = 1.03,label = c('ns'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3),y = 1.09,label = c('***'),size =6,fontface = "bold")+
  annotate("text",x = c(2.5),y = 1.15,label = c('ns'),size =6,fontface = "bold")+
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
ggsave(filename = 'Results/8. COX&clinical analysis/clin_Growth.pdf',width = 4.7,height = 9)

#Stage-----------------------------
clin$Stage <- ifelse(clin$Stage %in% c('I','II'),'I&II',
                     ifelse(clin$Stage %in% 'III','III','IV'))
table(clin$Stage)
plotdata <- clin
x2 <- aggregate(plotdata$Cluster,by = list(Cluster = plotdata$Cluster,
                                           Stage = plotdata$Stage),length)
data2 <- x2
p <- fisher.test(matrix(as.numeric(data2$x),nrow = 2),simulate.p.value=TRUE)$p.value
ct12 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C2'),]$x),nrow = 2))$p.value
ct23 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C3'),]$x),nrow = 2))$p.value
ct34 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C3','C4'),]$x),nrow = 2))$p.value
ct14 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C4'),]$x),nrow = 2))$p.value
ct13 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C3'),]$x),nrow = 2))$p.value
ct24 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C4'),]$x),nrow = 2))$p.value
data2$percent <- c(38/119,19/74,15/40,25/56,
                   31/119,30/74,15/40,18/56,
                   50/119,25/74,10/40,13/56)
data2$label <- round(data2$percent*100,1)

ggplot(data2, aes( x = Cluster, y = percent, fill = Stage))+
  geom_bar( stat = 'identity' ,position = position_stack())+
  geom_text(aes(label = label),position = position_stack(vjust = 0.5),size=5,color="black")+
  labs(title = 'Stage')+
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
  annotate("text",x = c(2),y = 1.03,label = c('ns'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3),y = 1.09,label = c('ns'),size =6,fontface = "bold")+
  annotate("text",x = c(2.5),y = 1.15,label = c('*'),size =6,fontface = "bold")+
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
ggsave(filename = 'Results/8. COX&clinical analysis/clin_Stage.pdf',width = 4.7,height = 9)


#Grade-----------------------------
data2 <- x3
p <- fisher.test(matrix(as.numeric(data2$x),nrow = 2))$p.value
ct12 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C2'),]$x),nrow = 2))$p.value
ct23 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C3'),]$x),nrow = 2))$p.value
ct34 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C3','C4'),]$x),nrow = 2))$p.value
ct14 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C4'),]$x),nrow = 2))$p.value
ct13 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C1','C3'),]$x),nrow = 2))$p.value
ct24 <- fisher.test(matrix(as.numeric(data2[data2$Cluster %in% c('C2','C4'),]$x),nrow = 2))$p.value
data2$percent <- c(115/117,74/74,37/40,44/56,2/117,3/40,12/56)
data2$label <- round(data2$percent*100,1)

ggplot(data2, aes( x = Cluster, y = percent, fill = Grade))+
  geom_bar( stat = 'identity' ,position = position_stack())+
  geom_text(aes(label = label),position = position_stack(vjust = 0.5),size=5,color="black")+
  labs(title = 'Grade')+
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
  annotate("text",x = c(1.5),y = -0.05,label = c("***"),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(2.5),y = -0.05,label = c('***'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3.5),y = -0.05,label = c('ns'),size =6,fontface = "bold")+
  annotate("text",x = c(2),y = 1.03,label = c('ns'),size =5,fontface = "bold",vjust  =0.1)+
  annotate("text",x = c(3),y = 1.09,label = c('***'),size =6,fontface = "bold")+
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
ggsave(filename = 'Results/8. COX&clinical analysis/clin_Grade.pdf',width = 4.7,height = 9)



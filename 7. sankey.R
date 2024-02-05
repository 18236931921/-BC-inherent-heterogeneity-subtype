
#NMIBC_sankey---------------------------------------------------------
rm(list=ls())
#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(classifyNMIBC)
load('Datasets/Bladder cancer.rda')
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

data <- total_expr_list[['TCGA_BLCA']][,colnames(data)]
range(data)
NMIBC_class <- classifyNMIBC(data)

mtcars <- merge(NMIBC_class,dd,by.x = 0,by.y = 1)
mtcars <- mtcars[,c(1,2,9)]
#mtcars <- mtcars%>%tibble::column_to_rownames('Row.names')
colnames(mtcars) <- c('sample','NMIBC_class','Pathway_based_BLCA_classification')  

df <- mtcars %>%
  make_long(NMIBC_class, Pathway_based_BLCA_classification)

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/NMIBC_sankey.pdf',width = 8,height = 10)

#六种转录组学BC_sankey---------------------------------------------------------
# require(devtools)
# devtools::install_github("cit-bioinfo/BLCAsubtyping", build_vignettes = TRUE)
#devtools::install_github("davidsjoberg/ggsankey")
rm(list=ls())
library(ggsankey)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(BLCAsubtyping)

load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

data <- as.matrix(data)
BC_class <- classify(data)


#Baylor.subtype
table(BC_class$Baylor.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,2,8)]
colnames(mtcars) <- c('sample','Baylor.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(Baylor.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/Baylor.subtype_sankey.pdf',width = 8,height = 10)

#UNC.subtype
table(BC_class$UNC.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,3,8)]
colnames(mtcars) <- c('sample','UNC.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(UNC.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/UNC.subtype_sankey.pdf',width = 8,height = 10)

#CIT.subtype
table(BC_class$CIT.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,4,8)]
colnames(mtcars) <- c('sample','CIT.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(CIT.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/CIT.subtype_sankey.pdf',width = 8,height = 10)

#Lund.subtype
table(BC_class$Lund.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,5,8)]
colnames(mtcars) <- c('sample','Lund.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(Lund.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/Lund.subtype_sankey.pdf',width = 8,height = 10)

#MDA.subtype
table(BC_class$MDA.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,6,8)]
colnames(mtcars) <- c('sample','MDA.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(MDA.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/MDA.subtype_sankey.pdf',width = 8,height = 10)

#TCGA.subtype
table(BC_class$TCGA.subtype)
mtcars <- merge(BC_class,dd,by= 1)
mtcars <- mtcars[,c(1,7,8)]
colnames(mtcars) <- c('sample','TCGA.subtype','Pathway_based_BLCA_classification')  
df <- mtcars %>%
  make_long(TCGA.subtype, Pathway_based_BLCA_classification)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/TCGA.subtype_sankey.pdf',width = 8,height = 10)


#MIBC_sankey---------------------------------------------------
# library(devtools)
# devtools::install_github("cit-bioinfo/consensusMIBC", build_vignettes = TRUE)
rm(list=ls())
library(ggsankey)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(consensusMIBC)
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

MIBC_class <- getConsensusClass(data,gene_id = c("entrezgene","ensembl_gene_id", "hgnc_symbol")[3])

mtcars <- merge(MIBC_class,dd,by.x = 0,by.y = 1)
mtcars <- mtcars[,c(1,2,11)]
colnames(mtcars) <- c('sample','MIBC_class','Pathway_based_BLCA_classification')  

df <- mtcars %>%
  make_long(MIBC_class, Pathway_based_BLCA_classification)

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
ggsave(filename = 'Results/7. sankey/MIBC_sankey.pdf',width = 8,height = 10)


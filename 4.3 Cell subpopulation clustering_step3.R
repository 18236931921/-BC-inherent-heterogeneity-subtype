rm(list=ls())
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(ConsensusClusterPlus)
library(heatmap3)
library(diceR)
library(clue)

ansMwwGST <- NULL
files <- list.files('Results/4.2 Cell subpopulation clustering_step2/Rdata2/')
for (i in 1:length(files)){
  load(paste0('Results/4.2 Cell subpopulation clustering_step2/Rdata2/',files[i]))
  names(aMwwGSTs) <- paste0(str_extract(files,"BC[0-9]+")[i],'_',names(aMwwGSTs))
  ansMwwGST <- c(ansMwwGST,aMwwGSTs)
}

source('4.3 SCC.step3.R')
SCC.step3(ansMwwGST,bestK.method = 'elbow', data.path = 'Results/4.3 Cell subpopulation clustering_step3/Rdata3',
          figure.path = 'Results/4.3 Cell subpopulation clustering_step3/Figure3')   #bestK需要自己根据聚类结果去挑

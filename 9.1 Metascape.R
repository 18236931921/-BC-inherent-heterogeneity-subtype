#-------------------------------------------------------------------------------
rm(list = ls())
load("Results/6. bulk_NTP&survival/markers.rda")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(base)

#Metascape
BCIS1 <- as.data.frame(markers[markers$class%in%"BCIS1",2])
BCIS2 <- as.data.frame(markers[markers$class%in%"BCIS2",2])
BCIS3 <- as.data.frame(markers[markers$class%in%"BCIS3",2])
BCIS4 <- as.data.frame(markers[markers$class%in%"BCIS4",2])
write.table(BCIS1,file = 'Results/9. Function/Metascape/BCIS1_markers.txt',quote = F,row.names = F)
write.table(BCIS2,file = 'Results/9. Function/Metascape/BCIS2_markers.txt',quote = F,row.names = F)
write.table(BCIS3,file = 'Results/9. Function/Metascape/BCIS3_markers.txt',quote = F,row.names = F)
write.table(BCIS4,file = 'Results/9. Function/Metascape/BCIS4_markers.txt',quote = F,row.names = F)




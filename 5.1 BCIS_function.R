rm(list = ls())
source('5.2 DEAgroups.R')
load("Results/4.3 Cell subpopulation clustering_step3/Rdata3/aConsensus_allClusters.rda", verbose = T)
load("Results/4.3 Cell subpopulation clustering_step3/Rdata3/aMwwGSTs_allClusters.rda", verbose = T)

levels(consClust)[levels(consClust) == "black"] <- "cyan"
levels(consClust)
levels(consClust) <- c("BCIS1", "BCIS2", "BCIS3", "BCIS4")

tmp <- as.character(consClust)
names(tmp) <- names(consClust)
consClust <- tmp

consClust <- consClust[consHc$order][length(consClust):1]

aDEA <- DEAgroups(ddata = NES_allClusters, groups = consClust, method = "MWW")
aDEA2 <- lapply(aDEA, function(x){
  x <- x[x$logFC > 0.5 & x$pValue < 0.01 , ]
  x <- x[order(x$pValue), ]
  x <- x[order(x$logFC, decreasing = T), ]
  return(x)
})
aDEA2 <- lapply(aDEA2, rownames)

NES_allClusters <- NES_allClusters[, names(consClust)]
save(NES_allClusters,aDEA2,file = 'Results/5. BCIS_function/NES_allcluster.rda') #为绘图做准备


#输出为xls格式----------------------
allPath <- as.character(unlist(aDEA2))
dup <- sort(unique(allPath[duplicated(allPath)]))
allPath <- allPath[!allPath %in% dup]
aDEA2 <- lapply(aDEA2, function(x) x[x %in% allPath])

outPath <- Reduce(c, aDEA2)
outPath <- data.frame(path = outPath, cluster = rep(names(unlist(lapply(aDEA2, length))), as.numeric(unlist(lapply(aDEA2, length)))))
outPath$pp <- gsub('_',' ', outPath$path)
outPath$pp <- tolower(outPath$pp)
table(outPath$cluster)
write.csv(outPath, file = 'Results/5. BCIS_function/cluster_outPath.csv')


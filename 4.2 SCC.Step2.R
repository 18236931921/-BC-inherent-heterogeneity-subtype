##### Step 2 Function #####


SCC.step2 <- function(dataList, pathways, data.path = 'Rdata/SCC.step2'){
  if(!file.exists(data.path)){dir.create(data.path,recursive = TRUE)}
  
  library(foreach)
  library(doParallel)
  
  
  for(i in 1:length(dataList)){ 
    G <- dataList[[i]]$bestK
    aConsensus <- dataList[[i]]$consRES[[G]]
    geData <- dataList[[i]]$geData
    obj <- aConsensus
    
    consHc <- obj$consensusTree
    consClust <- as.factor(cutree(consHc, k = G))
    names(consClust) <- colnames(geData)
    
    groups <- consClust
    levels(groups) <- paste0("c", 1:length(levels(groups)))
    
    consMatrix <- 1 - aConsensus$consensusMatrix
    rownames(consMatrix) <- colnames(consMatrix) <- colnames(geData)
    consMatrix <- as.dist(consMatrix)
    tmp <- as.numeric(consClust)
    names(tmp) <- names(consClust)
    
    library(cluster)
    aSil <- silhouette(x = tmp, dist = consMatrix)
    # plot(aSil, main = tumorID, col = levels(consClust))
    # abline(v=0.5)
    
    groups <- groups[aSil[, 3] > 0.5]
    
    nGG <- table(groups)
    toTake <- names(which(nGG > 10))
    if(length(toTake) == 0) next
    
    rankedLists <- vector("list", length(toTake))
    names(rankedLists) <- toTake
    for(gg in 1:length(rankedLists)){
      whichOfInterest <- names(groups)[groups == names(rankedLists)[gg]]
      theOthers <- setdiff(names(groups), whichOfInterest)
      
      ans <- apply(geData, 1, function(x) wilcox.test(x[whichOfInterest], x[theOthers]))
      rankedList <- unlist(sapply(ans, function(x) x$statistic))/length(whichOfInterest)/length(theOthers) 
      names(rankedList) <- gsub("\\.W", "", names(rankedList))
      rankedList <- log2(rankedList/(1-rankedList))
      rankedList <- sort(rankedList, decreasing = TRUE)
      
      
      rankedLists[[gg]] <- rankedList
      # print(gg)
    }
    
    library(yaGST)
    aMwwGSTs <- vector("list", length(rankedLists))
    names(aMwwGSTs) <- names(rankedLists)
    for(gg in 1:length(rankedLists)){
      rankedList <- rankedLists[[gg]]
      
      aMwwGST <- lapply(pathways, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
      aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
      
      originalGeneSetCount <- sapply(aMwwGST, function(x) x$originalGeneSetCount)
      actualGeneSetCount <- sapply(aMwwGST, function(x) x$actualGeneSetCount)
      NES <- sapply(aMwwGST, function(x) x$nes)
      odd_NES <- sapply(aMwwGST, function(x) x$pu)
      logit2_NES <- sapply(aMwwGST, function(x) x$log.pu)
      pValue <- sapply(aMwwGST, function(x) x$p.value)
      qValue <- p.adjust(pValue, method = "BH")
      
      aMwwGST <- cbind(originalGeneSetCount, actualGeneSetCount, NES, odd_NES, logit2_NES, pValue, qValue)
      aMwwGST <- as.data.frame(aMwwGST)
      
      aMwwGSTs[[gg]] <- aMwwGST
      # print(gg)
    }
    ffile <- paste0(data.path, "/rankedLists_aMwwGSTs_", names(dataList)[i], ".rda")
    save(rankedLists, aMwwGSTs, file = ffile)
    
    names(aMwwGSTs) <- paste0(names(dataList)[i], "_", names(aMwwGSTs))
    
    print(paste0("TumorID :", names(dataList)[i]))
  }
}

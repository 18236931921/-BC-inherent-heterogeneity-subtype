##### Step 3 Function #####

SCC.step3 <- function(ansMwwGST, bestK.method = 'PAC', data.path = 'Rdata/SCC.step3',
                      figure.path = 'Figures/SCC.step3'){
  if(!file.exists(data.path)){dir.create(data.path,recursive = TRUE)}
  if(!file.exists(figure.path)){dir.create(figure.path,recursive = TRUE)}
  
  rowNames <- sort(unique(unlist(lapply(ansMwwGST, rownames))))
  colNames <- names(ansMwwGST)
  
  NES_allClusters <- matrix(NA, nrow = length(rowNames), ncol = length(colNames))
  rownames(NES_allClusters) <- rowNames
  colnames(NES_allClusters) <- colNames
  pValue_allClusters <- FDR_allClusters <- NES_allClusters
  for(j in 1:length(ansMwwGST)){
    NES_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$logit2_NES
    pValue_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$pValue
    FDR_allClusters[rownames(ansMwwGST[[j]]), j] <- ansMwwGST[[j]]$qValue
  }
  
  pValue_allClusters <- pValue_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  FDR_allClusters <- FDR_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  NES_allClusters <- NES_allClusters[rowSums(is.na(NES_allClusters)) == 0, ]
  
  ffile <- paste0(data.path, "/aMwwGSTs_allClusters.rda")
  save(NES_allClusters, pValue_allClusters, FDR_allClusters, file = ffile)
  NES <- NES_allClusters
  range(NES)
  
  #将NES转化为二分类矩阵-----------
  NES[NES_allClusters < 0.58] <- 0
  NES[FDR_allClusters > 0.001] <- 0
  #NES[pValue_allClusters > 0.001] <- 0
  NES <- NES[which(rowSums(NES > 0) != 0), ]
  NES <- NES[, which(colSums(NES > 0) != 0)]
  
  jaccard <- matrix(0, nrow = ncol(NES), ncol = ncol(NES))
  rownames(jaccard) <- colnames(jaccard) <- colnames(NES)
  for(i in 1:ncol(NES)){
    for(j in 1:ncol(NES)){
      a <- rownames(NES[NES[, i] != 0, i, drop = F])
      b <- rownames(NES[NES[, j] != 0, j, drop = F])
      jaccard[i, j] <- length(intersect(a, b))/length(union(a, b))
      #print(i)
    }
  }
  jaccard <- 1 - jaccard
  
  #聚类-----------------
  maxK <- 10
  results <- ConsensusClusterPlus(d = jaccard,    #,t(scale(t(NES))),as.dist(jaccard)
                                  maxK=maxK,     #K最大值
                                  reps=1000,    #重抽样次数
                                  pItem=0.9,   #样品抽样比例
                                  pFeature=1,   #Feature抽样比例
                                  title= figure.path,
                                  #innerLinkage="ward.D2",
                                  #finalLinkage="ward.D2",
                                  clusterAlg="pam",  #聚类算法,pam
                                  distance="euclidean",  #计算距离的方法,manhattan,euclidean
                                  seed=123456,
                                  plot="pdf")   #不设置，图片结果仅输出到屏幕
  #aConsensus <- ConsensusClusterPlus(d = as.dist(jaccard), maxK = ifelse(nCls == 2, nCls + 1, nCls), reps = 10000, pItem = 0.7,innerLinkage = "ward.D2", finalLinkage = "ward.D2", distance = "euclidean")
  
  
  #识别最佳聚类数-----------------------------
  #PAC确定最佳聚类数K
  pac <- lapply(results[2:length(results)], function(x){
    PAC(x$consensusMatrix,lower = 0.1, upper = 0.9)
  })
  bestK1 <- which.min(as.numeric(pac))+1
  
  #elbow确定最佳聚类数
  triangle <- function (m, mode = 1) 
  {
    n = dim(m)[1]
    nm = matrix(0, ncol = n, nrow = n)
    fm = m
    nm[upper.tri(nm)] = m[upper.tri(m)]
    fm = t(nm) + nm
    diag(fm) = diag(m)
    nm = fm
    nm[upper.tri(nm)] = NA
    diag(nm) = NA
    vm = m[lower.tri(nm)]
    if (mode == 1) {
      return(vm)
    }
    else if (mode == 3) {
      return(fm)
    }
    else if (mode == 2) {
      return(nm)
    }
  }
  
  ml <- lapply(results[2:maxK], function(x){x$ml})
  areaK <- c()
  for (i in 1:length(ml)) {
    v = triangle(ml[[i]], mode = 1)
    h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1/100))
    h$counts = cumsum(h$counts)/sum(h$counts)
    thisArea = 0
    for (bi in 1:(length(h$breaks) - 1)) {
      thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
      bi = bi + 1
    }
    areaK = c(areaK, thisArea)
    if(i == 1){
      deltaK = areaK
    } else{
      deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i - 1])
    }
  }
  for (k in 1:length(deltaK)) {
    if(k == 1){
      ddeltaK <- c()
    }else if(k>=2 && k<length(deltaK)){
      ddeltaK = c(ddeltaK,((deltaK[k-1]-deltaK[k]) - abs((deltaK[k]-deltaK[k+1]))))
    }
  }
  bestK2 <- as.numeric(which.max(ddeltaK)+2)
  
  #calinsky确定最佳聚类数
  calinsky <- function(hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
    msg <- ""
    if (is.null(dist)) {
      require(clue)
      dist <- sqrt(as.cl_ultrametric(hhc))
      # message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
    } else if (attr(dist, "method") != "euclidean") {
      require(clue)
      dist <- sqrt(as.cl_ultrametric(hhc))
      # message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
    }
    dist <- as.matrix(dist)^2
    A <- -dist/2
    A_bar <- apply(A, 1, mean)
    totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
    n <- length(hhc$order)
    ans <- rep(0, gMax)
    for (g in 2:gMax) {
      cclust <- cutree(hhc, k = g)
      withinSum <- 0
      for (k in 1:g) {
        if (sum(cclust == k) == 1) 
          next
        A <- as.matrix(-dist/2)[cclust == k, cclust == k]
        A_bar <- apply(A, 1, mean)
        withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
                                       mean(A))
      }
      betweenSum <- totalSum - withinSum
      betweenSum <- betweenSum/(g - 1)
      withinSum <- withinSum/(n - g)
      ans[g] <- betweenSum/withinSum
    }
    class(ans) <- "calinsky"
    attr(ans, "message") <- msg
    return(ans)
  }
  
  sHc <- hclust(as.dist(jaccard), method = "average")
  aCalinsky <- calinsky(sHc, gMax = 10)
  bestK3 <- which.max(aCalinsky)
  
  #显示由不同的算法得出的最佳聚类数
  if(bestK.method == 'PAC'){
    nCls <- bestK1
    cat('Best K of PAC is',nCls)
  } else if(bestK.method == 'elbow'){
    nCls <- bestK2
    cat('Best K of elbow is',nCls)
  } else if(bestK.method == 'calinsky'){
    nCls <- bestK3
    cat('Best K of calinsky is',nCls)
  }
  
  #根据最佳聚类数提取结果---nCls=4--
  aConsensus <- results[[nCls]]
  
  obj <- aConsensus
  G <- nCls
  
  consHc <- obj$consensusTree
  consClust <- as.factor(cutree(consHc, k = G))
  names(consClust) <- colnames(NES)
  levels(consClust) <- palette()[1:G]
  
  col <- colorRampPalette(c("lemonchiffon2", "#053061"))(51)
  
  Colv  <- as.dendrogram(consHc)
  
  consMatrix <- obj$consensusMatrix
  colnames(consMatrix) <- colnames(NES)
  consMatrix <- consMatrix[, consHc$order]
  
  # RowSideColors <- as.character(consClust[colnames(consMatrix)])
  RowSideColors1 <- c('#436b95','#436b95','#436b95','#436b95','#436b95','#436b95','#436b95','#436b95',
                      '#b81d25','#b81d25','#b81d25','#b81d25','#b81d25','#b81d25','#b81d25',
                      '#2f7e2f', '#2f7e2f', '#2f7e2f', '#2f7e2f', '#2f7e2f', '#2f7e2f', 
                      '#800080','#800080','#800080','#800080','#800080')
  
  library(heatmap3)
  llabRow <- colnames(consMatrix)
  ffile <- paste0(figure.path, "/heatmap_allClusters.pdf")
  pdf(file = ffile, width = 12, height = 8)
  heatmap3(t(consMatrix), col = col, scale = "none",
           labCol = NA, labRow = llabRow, cexRow = 0.5,
           Colv = Colv, Rowv = NA,
           RowSideColors = RowSideColors1)
  dev.off()
  
  ffile <- paste0(data.path, "/aConsensus_allClusters.rda")
  save(aConsensus, consHc, consClust, file = ffile)
}

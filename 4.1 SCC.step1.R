##### Step 1 Function #####


SCC.step1 <- function(NESlist,data.path = 'Rdata/SCC.step1',
                      figure.path = 'Figures/SCC.step1'){
  ## creat folder
  if(!file.exists(data.path)){dir.create(data.path,recursive = TRUE)}
  if(!file.exists(figure.path)){dir.create(figure.path,recursive = TRUE)}
  
  library(foreach)
  library(doParallel)
  #设置并行
  cores <- detectCores()
  cl <- makeCluster(cores - 2) # not to overload your computer
  registerDoParallel(cl)
  
  foreach(ss = 1:length(NESlist),.packages = c("ConsensusClusterPlus")) %dopar% { 
    sn <- names(NESlist)[ss]
    dat <- NESlist[[sn]]
    dat.s <- scale(t(dat[,1:ncol(dat)]))
    sn.figure.path <- paste0(figure.path,'/',sn)
    #进行ConsensusCluster
    input <- t(dat.s)
    maxK=10
    results <- ConsensusClusterPlus(d = input,
                                    maxK=maxK,     #K最大值
                                    reps=1000,    #重抽样次数
                                    pItem=0.9,   #样品抽样比例
                                    pFeature=1,  #Feature抽样比例
                                    title= sn.figure.path,
                                    #innerLinkage="complete",
                                    #finalLinkage="complete",
                                    clusterAlg="pam",  #聚类算法
                                    distance="euclidean",    #计算距离的方法,manhattan,
                                    seed=123456,
                                    plot="png")   #如果不设置，图片结果仅输出到屏幕
    
    
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
    areaK <- deltaK <- ddeltaK <- c()
    axis <- data.frame()
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
      
      if(i <= 2){
        ddeltaK <- c()
      }else if(i>2){
        ddeltaK = c(ddeltaK,(deltaK[i]-deltaK[i-1])-(deltaK[i-1]-deltaK[i-2]))
      }
    }
    
    #保存结果
    bestK <- as.integer(which.max(ddeltaK)+2)
    save(results,bestK, file = paste0(data.path,'/',sn,'_Consensus.rda'))
    cat('The best K of',sn, 'is ',bestK)
  }
  stopCluster(cl)
}



### Differential expression/pathway analysis more than 2 groups
DEAgroups <- function(ddata, groups, method = c("MWW", "t-test")){
  gg <- sort(unique(groups))
  
  ans <- vector("list", length(gg))
  names(ans) <- gg
  for(i in 1:length(gg)){
    whichOfInterest <- names(groups)[groups == gg[i]]
    theOthers <- setdiff(colnames(ddata), whichOfInterest)
    diffActivity <- apply(ddata, 1, function(x){
      if(method == "MWW") suppressWarnings(a <- wilcox.test(x[whichOfInterest], x[theOthers]))
      if(method == "t-test") suppressWarnings(a <- t.test(x[whichOfInterest], x[theOthers]))
      a <- c(as.numeric(a$statistic), a$p.value)
      return(a)
    })
    diffActivity <- t(diffActivity)
    colnames(diffActivity) <- c("wTest", "pValue")
    fc <- rowMeans(ddata[, whichOfInterest]) - rowMeans(ddata[, theOthers])
    qValue <- p.adjust(diffActivity[, "pValue"], method = "fdr")
    
    diffActivity <- data.frame(statistic = diffActivity[, "wTest"], dm = fc, p.value = diffActivity[, "pValue"], fdr = qValue)
    diffActivity <- diffActivity[, -1]
    colnames(diffActivity) <- c("logFC", "pValue", "qValue")
    # diffActivity <- diffActivity[order(diffActivity$p.value), ]
    # diffActivity <- diffActivity[order(diffActivity$dm, decreasing = T), ]
    ans[[i]] <- diffActivity
  }
  return(ans)
}

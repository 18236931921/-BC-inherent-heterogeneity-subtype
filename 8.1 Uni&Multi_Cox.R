rm(list = ls())
library(data.table)
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)

#多重插补----
load('Results/6. bulk_NTP&survival/TCGA_cluster.rda')

  #插补前数据处理
  up_df1 <- column_to_rownames(dd,"ID")
  colnames(up_df1)[1] <- "Cluster"
  up_df1$DFS <- as.numeric(up_df1$DFS)
  up_df1$PFS <- as.numeric(up_df1$PFS)
  up_df1$DSS <- as.numeric(up_df1$DSS)
  up_df1$Age <- factor(up_df1$Age)
  up_df1$Gender <- factor(up_df1$Gender)
  up_df1$Stage <- factor(up_df1$Stage,levels = c("I","II","III","IV"),
                         labels = c("Stage I","Stage II","Stage III","Stage IV"))
  up_df1$Grade <- factor(up_df1$Grade,levels = c("High","Low"),
                         labels = c("Grade High","Grade Low"))
  up_df1$T <- factor(up_df1$T)
  up_df1$N <- factor(up_df1$N)
  up_df1$M <- factor(up_df1$M)
  up_df1$Subtype <- factor(up_df1$Subtype)
  up_df1 <- up_df1[,-c(4,5,8,9,18)]
  #缺失值插补
  str(up_df1)
  library(mice)
  impdatas<-mice(up_df1,seed=20202020)
  fit<-with(impdatas,glm(Cluster~.,family=binomial(),data = up_df1))
  impdatas$data
  pooled<-pool(fit)
  completedf <-complete(impdatas,2)
  table(is.na(completedf))
  str(completedf)
  
  #插补后进一步处理临床信息
  cc <- completedf
  cc$Stage <- ifelse(cc$Stage==c("Stage I","Stage II"),"Stage I|II","Stage III|IV")
  cc$Stage <- factor(cc$Stage,levels = c("Stage I|II","Stage III|IV"))
  
  cc$T <- ifelse(cc$T==c("T1","T2"),"T 1|2","T 3|4")
  cc$T <- factor(cc$T)
  
  cc$N <- ifelse(cc$N==c("N1","N2"),"N 1|2","N 3|4")
  cc$N <- factor(cc$N)
  
  str(cc)
  
  cc <- cc[,-1]
  
  save(cc,file ="Results/8. COX&clinical analysis/complete clinical data.RData")


#单因素及多因素------
rm(list = ls())
load("Results/8. COX&clinical analysis/complete clinical data.RData")
precox <- cc %>% select("OS","OS.time","PFS","PFS.time",everything())
#OS单因素

if(T){
  df <- data.frame()
  for (i in colnames(precox)[5:ncol(precox)]) {
    fit <- coxph(Surv(OS.time,OS)~get(i),data = precox)
    fit1 <- summary(fit)
    coef <- coef(fit);HR <- exp(coef(fit))
    HR_l <- exp(confint(fit))[1];HR_H <- exp(confint(fit))[2];
    p_val <- fit1$coefficients[5]
    df <- rbind(df,
                cbind(i,coef,HR,HR_l,HR_H,p_val))
  }
  name<- df$i
  unicox_OS <- df[,-1] %>% apply(2,as.numeric) %>% as.data.frame()
  rownames(unicox_OS) <- name
}

#OS多因素
i <- rownames(unicox_OS)[unicox_OS$p_val<0.05]
precox_multi <- precox[,c("OS.time","OS", "Age" ,"Stage","T","N",'M','Subtype')]
fit <- coxph(Surv(OS.time,OS)~.,data = precox_multi)
fit1 <- summary(fit)
multi_OS <- data.frame(coef = coef(fit),HR = exp(coef(fit)),
                       HR_l = exp(confint(fit))[,1],
                       HR_H = exp(confint(fit))[,2],
                       p_val = fit1$coefficients[,5])


#PFS单因素
if(T){
  df <- data.frame()
  for (i in colnames(precox)[5:ncol(precox)]) {
    fit <- coxph(Surv(PFS.time,PFS)~get(i),data = precox)
    fit1 <- summary(fit)
    coef <- coef(fit);HR <- exp(coef(fit))
    HR_l <- exp(confint(fit))[1];HR_H <- exp(confint(fit))[2];
    p_val <- fit1$coefficients[5]
    df <- rbind(df,
                cbind(i,coef,HR,HR_l,HR_H,p_val))
  }
  name<- df$i
  unicox_PFS <- df[,-1] %>% apply(2,as.numeric) %>% as.data.frame()
  rownames(unicox_PFS) <- name
}

#DFS多因素
i <- rownames(unicox_DFS)[unicox_DFS$p_val<0.05]
precox_multi <- precox[,c("PFS.time","PFS","Stage","T","N",'M','Subtype')]
fit <- coxph(Surv(PFS.time,PFS)~.,data = precox_multi)
fit1 <- summary(fit)
multi_PFS <- data.frame(coef = coef(fit),HR = exp(coef(fit)),
                        HR_l = exp(confint(fit))[,1],
                        HR_H = exp(confint(fit))[,2],
                        p_val = fit1$coefficients[,5])

save(unicox_PFS,unicox_OS,multi_PFS,multi_OS,file = "Results/8. COX&clinical analysis/BLCA cox result.RData")



##森林图---------------------------------------------
rm(list = ls())
load("Results/8. COX&clinical analysis/BLCA cox result.RData")

#森林图unicox_OS----
df <- unicox_OS[unicox_OS$p_val<0.05,]
df <- rownames_to_column(df,var = "id")
df <- df[,-2]
unicox <- df
colnames(unicox) <- c('id','HR','HRL','HRH','P')
unicox$P <- abs(unicox$P)
unicox[is.na(unicox)] <- 1
unicox[unicox == Inf] <- max(unicox$HRL)
# unicox$ll <- ifelse(unicox$P<0.0001,'****',ifelse(unicox$P<0.001,'***',ifelse(unicox$P<0.01,'**',ifelse(unicox$P<0.05,'*',''))))
cols <- c('#f2cc8f','#81b29a','#118ab2','#e07a5f')
cols2 <- c('#eef3fe','#fff4f5','#f9f5ec','#f0fefb')
# unicox$Cell <- factor(unicox$Cell,levels = rev(unique(unicox$Cell)))

# hrtable <- unicox[unicox$Cluster == 'CMS1',-c(2:3)]
hrtable <- unicox
hrtable$ss <- ifelse(hrtable$P < 0.05, '*',
                     ifelse(hrtable$P < 0.01, '**',
                            ifelse(hrtable$P < 0.001, '***','')))

tabletext <- cbind(c("id",hrtable$id),
                   c("HR (95% CI)",format(paste0(round(as.numeric(hrtable$HR),2),' [',round(as.numeric(hrtable$HRL),2),'-',round(as.numeric(hrtable$HRH),2),']'))),
                   # c("HRL",format(round(as.numeric(hrtable$HRL),3),nsmall = 3)),
                   # c("HRH",format(round(as.numeric(hrtable$HRH),3),nsmall = 3)),
                   c("pvalue",formatC(as.numeric(hrtable$P), format = "e", digits = 2)),
                   c("",format(hrtable$ss)))

pdf(file = 'Results/8. COX&clinical analysis/TCGA_uni_OS_clin_forestplot.pdf',height = 8, width = 10)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
           lower=c(NA,as.numeric(hrtable$HRL)), #log2(95%置信区间下限)
           upper=c(NA,as.numeric(hrtable$HRH)),#log2(95%置信区间上限)
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawNormalCI",#box类型选择钻石
           col=fpColors(box="#e15151", lines="#158bb8", zero = "black"),#box颜色
           boxsize=0.1,#box大小固定
           lwd.ci=2,
           ci.vertices.height = 0.1,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           clip = c(0.5,3),#横坐标刻度根据需要可随意设置
           xticks = c(0,0.5,1,1.5,2),
           lwd.xaxis=2, vertices = TRUE,
           # xlab=expression("HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2)),#第二行顶部加灰色虚线
           txt_gp=fpTxtGp(label=gpar(cex=1.5),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.5),
                          title=gpar(cex=2)),
           lineheight = unit(.85,"cm"),#固定行高
           colgap = unit(0.4,"cm"),
           mar=unit(rep(0.5, times = 4), "cm"),
           new_page = F,
           xlog = F,
           title = "TCGA_OS"
)

dev.off()



#森林图unicox_PFS----
df <- unicox_PFS[unicox_PFS$p_val<0.05,]
df <- rownames_to_column(df,var = "id")
df <- df[,-2]
unicox <- df
colnames(unicox) <- c('id','HR','HRL','HRH','P')
unicox$P <- abs(unicox$P)
unicox[is.na(unicox)] <- 1
unicox[unicox == Inf] <- max(unicox$HRL)
# unicox$ll <- ifelse(unicox$P<0.0001,'****',ifelse(unicox$P<0.001,'***',ifelse(unicox$P<0.01,'**',ifelse(unicox$P<0.05,'*',''))))
cols <- c('#f2cc8f','#81b29a','#118ab2','#e07a5f')
cols2 <- c('#eef3fe','#fff4f5','#f9f5ec','#f0fefb')
# unicox$Cell <- factor(unicox$Cell,levels = rev(unique(unicox$Cell)))

# hrtable <- unicox[unicox$Cluster == 'CMS1',-c(2:3)]
hrtable <- unicox
hrtable$ss <- ifelse(hrtable$P < 0.05, '*',
                     ifelse(hrtable$P < 0.01, '**',
                            ifelse(hrtable$P < 0.001, '***','')))

tabletext <- cbind(c("id",hrtable$id),
                   c("HR (95% CI)",format(paste0(round(as.numeric(hrtable$HR),2),' [',round(as.numeric(hrtable$HRL),2),'-',round(as.numeric(hrtable$HRH),2),']'))),
                   # c("HRL",format(round(as.numeric(hrtable$HRL),3),nsmall = 3)),
                   # c("HRH",format(round(as.numeric(hrtable$HRH),3),nsmall = 3)),
                   c("pvalue",formatC(as.numeric(hrtable$P), format = "e", digits = 2)),
                   c("",format(hrtable$ss)))

pdf(file = 'Results/8. COX&clinical analysis/TCGA_uni_PFS_clin_forestplot.pdf',height = 8, width = 10)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
           lower=c(NA,as.numeric(hrtable$HRL)), #log2(95%置信区间下限)
           upper=c(NA,as.numeric(hrtable$HRH)),#log2(95%置信区间上限)
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawNormalCI",#box类型选择钻石
           col=fpColors(box="#e15151", lines="#158bb8", zero = "black"),#box颜色
           boxsize=0.1,#box大小固定
           lwd.ci=2,
           ci.vertices.height = 0.1,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           clip = c(0.5,3),#横坐标刻度根据需要可随意设置
           xticks = c(0,0.5,1,1.5,2),
           lwd.xaxis=2, vertices = TRUE,
           # xlab=expression("HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2)),#第二行顶部加灰色虚线
           txt_gp=fpTxtGp(label=gpar(cex=1.5),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.5),
                          title=gpar(cex=2)),
           lineheight = unit(.85,"cm"),#固定行高
           colgap = unit(0.4,"cm"),
           mar=unit(rep(0.5, times = 4), "cm"),
           new_page = F,
           xlog = F,
           title = "TCGA_PFS"
)

dev.off()



#绘制森林图multi_OS----
df <- multi_OS[multi_OS$p_val<0.05,]
df <- rownames_to_column(df,var = "id")
df <- df[,-2]
df$id <- c('Age≤65','N 3|4','M1')
unicox <- df
colnames(unicox) <- c('id','HR','HRL','HRH','P')
unicox$P <- abs(unicox$P)
unicox[is.na(unicox)] <- 1
unicox[unicox == Inf] <- max(unicox$HRL)
# unicox$ll <- ifelse(unicox$P<0.0001,'****',ifelse(unicox$P<0.001,'***',ifelse(unicox$P<0.01,'**',ifelse(unicox$P<0.05,'*',''))))
cols <- c('#f2cc8f','#81b29a','#118ab2','#e07a5f')
cols2 <- c('#eef3fe','#fff4f5','#f9f5ec','#f0fefb')
# unicox$Cell <- factor(unicox$Cell,levels = rev(unique(unicox$Cell)))

# hrtable <- unicox[unicox$Cluster == 'CMS1',-c(2:3)]
hrtable <- unicox
hrtable$ss <- ifelse(hrtable$P < 0.05, '*',
                     ifelse(hrtable$P < 0.01, '**',
                            ifelse(hrtable$P < 0.001, '***','')))

tabletext <- cbind(c("id",hrtable$id),
                   c("HR (95% CI)",format(paste0(round(as.numeric(hrtable$HR),2),' [',round(as.numeric(hrtable$HRL),2),'-',round(as.numeric(hrtable$HRH),2),']'))),
                   # c("HRL",format(round(as.numeric(hrtable$HRL),3),nsmall = 3)),
                   # c("HRH",format(round(as.numeric(hrtable$HRH),3),nsmall = 3)),
                   c("pvalue",formatC(as.numeric(hrtable$P), format = "e", digits = 2)),
                   c("",format(hrtable$ss)))

pdf(file = 'Results/8. COX&clinical analysis/TCGA_multi_OS_clin_forestplot.pdf',height = 8, width = 10)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
           lower=c(NA,as.numeric(hrtable$HRL)), #log2(95%置信区间下限)
           upper=c(NA,as.numeric(hrtable$HRH)),#log2(95%置信区间上限)
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawNormalCI",#box类型选择钻石
           col=fpColors(box="#e15151", lines="#158bb8", zero = "black"),#box颜色
           boxsize=0.1,#box大小固定
           lwd.ci=2,
           ci.vertices.height = 0.1,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           clip = c(0.5,3),#横坐标刻度根据需要可随意设置
           xticks = c(0,0.5,1,1.5,2),
           lwd.xaxis=2, vertices = TRUE,
           # xlab=expression("HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2)),#第二行顶部加灰色虚线
                           # "4" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.5),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.5),
                          title=gpar(cex=2)),
           lineheight = unit(.85,"cm"),#固定行高
           colgap = unit(0.4,"cm"),
           mar=unit(rep(0.5, times = 4), "cm"),
           new_page = F,
           xlog = F,
           title = "TCGA_multi_OS"
)

dev.off()


#绘制森林图multi_PFS----
df <- multi_PFS[multi_PFS$p_val<0.05,]
df <- rownames_to_column(df,var = "id")
df <- df[,-2]
df$id <- c('N 3|4','M1')
unicox <- df
colnames(unicox) <- c('id','HR','HRL','HRH','P')
unicox$P <- abs(unicox$P)
unicox[is.na(unicox)] <- 1
unicox[unicox == Inf] <- max(unicox$HRL)
# unicox$ll <- ifelse(unicox$P<0.0001,'****',ifelse(unicox$P<0.001,'***',ifelse(unicox$P<0.01,'**',ifelse(unicox$P<0.05,'*',''))))
cols <- c('#f2cc8f','#81b29a','#118ab2','#e07a5f')
cols2 <- c('#eef3fe','#fff4f5','#f9f5ec','#f0fefb')
# unicox$Cell <- factor(unicox$Cell,levels = rev(unique(unicox$Cell)))

# hrtable <- unicox[unicox$Cluster == 'CMS1',-c(2:3)]
hrtable <- unicox
hrtable$ss <- ifelse(hrtable$P < 0.05, '*',
                     ifelse(hrtable$P < 0.01, '**',
                            ifelse(hrtable$P < 0.001, '***','')))

tabletext <- cbind(c("id",hrtable$id),
                   c("HR (95% CI)",format(paste0(round(as.numeric(hrtable$HR),2),' [',round(as.numeric(hrtable$HRL),2),'-',round(as.numeric(hrtable$HRH),2),']'))),
                   # c("HRL",format(round(as.numeric(hrtable$HRL),3),nsmall = 3)),
                   # c("HRH",format(round(as.numeric(hrtable$HRH),3),nsmall = 3)),
                   c("pvalue",formatC(as.numeric(hrtable$P), format = "e", digits = 2)),
                   c("",format(hrtable$ss)))

pdf(file = 'Results/8. COX&clinical analysis/TCGA_multi_PFS_clin_forestplot.pdf',height = 8, width = 10)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
           lower=c(NA,as.numeric(hrtable$HRL)), #log2(95%置信区间下限)
           upper=c(NA,as.numeric(hrtable$HRH)),#log2(95%置信区间上限)
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawNormalCI",#box类型选择钻石
           col=fpColors(box="#e15151", lines="#158bb8", zero = "black"),#box颜色
           boxsize=0.1,#box大小固定
           lwd.ci=2,
           ci.vertices.height = 0.1,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           clip = c(0.5,3),#横坐标刻度根据需要可随意设置
           xticks = c(0,0.5,1,1.5,2),
           lwd.xaxis=2, vertices = TRUE,
           # xlab=expression("HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2)),#第二行顶部加灰色虚线
                           # "4" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.5),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.5),
                          title=gpar(cex=2)),
           lineheight = unit(.85,"cm"),#固定行高
           colgap = unit(0.4,"cm"),
           mar=unit(rep(0.5, times = 4), "cm"),
           new_page = F,
           xlog = F,
           title = "TCGA_multi_PFS"
)

dev.off()


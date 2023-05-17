# 1. Cibersort
rm(list=ls())
symbol=data.table::fread(file="GSE157153_endometriosis_associated.txt")
colnames(symbol)[1]="Gene symbol"

exprSet=column_to_rownames(symbol,var = "Gene symbol")
write.table(exprSet,file="symbol.txt",sep="\t",quote=F,col.names=T,row.names = T)   
source("CIBERSORT.R")
results=CIBERSORT("LM22.txt", "symbol.txt", perm=1000, QN=FALSE)

write.table(results,file="cibersort_results.txt",sep="\t",quote=F,col.names=T,row.names = T)
#########################################################################################################
results1=read.table(file="cibersort_results.txt",header=T,sep="\t",check.names=F,row.names=1)
results_filter=results1[,c(-23:-25)]
###################################################绘制小提琴图#############################################
table(results_filter$group)
normal=37                                                   
tumor= 29   
rt=results_filter

library(vioplot)

pdf("vioplot.pdf",height=8,width=16)             

par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,66),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=15,
     col="white",
     xaxt="n")
legend("topright",                                        #图例位置为上方
       legend=c("EAOC","Endometriosis"),
       #ncol=4,
       # cex=0.8,
       #bty="n",
       col=c('#bc3c29',"#0072b5"),
       lty=1,lwd=2)

for(i in 1:ncol(rt)){
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = "#0072b5")##blue
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = '#bc3c29')###red
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,66,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2,family="serif")
}
dev.off()
###########################################免疫细胞相关性分析########################
library(corrgram)
library(corrplot)
library(corrplot)
library(ggplot2)
library(ggcorrplot)


results1=read.table(file="cibersort_results.txt",header=T,sep="\t",check.names=F,row.names=1)
results_filter=results1[,c(-23:-25)]
M <- cor(results_filter,method = "pearson")
head(M)
# 自定义渐变色板 colorRampPalette()
# 用法：
# 1.输入调色板的主要颜色，返回函数
# 2.在函数中输入想要获取的颜色数量，返回颜色数值
library(scales) #调用scales库中的show_col函数
mypalette<- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
mycolors<- mypalette(200)#200种颜色，返回颜色数值
show_col(mycolors)
corrplot(M, method = "square",##调整为正方形
         col=mycolors, #200种颜色，返回颜色数值
         tl.col ="black",# # 指定文本标签的颜色
         tl.srt = 45, #角度
         order = "AOE",##
         type="full",#保留全部  lower  upper full 
         tl.cex=0.8#字体大小
         )      
#addCoef.col = "black", #添加相关系数
#diag=FALSE #去除自身相关
##################################单基因与免疫细胞相关性分析######################################
# 样品1:expr_data
# 样品2:immu_data
# 打开两个想进行相关性分析的清洁数据
# 数据1为行名为样品名，列名为基因名
# 数据2为行名为样品名，列名为想做的和基因有相关性关系的，如免疫细胞

results1=read.table(file="cibersort_results.txt",header=T,sep="\t",check.names=F,row.names=1)
results_filter=results1[,c(-23:-25)]

immu_data=results_filter

symbol=read.table(file="symbol.txt",header=T,sep="\t",check.names=F,row.names=1)
a=as.data.frame(t(symbol))
#FAM83H SPINT1

expr_data=a

#绘制棒棒糖图
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(dplyr)
#进行相关性计算
gene <- "EDNRA"
y <- as.numeric(expr_data[,gene])

cor_data <- do.call(rbind,lapply(colnames(immu_data),function(x){
  dd <- cor.test(as.numeric(immu_data[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


# 画单基因与样本2中所有列名的相关性的棒棒糖图

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

dat <- cor_data
gene <- 'EDNRA'   #修改基因名
#根据pvalue定义圆圈的颜色
points.color = fcolor(x=dat$p.value,p.col=p.col)
dat$points.color = points.color

#根据相关系数定义圆圈的大小
points.cex = fcex(x=dat$cor)-0.8
dat$points.cex = points.cex
dat=dat[order(dat$cor),]

########绘制图形########
xlim = ceiling(max(abs(dat$cor))*10)/10         #x轴范围
pdf(file="EDNRA.pdf", width=9, height=7)      #输出图形  ##修改文件名
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(dat)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(dat),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=dat$cor,y0=1:nrow(dat),x1=0,y1=1:nrow(dat),lwd=4)
#绘制图形的圆圈
points(x=dat$cor,y = 1:nrow(dat),col = dat$points.color,pch=16,cex=dat$points.cex)+
  scale_size_continuous(range =c(2,4))
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(dat),dat$cell,adj=1,xpd=T,cex=1.5)
#展示pvalue
pvalue.text=ifelse(dat$p.value<0.001,'<0.001',sprintf("%.03f",dat$p.value))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(dat),pvalue.text,adj=0,xpd=T,col=ifelse(abs(dat$cor)>redcutoff_cor & dat$p.value<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()



########################################
results1=read.table(file="cibersort_results.txt",header=T,sep="\t",check.names=F,row.names=1)
results_filter=results1[,c(-23:-25)]

immu_data=results_filter

symbol=read.table(file="symbol.txt",header=T,sep="\t",check.names=F,row.names=1)
a=as.data.frame(t(symbol))
#FAM83H SPINT1

expr_data=a

#绘制棒棒糖图
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(dplyr)
#进行相关性计算
gene <- "OCLN"
y <- as.numeric(expr_data[,gene])

cor_data <- do.call(rbind,lapply(colnames(immu_data),function(x){
  dd <- cor.test(as.numeric(immu_data[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


# 画单基因与样本2中所有列名的相关性的棒棒糖图

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

dat <- cor_data
gene <- 'OCLN'   #修改基因名
#根据pvalue定义圆圈的颜色
points.color = fcolor(x=dat$p.value,p.col=p.col)
dat$points.color = points.color

#根据相关系数定义圆圈的大小
points.cex = fcex(x=dat$cor)-0.8
dat$points.cex = points.cex
dat=dat[order(dat$cor),]

########绘制图形########
xlim = ceiling(max(abs(dat$cor))*10)/10         #x轴范围
pdf(file="OCLN.pdf", width=9, height=7)      #输出图形  ##修改文件名
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(dat)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(dat),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=dat$cor,y0=1:nrow(dat),x1=0,y1=1:nrow(dat),lwd=4)
#绘制图形的圆圈
points(x=dat$cor,y = 1:nrow(dat),col = dat$points.color,pch=16,cex=dat$points.cex)+
  scale_size_continuous(range =c(2,4))
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(dat),dat$cell,adj=1,xpd=T,cex=1.5)
#展示pvalue
pvalue.text=ifelse(dat$p.value<0.001,'<0.001',sprintf("%.03f",dat$p.value))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(dat),pvalue.text,adj=0,xpd=T,col=ifelse(abs(dat$cor)>redcutoff_cor & dat$p.value<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()

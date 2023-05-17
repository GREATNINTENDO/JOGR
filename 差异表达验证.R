rm(list=ls())
#####157153
symbol=read.table(file="symbol.txt",header=T,sep="\t",check.names=F,row.names=1)
symbol=as.data.frame(t(symbol))
gene=c("EDNRA","OCLN")

exprset=symbol[,gene]

group=data.table::fread(file="group.txt")
exprset$group=group$type

#设置比较组
data=exprset

group=levels(factor(data$group))
data$group=factor(data$group, levels=c("Tumor", "Control"))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggpubr)
library(ggsci)
library(ggplot2)


library(ggpubr)
ggviolin(data, x="group", y="OCLN", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of OCLN",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN_157153.PDF", width=3.5,height=5)
print(boxplot)
dev.off()

ggviolin(data, x="group", y="EDNRA", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of EDNRA",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA_157153.PDF", width=3.5,height=5)
print(boxplot)
dev.off()





#绘制箱线图
boxplot=ggboxplot(data, x="group", y="EDNRA", fill="group",
                  xlab="",
                  ylab="Expression of EDNRA",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

boxplot=ggboxplot(data, x="group", y="OCLN", fill="group",
                  xlab="",
                  ylab="Expression of OCLN",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN.pdf",width=3.5,height=5)
print(boxplot)
dev.off()


###################################
#GSE18520
rm(list=ls())
symbol=read.table(file="exprSet_GSE18520.txt",header=T,sep="\t",check.names=F,row.names=1)
symbol=as.data.frame(t(symbol))
gene=c( "OCLN", "EDNRA" )

exprset=symbol[,gene]

group=factor(c(rep('Tumor',53),rep('Control',10)))

exprset$group=group

#exprset$group=as.factor(exprset$group)
data=exprset
group=levels(factor(data$group))
data$group=factor(data$group, levels=c("Tumor", "Control"))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggpubr)
library(ggsci)
library(ggplot2)

ggviolin(data, x="group", y="OCLN", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of OCLN",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN_18520.PDF", width=3.5,height=5)
print(boxplot)
dev.off()

ggviolin(data, x="group", y="EDNRA", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of EDNRA",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA_18520.PDF", width=3.5,height=5)
print(boxplot)
dev.off()






#绘制箱线图
boxplot=ggboxplot(data, x="group", y="EDNRA", fill="group",
                  xlab="",
                  ylab="Expression of EDNRA",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA_18520.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

boxplot=ggboxplot(data, x="group", y="OCLN", fill="group",
                  xlab="",
                  ylab="Expression of OCLN",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN_18520.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

#######################################
rm(list=ls())
symbol=read.table(file="exprSet_GSE26712.txt",header=T,sep="\t",check.names=F,row.names=1)
symbol=as.data.frame(t(symbol))
gene=c( "OCLN", "EDNRA")

exprset=symbol[,gene]

group=factor(c(rep('Control',10),rep('Cancer',185)))#查看标本分组信息，构建groupList分组
exprset$group=group

#exprset$group=as.factor(exprset$group)
data=exprset
group=levels(factor(data$group))
data$group=factor(data$group, levels=c("Cancer", "Control"))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggpubr)
library(ggsci)
library(ggplot2)
#绘制箱线图
boxplot=ggboxplot(data, x="group", y="EDNRA", fill="group",
                  xlab="",
                  ylab="Expression of EDNRA",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA_26712.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

boxplot=ggboxplot(data, x="group", y="OCLN", fill="group",
                  xlab="",
                  ylab="Expression of OCLN",
                  legend.title="",
                  #add = "jitter",
                  palette = c("#bc3c29","#0072b5") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN_26712.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

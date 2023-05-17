library(clusterProfiler)
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library('GSEABase')
library(fgsea)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
setwd("D:/Endometriosis_Cancer/10.GSEA")
#devtools::install_github("nicolash2/gggsea")
symbol=read.table(file="symbol.txt",header=T,sep="\t",check.names=F,row.names=1)
symbol=as.data.frame(t(symbol))

dat=read.table(file="symbol.txt",header=T,sep="\t",check.names=F,row.names=1)

risk=as.vector(ifelse(symbol$OCLN>median(symbol$OCLN), "High", "Low"))
risk=as.vector(ifelse(symbol$EDNRA>median(symbol$EDNRA), "High", "Low"))
Entire=cbind(risk, symbol)

groupList=as.factor(Entire$risk)#查看标本分组信息，构建groupList分组
table(groupList)

################################LIMMA差异分析##########
#2.构建实验设计矩阵design
library(limma)
design <- model.matrix(~0+factor(groupList))
colnames(design) <- levels(factor(groupList))#design的列名
rownames(design)=colnames(dat)                 #design的行名
#3.构建对比模型，比较两个实验条件下表达数据
contrast.matrix<-makeContrasts(High-Low,levels = design)
contrast.matrix##这个矩阵声明，我们要把cancer组跟control组进行差异分析比较

#4.差异分析
##step1 线性模型拟合
fit <- lmFit(dat,design)
##step2 根据对比模型进行差值计算
fit2 <- contrasts.fit(fit, contrast.matrix)
##step3 贝叶斯检验
fit2 <- eBayes(fit2)
##step4 生成基因的检验结果报告
allDiff=topTable(fit2,adjust='fdr',number=200000)
#write.table(allDiff,file="limmaTab.txt",sep="\t",quote=F)

##step5 按log(foldchange)排序,为RRA做准备
allLimma=allDiff
allLimma=allLimma[order(allLimma$logFC),]
head(allLimma)
library(tidyverse)
out=allLimma
geneList2=out$logFC
names(geneList2)=rownames(out)#提取相关基因及其logFC
geneList2 = sort(geneList2, decreasing = TRUE)#按logFC进行排序
gene_diff=as.data.frame(geneList2)
gene_diff=rownames_to_column(gene_diff,var="symbol")
colnames(gene_diff)[2]="logFC"
write.table(gene_diff, 'gsea_input.txt', sep = '\t', quote = FALSE,col.names = NA)

kegmt<-read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
#kegmt<-read.gmt("h.all.v7.4.symbols.gmt")

library(clusterProfiler)
library(org.Hs.eg.db)
gsea.re1<- clusterProfiler::GSEA(geneList2, pvalueCutoff = 1,TERM2GENE =  kegmt, pAdjustMethod = 'BH') #GSEA分析 pvalueCutoff = 1输出全部
result=gsea.re1 @result
head(result)#result即为GSEA富集的结果

#输出
write.table(gsea.re1, 'gsea_kegg.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#提取显著富集的基因集
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

## 这里去掉了基因集前缀
kegmt$term <- gsub('KEGG_','',kegmt$term)
#[[代表提取元素，比如数据框中data[1]提取第一列，数据类型还是data.frame，data[[1]]提取第一列,数据类型是character
#下面表示取第二列，并表示为character
kegmt.list <- kegmt %>% split(.$term) %>% lapply( "[[", 2)

gsea.re2 <- fgsea(pathways =kegmt.list,#基因集列表
                  stats = geneList2,#排序后的基因level,这里是logFC
                  nperm=1000,#置换检验的次数
                  minSize=1,#富集到某个条目的最小包含基因数，如果基因数小于该值则这个条目将被过滤
                  maxSize=10000)#富集到某个条目的最大包含基因数，如果基因数大于该值则这个条目将被过滤
                  #                 nproc = 0#如果不等于零，则将 BPPARAM 设置为使用的nproc （默认值 = 0）。
                  #                 gseaParam = 1,#GSEA 权重参数（0 为未加权，建议值为 1）
                  #                 BPPARAM  = NULL,#Bplapply中使用的并行化参数。可用于指定要运行的群集。如果没有显式初始化或通过设置“nproc”默认值“，bpparam()”被使用


#提取显著富集的基因集
colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$padj<0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]
#输出
save(gsea.re1,g1,gsea.re2,g2,file = 'gsea.RData')


#devtools::install_github("junjunlab/GseaVis")  #安装GseaVis代码
#install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock')) #Permission denied 代码
#devtools::install_github("junjunlab/jjAnno")
library(GseaVis)
library(jjAnno)
num1=10
#指定名称：图纸
p1<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1])
# 只保留曲线
p2<-gseaNb(object =  gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           subPlot = 1)
#曲线保留和热图retain curve and heatmap
p3<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           subPlot = 2)
# 名称太长截断换行
p4<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           subPlot = 3,
           termWidth = 30,
           curveCol = jjAnno::useMyCol('paired',10))
#标记里的一些名称：标记里
# add gene in specific pathway
mygene <- c("OCLN","EDNRA")

# plot
p5<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           subPlot = 2,
           addGene = mygene)

#基因改变标签颜色和箭头类型：
p6<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           subPlot = 2,
           addGene = mygene,
           arrowType = 'open',
           geneCol = 'black')

#剩下所有图形：
p7<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[c(7,8,9,12,16,19,25)],
           subPlot = 3,
           #addGene = mygene,
           #termWidth = 30,
           curveCol = jjAnno::useMyCol('paired',10),
           legend.position = c(0.9,0.82))
           #rmSegment = TRUE)

# clsaasic with pvalue
p8<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1],
           #addGene = mygene,
           addPval = T,
           pvalX = 0.75,pvalY = 0.8,
           pCol = 'black',
           pHjust = 0)



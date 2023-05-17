setwd("D:/Endometriosis_Cancer/06.Cancer_Endo_Validation")
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型

exprSet=data.table::fread(file = "GSE157153_endometriosis_associated.txt")
exprSet=column_to_rownames(exprSet,var = "gene")
dat=exprSet
group=read.table(file="group.txt",header = T,sep = "\t",check.names = F)
#group=column_to_rownames(group,var = "sample")

groupList=as.factor(group$type)
table(groupList)

#2.构建实验设计矩阵design
library(limma)
design <- model.matrix(~0+factor(groupList))
colnames(design) <- levels(factor(groupList))#design的列名
rownames(design)=colnames(dat)                 #design的行名
#3.构建对比模型，比较两个实验条件下表达数据
contrast.matrix<-makeContrasts(Cancer-Control,levels = design)
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
data=allLimma%>%filter(abs(logFC)>1)%>%filter(adj.P.Val<0.05)

write.table(allLimma,file="GSE157153.txt",sep="\t",quote = F,col.names = T)
write.table(data,file="GSE157153_DEGs.txt",sep="\t",quote = F,col.names = T)
write.table(rownames(data),file="GSE157153_DEGs_GENE.txt",sep="\t",quote = F,col.names = T)
allLimma=read.table(file="GSE157153.txt",header = T,sep = "\t",check.names = F)

#引用包
# install.packages("VennDiagram")
library(VennDiagram)               #引用包
outFile="intersectGenes.txt"       #输出交集基因文件
outPic="venn.pdf"                  #输出图片文件

files=dir()                        #获取目录下所有文件
files=grep("txt$",files,value=T)   #提取TXT结尾的文件
geneList=list()


#读取所有txt文件中的基因信息，保存到GENELIST              #默认读取第一列  为基因名
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)        #读取
  geneNames=as.vector(rt[,1])              #提取基因名
  geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
  uniqGene=unique(geneNames)               #基因取unique
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

head(geneList)


#绘制venn??
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file=outPic, width=5, height=5)
grid.draw(venn.plot)
dev.off()

a=read.table(file="Yellow_module.txt",header = F,sep = "\t",check.names = F)
b=read.table(file="GSE157153.txt",header = F,sep = "\t",check.names = F)
common=intersect(a$V1,b$V1)

write.table(common,file = "Gene_For_MachineLearning.txt",sep="\t",quote = F,col.names = T,row.names = F)



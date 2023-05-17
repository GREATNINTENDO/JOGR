setwd("D:/Endometriosis_Cancer/06.Cancer_Endo_Validation")
rm(list = ls())  ## ħ�ò�����һ�����~
options(stringsAsFactors = F)#�ڵ���as.data.frame��ʱ����stringsAsFactors����ΪFALSE���Ա���character�����Զ�ת��Ϊfactor����

exprSet=data.table::fread(file = "GSE157153_endometriosis_associated.txt")
exprSet=column_to_rownames(exprSet,var = "gene")
dat=exprSet
group=read.table(file="group.txt",header = T,sep = "\t",check.names = F)
#group=column_to_rownames(group,var = "sample")

groupList=as.factor(group$type)
table(groupList)

#2.����ʵ����ƾ���design
library(limma)
design <- model.matrix(~0+factor(groupList))
colnames(design) <- levels(factor(groupList))#design������
rownames(design)=colnames(dat)                 #design������
#3.�����Ա�ģ�ͣ��Ƚ�����ʵ�������±�������
contrast.matrix<-makeContrasts(Cancer-Control,levels = design)
contrast.matrix##�����������������Ҫ��cancer���control����в�������Ƚ�

#4.�������
##step1 ����ģ�����
fit <- lmFit(dat,design)
##step2 ���ݶԱ�ģ�ͽ��в�ֵ����
fit2 <- contrasts.fit(fit, contrast.matrix)
##step3 ��Ҷ˹����
fit2 <- eBayes(fit2)
##step4 ���ɻ���ļ���������
allDiff=topTable(fit2,adjust='fdr',number=200000)
#write.table(allDiff,file="limmaTab.txt",sep="\t",quote=F)

##step5 ��log(foldchange)����,ΪRRA��׼��
allLimma=allDiff
allLimma=allLimma[order(allLimma$logFC),]
head(allLimma)
data=allLimma%>%filter(abs(logFC)>1)%>%filter(adj.P.Val<0.05)

write.table(allLimma,file="GSE157153.txt",sep="\t",quote = F,col.names = T)
write.table(data,file="GSE157153_DEGs.txt",sep="\t",quote = F,col.names = T)
write.table(rownames(data),file="GSE157153_DEGs_GENE.txt",sep="\t",quote = F,col.names = T)
allLimma=read.table(file="GSE157153.txt",header = T,sep = "\t",check.names = F)

#���ð�
# install.packages("VennDiagram")
library(VennDiagram)               #���ð�
outFile="intersectGenes.txt"       #������������ļ�
outPic="venn.pdf"                  #���ͼƬ�ļ�

files=dir()                        #��ȡĿ¼�������ļ�
files=grep("txt$",files,value=T)   #��ȡTXT��β���ļ�
geneList=list()


#��ȡ����txt�ļ��еĻ�����Ϣ�����浽GENELIST              #Ĭ�϶�ȡ��һ��  Ϊ������
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)        #��ȡ
  geneNames=as.vector(rt[,1])              #��ȡ������
  geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
  uniqGene=unique(geneNames)               #����ȡunique
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

head(geneList)


#����venn??
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file=outPic, width=5, height=5)
grid.draw(venn.plot)
dev.off()

a=read.table(file="Yellow_module.txt",header = F,sep = "\t",check.names = F)
b=read.table(file="GSE157153.txt",header = F,sep = "\t",check.names = F)
common=intersect(a$V1,b$V1)

write.table(common,file = "Gene_For_MachineLearning.txt",sep="\t",quote = F,col.names = T,row.names = F)


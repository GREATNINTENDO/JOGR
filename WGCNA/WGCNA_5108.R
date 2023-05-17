library("WGCNA")
library(tidyverse)
rm(list=ls())
AssayData=data.table::fread("exprset_GSE5108.txt")
AssayData=column_to_rownames(AssayData,var = "V1")
WGCNA_matrix = t(AssayData[order(apply(AssayData,1,mad), decreasing = T)[1:20000],])#mad代表绝对中位差
datExpr0 = WGCNA_matrix

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(datExpr0), method = "average")


par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

abline(h = 200, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 190, minSize = 10)
table(clust)
keepSamples = (clust==1)#1保留   0  删除
datExpr0 = datExpr0[keepSamples, ]

#traitData=data.frame(Endometriosis=c(rep(1,10),rep(0,10)),
 #                    Control=c(rep(0,10),rep(1,10)))

traitData=read.table(file="traitdata_5108.txt",header = T,check.names = F)


dim(traitData)
row.names(traitData)=colnames(AssayData)


sameSample=intersect(rownames(datExpr0), rownames(traitData))
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]


sampleTree2 = hclust(dist(datExpr0), method = "average")#去除离群值后再次聚类
plot(sampleTree2)

traitColors = numbers2colors(datTraits, signed = FALSE)
##画图
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


powers = c(c(1:10), seq(from = 12, to = 20, by = 2))     #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)


par(mfrow = c(2,2))#实现一页多图的功能
par(mar = c(4,6,2,10))#下 左 上  右
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
adjacency = adjacency(datExpr0, power =softPower )


TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

gc() 

memory.limit(size = 35000) 

geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


minModuleSize =50      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

gene_count=table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


MEDissThres = 0.15 #剪切高度可修改  MEDissThres =0.15#0.15剪切高度可修改 ###可以完成相似模块的合并,剪切高度是0.15,
#也就是将相似性高于0.85的模块进行了合并
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



moduleColors = mergedColors
table(moduleColors)
write.table(moduleColors,file = "gene_count_module.txt",sep="\t",quote = F,col.names = T)
#moduleColors
#black         blue         cyan        green         grey       grey60 midnightblue         pink 
#1052         3756         5139         5646         2222          177          162          626 
#purple          red 
#454          766 

colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 9, 3, 3))  #mar：以数值向量表示边界大小，顺序为"下、左、上、右"，单位为英分，默认值c(5, 4, 4, 2)+0.1
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#1、计算模块与基因之间的相关性矩阵



modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

a=MMPvalue%>% filter(as.numeric(MMPvalue$p.MMbrown)>0.80)
b=MMPvalue%>% filter(as.numeric(MMPvalue$p.MMmagenta)>0.80)

d=rbind(a,b)
MM=rownames(d)
#计算临床性状与基因之间的相关性矩阵



traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

b=GSPvalue%>% filter(as.numeric(GSPvalue$p.GS.Cancer)>0.20)
GS=rownames(b)

WGCNA_GREEN_GENE=intersect(MM,GS)
write.table(WGCNA_GREEN_GENE,file = "WGCNA_GREEN_GENE.txt",sep="\t",quote = F,col.names = T,row.names = F)

for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}


#最后，我们需要把每个模块里面的基因进行输出，以用于后续的分析。



for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}



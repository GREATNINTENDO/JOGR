library(e1071)
set.seed(12345)
library(tidyverse)
library(parallel)
cl.cores <- detectCores()
cl <- makeCluster(cl.cores)
setwd("D:/Endometriosis_Cancer/07.Machine_Learning")
expreset=read.table(file = "GSE157153_endometriosis_associated.txt",sep = "\t",header = T,check.names = F)
AssayData=column_to_rownames(expreset,var = "gene")

gene=read.table(file = "hubgene.txt",sep = "\t",header = T,check.names = F)
AssayData1=AssayData[gene$gene,]

expreset=as.data.frame(t(AssayData1))
group=read.table(file = "group.txt",sep = "\t",header = T,check.names = F)
expreset$group=group$type
num=as.numeric(ncol(expreset))
expreset=expreset[,c(num,1:(num-1))] #将group放在第一列

expreset$group = ifelse(str_detect(expreset$group,"Cancer"),"1","0") 

expreset$group=as.factor(expreset$group) #行名为样本名，第一列为分组信息，列名为基因名，对应的就是表达情况。
write.table(expreset, 'expreset.txt', sep = '\t', quote = FALSE,row.names = T)
input=expreset
source('msvmRFE.R')  #加载svmrf代码  需要此文件
svmRFE(input, k=10, halve.above=100)

nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)
featsweep = lapply(1:15, FeatSweep.wrap, results, input)
no.info = min(prop.table(table(input[,1])))

errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')

PlotErrors(errors, no.info=no.info)

dev.off()
write.table(top.features[1:4, ], 'svmrfe_top4.txt', sep = '\t', col.names = NA, quote = FALSE)


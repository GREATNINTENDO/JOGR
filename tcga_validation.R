library(data.table)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(rtracklayer)
library(limma)
library(survival)
library(edgeR)
memory.limit(100000000)
#读入表型文件
pd <- fread("TcgaTargetGTEX_phenotype.txt.gz")  #数据集中既有卵巢ca，又有正常卵巢组织，427+88=515
pd <- as.data.frame.matrix(pd)
table(pd[,4])
#甲状腺正常组织样品ID
ov <- pd[pd[,4] == "Ovary",]
ov_all <- as.character(ov$sample)
#读入表达矩阵文件
exp <- fread("TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz")
exp <- as.data.frame.matrix(exp)
#表达矩阵预处理
exp[1:5,1:5]
rownames(exp) <- exp[,1]
exp[1:5,1:5]
exp <- exp[,2:ncol(exp)]
exp[1:5,1:5]

#获取正常组织和肿瘤组织表达矩阵
exp_all <- exp[,colnames(exp) %in% ov_all]
#exp_all=exp
#表达矩阵数值转换
exp_nc <- (2^exp_all - 1)
exp_nc[1:5,1:5]
exp_nc <- round(exp_nc,0)
exp_nc[1:5,1:5]

#####注释ensmbleID#####
expSet=rownames_to_column(exp_nc,var = "gene_id")
expSet=tidyr::separate(expSet,gene_id,into = c("gene_id"),sep="\\.")
gtf<- rtracklayer::import("D:\\代码合集\\gencode.v36.annotation.gtf") 
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df %>% 
  tidyr::separate(gene_id,into = c("gene_id"),sep="\\.") %>% 
  select("gene_id","gene_name")
index <- duplicated(gtf_df[,2])
geneid_df <- gtf_df[! index, ]

expreset <- inner_join(geneid_df, expSet,by = "gene_id")%>%
  dplyr::select(-1)

expreset=column_to_rownames(expreset,var = "gene_name")
write.csv(expreset, "TCGA_GTEx_OV.csv", quote = F, row.names = T)
#DESeq2包进行差异分析
#生成样本信息的factor
condition <- factor(ifelse(substr(colnames(expreset),14,15) == "01"|substr(colnames(expreset),14,15) == "02","Tumor","Control"))
#定义表达矩阵的信息
colData <- data.frame(row.names = colnames(expreset),condition)

gene=c("OCLN","EDNRA")

hub_expr=expreset[gene,]
hub_expr_t=as.data.frame(t(hub_expr))
hub_expr_t=log10(hub_expr_t)
com=cbind(colData,hub_expr_t)
colnames(com)[1]="group"



#设置比较组
data=com

group=levels(factor(data$group))
data$group=factor(data$group, levels=c("Tumor", "Control"))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggpubr)
library(ggsci)
library(ggplot2)
##########################################小提琴图###
library(ggpubr)
ggviolin(data, x="group", y="OCLN", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of OCLN",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="OCLN_TCGA.PDF", width=3.5,height=5)
print(boxplot)
dev.off()


ggviolin(data, x="group", y="EDNRA", fill = "group",
         palette =c("#bc3c29","#0072b5" ),
         add = "boxplot",add.params = list(fill="white"),
         ylab="Expression of EDNRA",
)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="EDNRA_TCGA.PDF", width=3.5,height=5)
print(boxplot)
dev.off()


##############################TCGA-ROC########################################
expr=com
expr$group = ifelse(str_detect(expr$group,"Tumor"),"1","0") 

expr$group=as.numeric(expr$group)

#dat=expr%>%select("group", "OCLN", "EDNRA" )


dat=expr
fit1 <- glm(group ~ OCLN+EDNRA,
            data=dat,
            family = binomial())  
summary(fit1)

dat$prob <- predict(fit1, 
                     newdata=dat, 
                     type="response")
head(dat)


#roc1 <- roc(data$group, data$EPCAM) 
roc2 <- roc(dat$group, dat$OCLN)
roc3 <- roc(dat$group, dat$EDNRA)
#roc4<- roc(data$group, data$MAL)

roc5 <- roc(dat$group, dat$prob)
roc1;
roc2;
roc3;
roc4;
roc5

plot(roc1,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.5, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="skyblue",  # 设置多边形的填充颜色
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # 设置AUC文本的颜色# 设置曲线颜色# 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度

plot(roc2,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.5, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="skyblue",  # 设置多边形的填充颜色
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # 设置AUC文本的颜色# 设置曲线颜色# 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度
plot(roc3,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.5, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="skyblue",  # 设置多边形的填充颜色
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # 设置AUC文本的颜色# 设置曲线颜色# 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度

plot(roc5,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.5, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="skyblue",  # 设置多边形的填充颜色
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # 设置AUC文本的颜色# 设置曲线颜色# 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度

################################################################

mRNA_fpkm_exprset_OV=data.table::fread(file="mRNA_fpkm_exprset_OV.csv")
a=mRNA_fpkm_exprset_OV[1:5,1:5]
expre=column_to_rownames(a,var = "V1")
############################################################################################

expre_survival =data.table::fread('exprset_survival_ov.txt')
a=expre_survival[1:5,1:5]

gene=c("OCLN","EDNRA")

data=expre_survival%>% 
  dplyr::select(c(V1,fustate,futime,"OCLN","EDNRA"))

data=column_to_rownames(data,var = "V1")

library(survival)
library(survminer)

mydata=data
res.cut <- surv_cutpoint(mydata, time = "futime", event = "fustate",#获取TMB最佳cutpoint值，使high lowTMB生存有差异
                         variables = c("OCLN"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
diff=survdiff(Surv(futime,fustate) ~OCLN,data = res.cat)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustate) ~OCLN, data = res.cat)
pdf(file="survival_OCLN_TCGA.pdf",onefile = FALSE,
    width = 4,            
    height =5)            
ggsurvplot(fit, 
           data=res.cat,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=FALSE,
           legend.labs=c("High", "Low"),
           legend.title="OCLN",
           xlab="Time(years)",
           break.time.by = 5,
           risk.table.title="",
           surv.median.line = "hv",####标注中位生存  横线+垂线##
           palette=c("#bc3c29","#0072b5"),
           risk.table.height=.25)
dev.off()

#########################################################################
res.cut <- surv_cutpoint(mydata, time = "futime", event = "fustate",#获取TMB最佳cutpoint值，使high lowTMB生存有差异
                         variables = c("EDNRA"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
diff=survdiff(Surv(futime,fustate) ~EDNRA,data = res.cat)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustate) ~EDNRA, data = res.cat)
pdf(file="survival_EDNRA_TCGA.pdf",onefile = FALSE,
    width = 4,            
    height =5)            
ggsurvplot(fit, 
           data=res.cat,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=FALSE,
           legend.labs=c("High", "Low"),
           legend.title="EDNRA",
           xlab="Time(years)",
           break.time.by = 5,
           risk.table.title="",
           surv.median.line = "hv",####标注中位生存  横线+垂线##
           palette=c("#bc3c29","#0072b5"),
           risk.table.height=.25)
dev.off()

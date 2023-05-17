
#############################机器学习中TOP30的基因###1.三种机器学习算法前30交集 2.AUC最高取前20？ 10？
# 生存信息
meta = read.table("clinical.txt",header = T,row.names = 1)

expreset=read.table(file = "exprSet_GSE26712.txt",sep = "\t",header = T,check.names = F)
expreset=as.data.frame(t(expreset))
symbol=expreset[,c("OCLN","EDNRA")]

comm=intersect(rownames(meta),rownames(symbol))

meta1=meta[comm,]
symbol1=symbol[comm,]
toge=cbind(meta1,symbol1)


expreset_hub_gene=expreset

groupList=factor(c(rep('Normal',10),rep('Cancer',185)))#查看标本分组信息，构建groupList分组
table(groupList)


expreset_hub_gene$Group=groupList
rt=expreset_hub_gene
#设置比较组
rt$Group=factor(rt$Group, levels=c("Normal", "Cancer"))
type=levels(factor(rt[,"Group"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#########################################violin plot#############################################
library(ggpubr)

boxplot=ggboxplot(rt, x="Group", y="TK1", color="Group",
                  xlab="",
                  ylab="Expression of TK1",
                  legend.title="",
                  add = "jitter",
                  palette = c("#0072b5","#bc3c29") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="TK1_gse26712.pdf",width=3.5,height=5)
print(boxplot)
dev.off()

boxplot=ggboxplot(rt, x="Group", y="PKMYT1", color="Group",
                  xlab="",
                  ylab="Expression of PKMYT1",
                  legend.title="",
                  add = "jitter",
                  palette = c("#0072b5","#bc3c29") )+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
pdf(file="PKMYT1_gse26712.pdf",width=3.5,height=5)
print(boxplot)
dev.off()


#GSE26712生存验证
library(survival)
library(survminer)

mydata=toge
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
pdf(file="survival_OCLN_GSE26712.pdf",onefile = FALSE,
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

############
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
pdf(file="survival_EDNRA_GSE26712.pdf",onefile = FALSE,
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



#################################TCGA生存验证########################################
library(survival)
library(survminer)

sur=data.table::fread(file="exprset_survival_ov.txt")
sur=column_to_rownames(sur,var = "V1")
clin=sur[,1:2]
sur[1:5,1:5]
symbol=sur[,c("PKMYT1","TK1")]
toge=cbind(clin,symbol)

mydata=toge
res.cut <- surv_cutpoint(mydata, time = "futime", event = "fustate",#获取TMB最佳cutpoint值，使high lowTMB生存有差异
                         variables = c("TK1"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
diff=survdiff(Surv(futime,fustate) ~TK1,data = res.cat)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustate) ~TK1, data = res.cat)
pdf(file="survival_TK1_tcga.pdf",onefile = FALSE,
    width = 4,            
    height =5)            
ggsurvplot(fit, 
           data=res.cat,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=FALSE,
           legend.labs=c("High", "Low"),
           legend.title="TK1",
           xlab="Time(years)",
           break.time.by = 5,
           risk.table.title="",
           surv.median.line = "hv",####标注中位生存  横线+垂线##
           palette=c("#bc3c29","#0072b5"),
           risk.table.height=.25)
dev.off()


res.cut <- surv_cutpoint(mydata, time = "futime", event = "fustate",#获取TMB最佳cutpoint值，使high lowTMB生存有差异
                         variables = c("PKMYT1"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
diff=survdiff(Surv(futime,fustate) ~PKMYT1,data = res.cat)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustate) ~PKMYT1, data = res.cat)
pdf(file="survival_PKMYT1_tcga.pdf",onefile = FALSE,
    width = 4,            
    height =5)            
ggsurvplot(fit, 
           data=res.cat,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=FALSE,
           legend.labs=c("High", "Low"),
           legend.title="PKMYT1",
           xlab="Time(years)",
           break.time.by = 5,
           risk.table.title="",
           surv.median.line = "hv",####标注中位生存  横线+垂线##
           palette=c("#bc3c29","#0072b5"),
           risk.table.height=.25)
dev.off()



res.cut <- surv_cutpoint(mydata, time = "futime", event = "fustate",#获取TMB最佳cutpoint值，使high lowTMB生存有差异
                         variables = c("SYNE1"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)
diff=survdiff(Surv(futime,fustate) ~SYNE1,data = res.cat)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustate) ~SYNE1, data = res.cat)
pdf(file="survival_SYNE1_tcga.pdf",onefile = FALSE,
    width = 4,            
    height =5)            
ggsurvplot(fit, 
           data=res.cat,
           conf.int=FALSE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=FALSE,
           legend.labs=c("High", "Low"),
           legend.title="SYNE1",
           xlab="Time(years)",
           break.time.by = 5,
           risk.table.title="",
           surv.median.line = "hv",####标注中位生存  横线+垂线##
           palette=c("#bc3c29","#0072b5"),
           risk.table.height=.25)
dev.off()
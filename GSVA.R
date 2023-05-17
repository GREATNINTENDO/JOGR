## 为正常和肿瘤的内皮细胞的基因表达矩阵
library(readxl)
library(dplyr)
a=read.table("symbol.txt",header=T,sep="\t",check.names=F)
#b=read.table("group.txt",header=T,sep="\t",check.names=F)#所有的差异基因

dat <- a
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA") # R 4.1.2 注意版本
library('GSEABase')
library(GSVA)
#geneSets <- getGmt('c2.cp.kegg.v7.4.symbols.gmt')    ###下载的基因集
#geneSets <- getGmt('c5.go.bp.v7.4.symbols.gmt') 
geneSets <- getGmt('h.all.v7.4.symbols.gmt') 
GSVA_hall <- gsva(expr=as.matrix(dat), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=10) # 并行线程数目
head(GSVA_hall)


##########################################LIMMA通路分析##############
## limma
#BiocManager::install('limma')
library(limma)
# 设置或导入分组
group <- factor(c(rep("Control", 37), rep("Cancer", 29)), levels = c('Cancer', 'Control'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design
# Tunor VS Normal
compare <- makeContrasts(Cancer - Control, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)

Diff <- topTable(fit3, coef=1, n = 200)
head(Diff)
#Diff=filter(Diff,adj.P.Val<0.01) #选择上下调前15通路
####################################################################################################
#### 发散条形图绘制 ####
## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
#dat_plot$id <- str_replace(dat_plot$id , "KEGG_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)
library(tidyverse)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, EAOC versus Endometriosis') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签
# 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
ggsave("gsva_bar.pdf",p,width = 8,height  = 8)


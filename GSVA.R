## Ϊ��������������Ƥϸ���Ļ���������
library(readxl)
library(dplyr)
a=read.table("symbol.txt",header=T,sep="\t",check.names=F)
#b=read.table("group.txt",header=T,sep="\t",check.names=F)#���еĲ������

dat <- a
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA") # R 4.1.2 ע��汾
library('GSEABase')
library(GSVA)
#geneSets <- getGmt('c2.cp.kegg.v7.4.symbols.gmt')    ###���صĻ���
#geneSets <- getGmt('c5.go.bp.v7.4.symbols.gmt') 
geneSets <- getGmt('h.all.v7.4.symbols.gmt') 
GSVA_hall <- gsva(expr=as.matrix(dat), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # ����Ϊ��̬�ֲ���T��˫����F
                  kcdf="Gaussian", #CPM, RPKM, TPM���ݾ���Ĭ��ֵ"Gaussian"�� read count������Ϊ"Poisson"��
                  parallel.sz=10) # �����߳���Ŀ
head(GSVA_hall)


##########################################LIMMAͨ·����##############
## limma
#BiocManager::install('limma')
library(limma)
# ���û������
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
#Diff=filter(Diff,adj.P.Val<0.01) #ѡ�����µ�ǰ15ͨ·
####################################################################################################
#### ��ɢ����ͼ���� ####
## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# ȥ��"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
#dat_plot$id <- str_replace(dat_plot$id , "KEGG_","")
# ����һ�� ����t��ֵ����
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# ����
dat_plot <- dat_plot %>% arrange(t)
# �����������
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# ����
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
  ylab('t value of GSVA score, EAOC versus Endometriosis') + #ע����������ת��
  guides(fill=F)+ # ����ʾͼ��
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# ���ӱ�ǩ
# �˴��ο��ˣ�https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# С��-2������
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# С��0������
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# С��2������
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# �ܵ���������
high1 <- nrow(dat_plot)

# ���δ��µ������ӱ�ǩ
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # С��-1��Ϊ��ɫ��ǩ
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # ��ɫ��ǩ
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # ��ɫ��ǩ
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # ����1��Ϊ��ɫ��ǩ
ggsave("gsva_bar.pdf",p,width = 8,height  = 8)

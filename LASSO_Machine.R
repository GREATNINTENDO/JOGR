library(pROC)
library(ggplot2)
library("glmnet")
library(rms)
library(tidyverse)
expreset=read.table("expreset.txt",sep="\t",check.names = F,header = T)
x=expreset[,-1]
y=as.numeric(expreset$group)
model_lasso=glmnet(x,y,nlambda = 10,alpha = 1)
print(model_lasso)
pdf("lasso.lambda.pdf")
plot(model_lasso, xvar = "lambda", label = TRUE)
dev.off()
x=as.matrix(x)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 200,alpha = 1)
plot(cv_fit)
 
#�������߷ֱ�ָʾ����������Ħ�ֵ,һ����lambda.min,һ����lambda.1se,
#������ֵ֮���lambda����Ϊ�Ǻ��ʵġ�lambda.1se������ģ����򵥣�
#��ʹ�õĻ��������٣���lambda.min��׼ȷ�ʸ���һ�㣬ʹ�õĻ�����������һ�㡣
#### 2.2 ����������ֵ���½�ģ
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
#������ֵ�����ڲ���lambda�ϡ�����ģ�ͣ����Խ�ɸѡ�Ļ����������ˡ����л�������ģ�͵��Ӽ�beta�У�
#�õ��Ļ�����һ��s0ֵ��û�õĻ���ֻ��¼�ˡ�.�������Կ�����������������õ��Ļ���

head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)#28
length(choose_gene_1se)#12

write.table(choose_gene_min,file = "LASSO_gene.txt",sep = "\t",row.names = F,quote=F)



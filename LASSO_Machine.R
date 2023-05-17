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
 
#两条虚线分别指示了两个特殊的λ值,一个是lambda.min,一个是lambda.1se,
#这两个值之间的lambda都认为是合适的。lambda.1se构建的模型最简单，
#即使用的基因数量少，而lambda.min则准确率更高一点，使用的基因数量更多一点。
#### 2.2 用这两个λ值重新建模
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
#这两个值体现在参数lambda上。有了模型，可以将筛选的基因挑出来了。所有基因存放于模型的子集beta中，
#用到的基因有一个s0值，没用的基因只记录了“.”，所以可以用下面代码挑出用到的基因。

head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)#28
length(choose_gene_1se)#12

write.table(choose_gene_min,file = "LASSO_gene.txt",sep = "\t",row.names = F,quote=F)




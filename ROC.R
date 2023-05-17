#install.packages("pROC") # 安装包
#install.packages("glmnet") # 建立Logistic回归模型需要
library(pROC) # 加载包
library(glmnet)
library(tidyverse)
expreset=read.table(file='expreset.txt',sep = "\t",header = T,check.names = F)
data=expreset%>%dplyr::select("group", "OCLN", "EDNRA")


fit1 <- glm(group ~ OCLN+EDNRA,
            data=data,
            family = binomial())  
summary(fit1)

data$prob <- predict(fit1, 
                     newdata=data, 
                     type="response")
head(data)





#roc1 <- roc(data$group, data$EPCAM) 
roc2 <- roc(data$group, data$OCLN)
roc3 <- roc(data$group, data$EDNRA)
#roc4<- roc(data$group, data$MAL)

roc5 <- roc(data$group, data$prob)

roc1;
roc2;
roc3;
roc4;

roc5

roc_OCLN <- plot.roc(data$group, data$OCLN,
                 ci=TRUE, print.auc=TRUE) 
roc_EDNRA <- plot.roc(data$group, data$EDNRA,
                     ci=TRUE, print.auc=TRUE) 
roc_prob <- plot.roc(data$group, data$prob,
                      ci=TRUE, print.auc=TRUE) 

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



plot(roc4,  # 前面构建的ROC对象
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

plot(roc6,  # 前面构建的ROC对象
     print.auc=TRUE, # 图像上输出AUC值
     print.auc.x=0.5, print.auc.y=0.5, # AUC值坐标为（x，y）
     auc.polygon=TRUE, # 将ROC曲线下面积转化为多边形
     auc.polygon.col="skyblue",  # 设置多边形的填充颜色
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # 设置AUC文本的颜色# 设置曲线颜色# 不显示网格背景线
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度

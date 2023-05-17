#install.packages("pROC") # ��װ��
#install.packages("glmnet") # ����Logistic�ع�ģ����Ҫ
library(pROC) # ���ذ�
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

plot(roc1,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����

plot(roc2,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����


plot(roc3,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����



plot(roc4,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����

plot(roc5,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����

plot(roc6,  # ǰ�湹����ROC����
     print.auc=TRUE, # ͼ�������AUCֵ
     print.auc.x=0.5, print.auc.y=0.5, # AUCֵ����Ϊ��x��y��
     auc.polygon=TRUE, # ��ROC���������ת��Ϊ�����
     auc.polygon.col="skyblue",  # ���ö���ε������ɫ
     grid= FALSE, 
     col = "red", 
     print.auc.col = "red", # ����AUC�ı�����ɫ# ����������ɫ# ����ʾ���񱳾���
     legacy.axes=TRUE)  # ʹ�����0��1����ʾΪ1-�����
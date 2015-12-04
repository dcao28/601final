setwd("~/601final")
library(xlsx)
set.b<-read.csv(file = "data/Set_b.csv",header = T)
set.c<-read.csv(file = "data/Set_c.csv",header = T)
yy<- read.xlsx("data/TP53.xlsx",sheetIndex = 1)##TP53

View(set.b);View(set.c)
names(set.b)[names(set.b) %in% names(set.c)]

library(MASS)

##1.transformation  
library(bestglm)
library(dplyr)
dd<-mutate(y=log(response),.data = set.b)#new y
dd<- select(dd,-response)
dd<-filter(dd,!(is.infinite(dd$y)))#new data  

##2. linear model selection(AIC,BIC,CP,CV)
aic<-bestglm(Xy = dd,
             family = gaussian,
             IC = "AIC", 
             method = "backward") 
aic$BestModel 

bic<-bestglm(Xy = dd,
             family = gaussian,
             IC = "BIC", 
             method = "forward") 
bic$BestModel

##3. PCR
library(pls)
pca<-pcr(y ~ ., ncomp = 10,scale = TRUE,
    data = dd, validation = "CV")
summary(pca)

##4.PLS
pls<-plsr(y ~ .,ncomp = 10,scale = TRUE,
          data = dd, validation = "CV",method="simpls")
summary(pls)
##5. GMC

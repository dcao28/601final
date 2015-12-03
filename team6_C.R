library(xlsx)
set.b<-read.csv2(file = "data/Set_b.csv",header = T,sep = ",")
set.c<-read.csv2(file = "data/Set_c.csv",header = T,sep = ",")
yy<- read.xlsx("data/TP53.xlsx",sheetIndex = 1)##TP53

View(set.b);View(set.c)
names(set.b)[names(set.b) %in% names(set.c)]

lm1<-lm(response~.,set.b)

drop1(lm1)

rm(list = ls())
gmc.forward<-list.files("/Users/CDX/601final/cdx/Rdata/GMC")
for(i in gmc.forward ){
  cat(i,"\n")
  i<-paste("/Users/CDX/601final/cdx/Rdata/GMC/",i,sep = "")
  load(i)
  print(stepBestGLM)
  cat("##################################################","\n")
}

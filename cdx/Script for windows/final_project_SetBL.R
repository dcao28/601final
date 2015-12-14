
### data 1:6####
rm(list = ls())
suppressMessages(require(parallel))
source("../Script/GMC.R")

load("../data/set_b_log.RData")
set_b_log$y <- exp(set_b_log$y)

# load("../data/set_c_log.RData")
# set_c_log$y <- exp(set_c_log$y)
#load("../data/set_combine_log.RData")

cat("read BL data sets \n")
### link function G1：4 ####  
g1<-function(x){
  return(x)
}

g2<-function(x){
  return(exp(x))
}

g3<-function(x){
  return(x^2)
}

g4<-function(x){
  return(x^3)
}
cat("set link function \n")

### optimiz func 1,2####
findopt1 <- function(lambda,g,model){#with linear assumption
  lambda1<-lambda[1];lambda2<-lambda[2]
  fn<-function(beta){
    xb= as.matrix(model$model[-1]) %*% beta[-1]+beta[1]
    fit.n<-g(xb)
    resid<-model$model$y- fit.n
    vg<-var(fit.n)
    return(-(vg/(vg+var(resid)))+ lambda1*abs(var(fit.n,resid))+lambda2*sum(abs(beta)))
  }
  initials <- model$coefficients
  opt<-optim(par = initials, fn =fn, method = "Nelder-Mead" )
  return(c(lambda,opt$par,-opt$value))
}

findopt2 <- function(lambda,g,model){#no linear assumption
  fn<-function(beta){
    xb= as.matrix(model$model[-1]) %*% beta[-1]+beta[1]
    fit.n<-g(xb)
    resid<-model$model$y- fit.n
    data<-cbind(fit.n,model$model$y)#X,Y
    return(-GMC(data)[2]+lambda*sum(abs(beta)))
  }
  opt<-optim(par = initials, fn =fn, method = "Nelder-Mead" )
  return(c(lambda,opt$par,-opt$value))
}
cat("set optimization function\n")

### initial lambdas ####

lambda2<- rnge<- seq(0,3,by = 0.5)#lasso lambda range
lambda1<- 1 # limit the variance between g(x) and residuals
lambda.set <- rbind(lambda1, lambda2)
lambda.set <- as.data.frame(lambda.set)#;dim(lambda.set)
cat("set lambda range\n")
### compute GMC ####
gmcPerLamda<-function(optgmc,model,g){
  xb<- as.matrix(model$model[-1]) %*% optgmc[c(-1,-2,-3,-length(optgmc))] + optgmc[3]
  fit.n<-g(xb)
  data<-cbind(fit.n, model$model$y)#X,Y
  GMC(data)[2]
}

compute_GMC<-function(lm,g,findopt){
  n.cores = detectCores()
  cluster = makePSOCKcluster(names=n.cores)
  #clusterEvalQ(cl=cluster, expr)
  OptEveryLambda1_b<-parLapply(cl=cluster, X=lambda.set, findopt,g=g,model=lm) #instead of mclapply()
  stopCluster(cl=cluster)
  
  #OptEveryLambda1_b<-mclapply(lambda.set, FUN = findopt,g=g,model=lm,mc.cores=4)
  GmcEverLamda<-sapply(OptEveryLambda1_b, gmcPerLamda,model=lm,g=g)#GMCs
  c(max(GmcEverLamda),OptEveryLambda1_b[[which.max(GmcEverLamda)]])
}# GMC, lambda1,lambda2,betas,opt.value


# compute_GMC<-function(lm,g,findopt){
#   OptEveryLambda1_b<-mclapply(lambda.set, FUN = findopt,g=g,model=lm,mc.cores=4)
#   GmcEverLamda<-sapply(OptEveryLambda1_b, gmcPerLamda,model=lm,g=g)#GMCs
#   c(max(GmcEverLamda),OptEveryLambda1_b[[which.max(GmcEverLamda)]])
# }# GMC, lambda1,lambda2,betas,opt.value
##########################################################

### forward.gmc set_b_log G1 ####
addFormula <- function(add1="") {
  text <- paste("~ .",paste("+", add1, collapse = " "))
  eval(parse(text = text))
}

cat("begin gmc forward selection in Set.B.Log with \n")
pre.m <- lm(y ~ 1, data = set_b_log)
stepBestGLM<- vector(mode="list")
stepBestLM <- vector(mode="list")
gmcs <- c()
cat("step limit is 10 \n")

start<- Sys.time()
for(i in seq(10)){# step of forward selection
  cat("step",i,"\n")
  factors <- attr(terms(as.formula(pre.m)),"factors")
  pre.xy <- rownames(factors)
  if(is.null(pre.xy)){
    avail.x <- setdiff(names(set_b_log),"y")
  }else{avail.x <- setdiff(names(set_b_log),pre.xy)}
  
  new.m <-lapply(avail.x, function(x) update(pre.m,addFormula(add1 = x)))
  #compute_GMC(new.m[[1]],g = g1,findopt = findopt1)
  gmc <- sapply(new.m, FUN = compute_GMC, g=g1,findopt=findopt1)#row one is GMC
  rownames(gmc)[1:3]<-c("GMC","lambda1","lambda2")
  stepBestGLM[[i]] <- gmc[,which.max(gmc[1,])]
  stepBestLM[[i]]<- new.m[[which.max(gmc[1,])]]#lm
  gmcs<-c(gmcs,max(gmc[1,]))
  pre.m <- stepBestLM[[i]]
}
print(stepBestGLM[[which.max(gmcs)]])
end<-Sys.time()
cat("time cost",end-start,"minites")
save(stepBestGLM,file = "../Rdata/GMC/stepBestGlM_Blog_G1.RData")

### forward.gmc set_b_log G2 ####
addFormula <- function(add1="") {
  text <- paste("~ .",paste("+", add1, collapse = " "))
  eval(parse(text = text))
}

cat("begin gmc forward selection in Set.B.Log with \n")
pre.m <- lm(y ~ 1, data = set_b_log)
stepBestGLM<- vector(mode="list")
stepBestLM <- vector(mode="list")
gmcs <- c()
cat("step limit is 10 \n")

start<- Sys.time()
for(i in seq(10)){# step of forward selection
  cat("step",i,"\n")
  factors <- attr(terms(as.formula(pre.m)),"factors")
  pre.xy <- rownames(factors)
  if(is.null(pre.xy)){
    avail.x <- setdiff(names(set_b_log),"y")
  }else{avail.x <- setdiff(names(set_b_log),pre.xy)}
  
  new.m <-lapply(avail.x, function(x) update(pre.m,addFormula(add1 = x)))
  
  gmc <- sapply(new.m, FUN = compute_GMC, g=g2,findopt=findopt1)#row one is GMC
  rownames(gmc)[1:3]<-c("GMC","lambda1","lambda2")
  stepBestGLM[[i]] <- gmc[,which.max(gmc[1,])]
  stepBestLM[[i]]<- new.m[[which.max(gmc[1,])]]#lm
  gmcs<-c(gmcs,max(gmc[1,]))
  pre.m <- stepBestLM[[i]]
}
print(stepBestGLM[[which.max(gmcs)]])
end<-Sys.time()
cat("time cost",end-start,"minites")
save(stepBestGLM,file = "../Rdata/GMC/stepBestGlM_Blog_G2.RData")

### forward.gmc set_b_log G3 ####
addFormula <- function(add1="") {
  text <- paste("~ .",paste("+", add1, collapse = " "))
  eval(parse(text = text))
}

cat("begin gmc forward selection in Set.B.Log with \n")
pre.m <- lm(y ~ 1, data = set_b_log)
stepBestGLM<- vector(mode="list")
stepBestLM <- vector(mode="list")
gmcs <- c()
cat("step limit is 10 \n")

start<- Sys.time()
for(i in seq(10)){# step of forward selection
  cat("step",i,"\n")
  factors <- attr(terms(as.formula(pre.m)),"factors")
  pre.xy <- rownames(factors)
  if(is.null(pre.xy)){
    avail.x <- setdiff(names(set_b_log),"y")
  }else{avail.x <- setdiff(names(set_b_log),pre.xy)}
  
  new.m <-lapply(avail.x, function(x) update(pre.m,addFormula(add1 = x)))
  
  gmc <- sapply(new.m, FUN = compute_GMC, g=g3,findopt=findopt1)#row one is GMC
  rownames(gmc)[1:3]<-c("GMC","lambda1","lambda2")
  stepBestGLM[[i]] <- gmc[,which.max(gmc[1,])]
  stepBestLM[[i]]<- new.m[[which.max(gmc[1,])]]#lm
  gmcs<-c(gmcs,max(gmc[1,]))
  pre.m <- stepBestLM[[i]]
}
print(stepBestGLM[[which.max(gmcs)]])
end<-Sys.time()
cat("time cost",end-start,"minites")
save(stepBestGLM,file = "../Rdata/GMC/stepBestGlM_Blog_G3.RData")

### forward.gmc set_b_log G4 ####
addFormula <- function(add1="") {
  text <- paste("~ .",paste("+", add1, collapse = " "))
  eval(parse(text = text))
}

cat("begin gmc forward selection in Set.B.Log with \n")
pre.m <- lm(y ~ 1, data = set_b_log)
stepBestGLM<- vector(mode="list")
stepBestLM <- vector(mode="list")
gmcs <- c()
cat("step limit is 10 \n")

start<- Sys.time()
for(i in seq(10)){# step of forward selection
  cat("step",i,"\n")
  factors <- attr(terms(as.formula(pre.m)),"factors")
  pre.xy <- rownames(factors)
  if(is.null(pre.xy)){
    avail.x <- setdiff(names(set_b_log),"y")
  }else{avail.x <- setdiff(names(set_b_log),pre.xy)}
  
  new.m <-lapply(avail.x, function(x) update(pre.m,addFormula(add1 = x)))
  
  gmc <- sapply(new.m, FUN = compute_GMC, g=g4,findopt=findopt1)#row one is GMC
  rownames(gmc)[1:3]<-c("GMC","lambda1","lambda2")
  stepBestGLM[[i]] <- gmc[,which.max(gmc[1,])]
  stepBestLM[[i]]<- new.m[[which.max(gmc[1,])]]#lm
  gmcs<-c(gmcs,max(gmc[1,]))
  pre.m <- stepBestLM[[i]]
}
print(stepBestGLM[[which.max(gmcs)]])
end<-Sys.time()
cat("time cost",end-start,"minites")
save(stepBestGLM,file = "../Rdata/GMC/stepBestGlM_Blog_G4.RData")








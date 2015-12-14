lambda2<- rnge<- seq(0,1,by = 0.5)#lasso lambda range
lambda1<- 1e-10 # limit the variance between g(x) and residuals
lambda.set <- rbind(lambda1, lambda2)
lambda.set <- as.data.frame(lambda.set)#;dim(lambda.set)
cat("lambda1 is ",lambda1,"\n")
cat("lambda2 is ",lambda2,"\n")
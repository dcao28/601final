########################################################
#GMC: input a matrix consists of (X,Y), output a vector#
# consists of (GMC(X|Y),GMC(Y|X)).					   #	
########################################################
GMC<-function(data,nlocal=25){#nlocal is used to calculate conditional expectation
	x=data[,1];y=data[,2]
	n=length(x);ln=nlocal

	xdata=data[order(x),];ydata=data[order(y),]
	X=t(matrix(rep(xdata[,1],n),ncol=n)); X=1/(1+abs(X-t(X)))
	Y=t(matrix(rep(ydata[,2],n),ncol=n)); Y=1/(1+abs(Y-t(Y)))
	E_xy=rep(0,n);E_yx=rep(0,n)
	for(i in 1:n){
		li=max(1,i-ln)
		ui=min(n,i+ln)
		E_yx[i]=sum(X[i,li:ui]*xdata[li:ui,2])/sum(X[i,li:ui]); #sum(Y*P_y|xi)
		E_xy[i]=sum(Y[i,li:ui]*ydata[li:ui,1])/sum(Y[i,li:ui])
	}
	
	GMC=c(var(E_xy)/var(x),var(E_yx)/var(y))

	return(GMC)
}
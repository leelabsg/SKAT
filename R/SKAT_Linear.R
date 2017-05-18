
SKAT.linear.Linear = function(res,Z,X1, kernel, weights = NULL, s2, method,res.out,n.Resampling,r.corr, IsMeta=FALSE){
	
	if(length(r.corr) > 1 && dim(Z)[2] == 1){
		r.corr=0
	}
	
	if(IsMeta){
		re = SKAT_RunFrom_MetaSKAT(res=res,Z=Z, X1=X1, kernel=kernel, weights=weights, s2=s2
		, out_type="C", method=method, res.out=res.out, n.Resampling=n.Resampling, r.corr=r.corr)
		
	
	} else if(length(r.corr) == 1 ){

		re = KMTest.linear.Linear(res,Z,X1, kernel, weights, s2, method
		, res.out, n.Resampling, r.corr)

	} else {



		re =SKAT_Optimal_Linear(res, Z, X1, kernel, weights, s2, method
		, res.out, n.Resampling, r.corr)
	}
	return(re)
}


#
#	Modified by Seunggeun Lee - Ver 0.3
#

KMTest.linear.Linear = function(res,Z,X1, kernel, weights, s2, method,res.out,n.Resampling,r.corr){


  n<-nrow(Z)
  # Weighted Linear Kernel 
  if (kernel == "linear.weighted") {
    Z = t(t(Z) * (weights))
  }
  
  if(r.corr == 1){
  	Z<-cbind(rowSums(Z))
  } else if(r.corr > 0){

   	p.m<-dim(Z)[2]	
	R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
	L<-chol(R.M,pivot=TRUE)
	Z<- Z %*% t(L) 
  }


  # get Q
  Q.Temp = t(res)%*%Z

  Q = Q.Temp %*% t(Q.Temp)/s2/2

  Q.res = NULL
  if(n.Resampling > 0){
  	Q.Temp.res = t(res.out)%*%Z

  	Q.res = rowSums(rbind(Q.Temp.res^2))/s2/2
  }

  W.1 = t(Z) %*% Z - (t(Z) %*%X1)%*%solve(t(X1)%*%X1)%*% (t(X1) %*% Z ) # t(Z) P0 Z


  if( method == "liu" ){

	out<-Get_Liu_PVal(Q, W.1, Q.res)    
	pval.zero.msg=NULL

  } else if( method == "liu.mod" ){

	out<-Get_Liu_PVal.MOD(Q, W.1, Q.res)    
	pval.zero.msg = NULL

  } else if( method == "davies"  ){

	out<-Get_Davies_PVal(Q, W.1, Q.res)    
	pval.zero.msg = out$pval.zero.msg

  } else {
	stop("Invalid Method!")
  }


  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling
  , Test.Type = method, Q = Q, param=out$param ,pval.zero.msg=pval.zero.msg )  

  return(re)


}



SKAT.linear.Other = function(res,Z,X1, kernel, weights = NULL, s2, method,res.out,n.Resampling){

  
  n<-nrow(Z)
  m = ncol(Z) 
  if (class(kernel) == "matrix") {

    K = kernel

  } else {

    K = lskmTest.GetKernel(Z, kernel, weights,n,m)

  }


  Q = t(res)%*%K%*%res/(2*s2)

  Q.res = NULL
  if(n.Resampling > 0){
	Q.res<-rep(0,n.Resampling)
	for(i in 1:n.Resampling){

  		Q.res[i] = t(res.out[,i])%*%K%*%res.out[,i]/(2*s2)
  	}
  }


  W = K - X1%*%solve( t(X1)%*%X1)%*%( t(X1) %*% K)	# W = P0 K

  if(method == "davies"){
  	# P0_half = P0
	W1 = W - (W %*% X1) %*%solve( t(X1)%*%X1)%*% t(X1)
  } 


  if( method == "liu" ){

	out<-Get_Liu_PVal(Q, W,Q.res)    
	pval.zero.msg=NULL

  } else if( method == "liu.mod" ){

	out<-Get_Liu_PVal.MOD(Q, W,Q.res)   
	pval.zero.msg = NULL 

  } else if( method == "davies" ){

	out<-Get_Davies_PVal(Q, W1,Q.res)  
	pval.zero.msg = out$pval.zero.msg  

  } else {
	stop("Invalid Method!")
  }

  

  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling
  , Test.Type = method, Q = Q, param=out$param, pval.zero.msg=pval.zero.msg )  

  return(re)

 

}




SKAT.linear = function(Z,y, X1, kernel = "linear", weights = NULL, method="liu", res.out=NULL, n.Resampling = 0, r.corr=r.corr){


  n = length(y) 
  m = ncol(Z) 


  mod = lm(y~X1 -1)

  s2 = summary(mod)$sigma**2

  res = mod$resid



  # If m >> p and ( linear or linear.weight) kernel than call 

  # Linear function

  if( (kernel =="linear" || kernel == "linear.weighted") && n > m){

    re = SKAT.linear.Linear(res,Z,X1, kernel, weights,s2,method,res.out,n.Resampling,r.corr)

  } else {  

    re = SKAT.linear.Other(res,Z,X1, kernel, weights,s2,method,res.out,n.Resampling)  

  }



  return(re)

}






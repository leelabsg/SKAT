
SKAT.logistic.Linear = function(res,Z,X1, kernel, weights = NULL, pi_1, method,res.out,n.Resampling,r.corr, IsMeta=FALSE){

	
	if(length(r.corr) > 1 && dim(Z)[2] == 1){
		r.corr=0
	}

	if(IsMeta){
	
		re = SKAT_RunFrom_MetaSKAT(res=res,Z=Z, X1=X1, kernel=kernel, weights=weights, pi_1=pi_1
		, out_type="D", method=method, res.out=res.out, n.Resampling=n.Resampling, r.corr=r.corr)
	
	} else if(length(r.corr) == 1 ){

		re = KMTest.logistic.Linear(res,Z,X1, kernel, weights, pi_1, method
		, res.out, n.Resampling, r.corr)

	} else {

		
		re =SKAT_Optimal_Logistic(res, Z, X1, kernel, weights, pi_1, method
		, res.out, n.Resampling, r.corr)

	}

	return(re)
}


#
#	Modified by Seunggeun Lee - Ver 0.1
#

KMTest.logistic.Linear = function(res, Z, X1, kernel, weights = NULL, pi_1, method,res.out,n.Resampling,r.corr){

  # Weighted Linear Kernel 
  if (kernel == "linear.weighted") {
    Z = t(t(Z) * (weights))
  }

  # r.corr
  if(r.corr == 1){
  	Z<-cbind(rowSums(Z))
  } else if(r.corr > 0){

   	p.m<-dim(Z)[2]	
	R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
	L<-chol(R.M,pivot=TRUE)
	Z<- Z %*% t(L) 
  }

  # Get temp
  Q.Temp = t(res)%*%Z
  Q = Q.Temp %*% t(Q.Temp)/2

  Q.res = NULL
  if(n.Resampling > 0){
  	Q.Temp.res = t(res.out)%*%Z
  	Q.res = rowSums(rbind(Q.Temp.res^2))/2
  }


  #gg = X1%*%solve(t(X1)%*%(X1 * pi_1))%*%t(X1 * pi_1)  ### Just a holder... not all that useful by itself
  #P0 = D-(gg * pi_1)      ### This is the P0 or P in Zhang and Lin
  # P0 = D-D%*%gg  
  # P0 = D- D%*%X1%*%solve(t(X1)%*%(X1 * pi_1))%*%t(X1) %*% D


 
  #W = P0%*%K
  #W = K * pi_1 - (X1 *pi_1) %*%solve(t(X1)%*%(X1 * pi_1))%*% ( t(X1 * pi_1) %*% K) 
  #muq  = sum(diag(W))/2   # this is the same as e-tilde

  # tr(W W) = tr(P0 K P0 K ) = tr ( Z^T P0 Z Z^T P0 Z ) = tr( P0 Z Z^T P0 Z Z^T )
  # tr(P0 K P0)
  # tr(A B) = tr(A * t(B))
  
  W.1 = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*%X1)%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1)) # t(Z) P0 Z


  if( method == "liu" ){
	out<-Get_Liu_PVal(Q, W.1, Q.res)    
  } else if( method == "liu.mod" ){
	out<-Get_Liu_PVal.MOD(Q, W.1, Q.res)    
  } else if( method == "davies" ){
	out<-Get_Davies_PVal(Q, W.1, Q.res)    
  } else {
	stop("Invalid Method!")
  }


  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling, Test.Type = method, Q = Q,  Q.resampling = Q.res, param=out$param )  
  return(re)
}



SKAT.logistic.Other = function(res, Z, X1, kernel , weights = NULL, pi_1, method,res.out,n.Resampling){
  
  n = nrow(Z) 
  m = ncol(Z)   
  
  # If m >> p and ( linear or linear.weight) kernel than call 
  # Linear function

  if (class(kernel) == "matrix") {
    K = kernel
  } else {
    K = lskmTest.GetKernel(Z, kernel, weights,n,m)
  }


  Q = t(res)%*%K%*%res/2
  Q.res = NULL
  if(n.Resampling > 0){
	Q.res<-rep(0,n.Resampling)
	for(i in 1:n.Resampling){
  		Q.res[i] = t(res.out[,i])%*%K%*%res.out[,i]/2
  	}
  }

  D  = diag(pi_1)   
  gg = X1%*%solve(t(X1)%*%(X1 * pi_1))%*%t(X1 * pi_1)  ### Just a holder... not all that useful by itself
  P0 = D-(gg * pi_1)      ### This is the P0 or P in Zhang and Lin
  # P0 = D-D%*%gg  

  if(method == "davies"){
  	P0_half = Get_Matrix_Square.1(P0)
	#print(dim(P0_half))
	W1 = P0_half %*% K %*% t(P0_half)
  } else {
	#W    = P0%*%K
  	W = K * pi_1 - (X1 *pi_1) %*%solve(t(X1)%*%(X1 * pi_1))%*% ( t(X1 * pi_1) %*% K) 
  }

  if( method == "liu" ){
	out<-Get_Liu_PVal(Q, W, Q.res)    
  } else if( method == "liu.mod" ){
	out<-Get_Liu_PVal.MOD(Q, W, Q.res)   
  } else if( method == "davies" ){
	out<-Get_Davies_PVal(Q, W1, Q.res)    
  } else {
	stop("Invalid Method!")
  }
  
  re<-list(p.value = out$p.value, p.value.resampling = out$p.value.resampling, Test.Type = method, Q = Q, param=out$param )  
  return(re)
 
}



#
#	Modified by Seunggeun Lee - Ver 0.3
#
# 	method : satterth, liu

SKAT.logistic = function(Z,y,X1, kernel = "linear", weights = NULL, method="liu"
, res.out=NULL, n.Resampling = 0, r.corr=r.corr){


	n = length(y) 
	m = ncol(Z) 

	glmfit= glm(y~X1 -1, family = "binomial")
 	betas = glmfit$coef
  	mu    = glmfit$fitted.values
  	eta   = glmfit$linear.predictors

	mu    = glmfit$fitted.values  
	pi_1 = mu*(1-mu)
  	res = y- exp(eta)/(1+exp(eta))

	if(method=="var.match"){
		re = KMTest.logistic.Linear.VarMatching(res, Z, X1, kernel, weights, pi_1, method,res.out,n.Resampling,r.corr, mu)
		return(re)
	}

	# If m >> p and ( linear or linear.weight) kernel than call 
	# Linear function
	if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		re = SKAT.logistic.Linear(res,Z,X1, kernel, weights , pi_1,method,res.out,n.Resampling,r.corr=r.corr)
	} else {  
		re = SKAT.logistic.Other(res,Z,X1, kernel, weights, pi_1, method,res.out,n.Resampling)  
	}

	return(re)
}


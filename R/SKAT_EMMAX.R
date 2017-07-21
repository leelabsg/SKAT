#######################################
#
#	Codes from EMMA

#
#	log restricted log-likelihood function
	
emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
	nq <- length(etas)
	delta <-  exp(logdelta)
	return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}


emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#######################
#
# Changed by SLEE

SKAT.emma.eigen.R.wo.Z <- function(K, X) {
  
  n <- nrow(X)
  q <- ncol(X)
  XX1<-X %*% (solve(crossprod(X,X)))
  K1<-t(X) %*% K
  K2<-K %*% X
  
  #Mat<-K - XX1 %*% K1 - K2 %*% t(XX1) + XX1 %*% ((K1 %*% X) %*% t(XX1)) - XX1 %*% t(X)
  Mat<-K - K2 %*% t(XX1) - XX1 %*% (K1  +  ((K1 %*% X) %*% t(XX1)) -  t(X))
  
  diag(Mat) = diag(Mat) +1
  eig <- eigen(Mat,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#
#	X : covariates
#	Z = NULL
#	Updated by SLEE 07/21/2017
SKAT_NULL_emmaX<- function(formula, data=NULL, K=NULL, Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE) {
  
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n1<-dim(obj2)[1]
	n<-dim(obj1)[1]
	id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)
	X<-model.matrix(formula,data=data)
	X<-Get_SKAT_Residuals.Get_X1(X)
	
	y<-model.response(obj1)
	q <- ncol(X)
	# n and n1 are opposite (compare to SKAT)
	if(n1 - n > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n1 - n)
		warning(MSG,call.=FALSE)
	}

	stopifnot(nrow(X) == n)
	
	########################################################
	# Read kinship 
	
	if(is.null(K) && is.null(Kin.File)){
		stop("both K and Kin.File are NULL!")
	} else if(is.null(K)){
		if(!file.exists(Kin.File)){
			msg<-sprintf("File %s is not exist!", Kin.File)
			stop(msg)
		}
		cat("Read ", Kin.File, "\n")
		K = as.matrix(read.delim(Kin.File, header=FALSE, colClasses="numeric", sep = "\t"))
		t <- nrow(K)
		cat("Read complete. ", Kin.File, " has ", t, " rows! \n")
	}
	
	if(class(K) != "matrix"){
		stop("K is not a matrix!")
	}
	
	if(n1 - n > 0){
		K<-K[id_include, id_include]
	}
	t <- nrow(K)
	stopifnot(ncol(K) == t)
		
	#######################################################
	# Estimate parameters
			
    eig.R <- emma.eigen.R.wo.Z(K,X)
    #eig.R<-SKAT:::SKAT.emma.eigen.R.wo.Z(K,X) # This code sometimes makes an error!
	
    etas <- crossprod(eig.R$vectors,y)
    #etas<-crossprod(eig.R$vectors,res)

    
  	logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      	optlogdelta <- append(optlogdelta, llim)
      	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      	optlogdelta <- append(optlogdelta, ulim)
      	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) ){
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ){
          	r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          	optlogdelta <- append(optlogdelta, r$root)
          	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
	}

	#############################################################################
  	# variance term : maxva K + maxve I
  	# eig.R$vectors $*$ diag(maxva * eig.R$values  + maxdelta) %*% t(eig.R$vectors)
  	
  	
	maxdelta <- exp(optlogdelta[which.max(optLL)])
  	maxLL <- max(optLL)
	maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    # additive effect, genetic

  	maxve <- maxva*maxdelta	# noise
  	 
	#############################################################################
  	# Get V^-1 (y - X \beta)
  	# = V-1 y - V-1 X (X'V-1X)-1 X'V-1 y
  	# = P y
  	# P = V-1 -  V-1 X (X'V-1X)-1 X'V-1

 	va=maxva
 	ve=maxve 

	#if(Is.EPACTS){
 	#	idx.lambda<-which(abs(eig.R$values) > mean(abs(eig.R$values)) / 10000)
 	#  
    #	lambda_inv<-1/(va * eig.R$values[idx.lambda] + ve)
  	#
  	#	P = eig.R$vectors[,idx.lambda] %*% (t(eig.R$vectors[,idx.lambda]) * (lambda_inv))
  	#	res = P %*% y
  	#	

  		
  	V = va * K
  	diag(V) = diag(V) + ve
    	
    # numerical reason
    # added
    if(ve/va < 0.01){
    	eig.V<-eigen(V, symmetric=TRUE)
		idx.lambda<-which(abs(eig.V$values) > mean(abs(eig.V$values)) / 10000)
		V_inv<- eig.V$vectors[,idx.lambda] %*% (t(eig.V$vectors[,idx.lambda]) * (1/eig.V$values[idx.lambda]) )
    	
    } else {
    	V_inv<- solve(V)
  	}
  
  	XVX = t(X) %*% (V_inv %*% X)
  	XVX_inv<-solve(XVX)
  	XV_inv = t(X) %*% V_inv
  	P = V_inv -  t(XV_inv) %*% (XVX_inv %*% XV_inv)
  	res = P %*% y
  	
  	re<-list( LL=maxLL, va=va, ve=ve, P=P, res=res, id_include=id_include)
  	if(Is.GetEigenResult){
  	  re$eig.R=eig.R
  	}
  	
  	class(re)<-"SKAT_NULL_Model_EMMAX"
  	return (re)
}


SKAT_emmaX = function( Z, obj, kernel= "linear.weighted", method="davies", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=max_maf, estimate_MAF=1, SetID = NULL){


	if(class(obj) != "SKAT_NULL_Model_EMMAX"){
		stop("ERROR: obj is not an returned object from SKAT.emmaX.null")
	} 
	#if(method=="optimal.adj"){
	#	stop("SKAT-O is not implemented for SKAT_EmmaX yet!")
	#}
	#if(length(r.corr) > 1){
	#	stop("SKAT-O is not implemented for SKAT_EmmaX yet! r.corr should be a scalar.")
	#}
	

	m = ncol(Z)
	n = nrow(Z)

	# Added by SLEE 4/24/2017
	out.method<-SKAT_Check_Method(method,r.corr, n=n, m=m)
	method=out.method$method
	r.corr=out.method$r.corr
	IsMeta=out.method$IsMeta
	
	if(method=="optimal.adj"){
		IsMeta=TRUE
	}
	
	SKAT_Check_RCorr(kernel, r.corr)
	#
	
	#####################################
	# Check genotypes and parameters
	out.z<-SKAT_MAIN_Check_Z(Z, n, obj$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype
	, is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}
	Z = out.z$Z.test
	weights = out.z$weights
	res = obj$res
	
	if(!IsMeta){
		re = SKAT_emmaX_work(Z=Z, obj=obj, kernel=kernel, method=method, weights=weights, r.corr=r.corr)
	} else {
		re = SKAT_RunFrom_MetaSKAT(res=obj$res,Z=Z, X1=NULL, kernel=kernel, weights=weights
		, P0=obj$P, out_type="V", method=method, res.out=NULL, n.Resampling=0, r.corr=r.corr)

	}

  	re$param$n.marker<-m
  	re$param$n.marker.test<-ncol(Z)
	
	return(re)
}

SKAT_emmaX_work = function( res, Z, obj, kernel, method, weights=NULL, r.corr=0){

	m = ncol(Z)
	n = nrow(Z)
	
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
  	Q.Temp = t(obj$res)%*%Z
  	Q = Q.Temp %*% t(Q.Temp)/2
  	Q.res=NULL

	# Get Z' P0 Z
  	W.1= t(Z) %*% (obj$P %*% Z) # t(Z) P0 Z

 	if( method == "liu.mod" ){

		out<-Get_Liu_PVal.MOD(Q, W.1, Q.res)    

  	} else if( method == "davies" ){

		out<-Get_Davies_PVal(Q, W.1, Q.res)    

  	} else {
		stop("Invalid Method!")
  	}

  
  re<-list(p.value = out$p.value, Test.Type = method, Q = Q, param=out$param ) 

 
  return(re)
}



#
# x is either y or SKAT_NULL_Model 
#
SKAT_emmaX.SSD.OneSet = function(SSD.INFO, SetID, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	Set_Index<-SSD.INFO$SetInfo$SetIndex[id1]

	Z<-Get_Genotypes_SSD(SSD.INFO, Set_Index)
	re<-SKAT_emmaX(Z, obj, ...)
	
	return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKAT_emmaX.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]


	Z<-Get_Genotypes_SSD(SSD.INFO, SetIndex)
	re<-SKAT_emmaX(Z, obj, ...)
	return(re)
}


#
# Only SKAT_Null_Model obj can be used
#
SKAT_emmaX.SSD.All = function(SSD.INFO, obj, ...){
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	Is.Resampling = FALSE
	n.Resampling = 0
	
	pb <- txtProgressBar(min=0, max=N.Set, style=3)
	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1<-try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
		if(class(try1) != "try-error"){
			Z<-try1
			Is.Error<-FALSE
			
			
		} else {
			err.msg<-geterrmessage()
			msg<-sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
		}
	
		if(!Is.Error){
			Is.Error<-TRUE
			try2<-try(SKAT_emmaX(Z,obj, ...),silent = TRUE)
			
			if(class(try2) != "try-error"){
				re<-try2
				Is.Error<-FALSE
			} else {

				err.msg<-geterrmessage()
				msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
				warning(msg,call.=FALSE)
			}
		}
		
		if(!Is.Error){

			OUT.Pvalue[i]<-re$p.value
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
		}
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl)
	class(re)<-"SKAT_SSD_ALL"

	return(re)	
}




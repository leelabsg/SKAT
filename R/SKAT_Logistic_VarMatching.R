

#
#	Modified by Seunggeun Lee - Ver 0.1
#
KMTest.logistic.Linear.VarMatching = function(res, Z, X1, kernel, weights = NULL, pi_1, method, res.out,n.Resampling
, r.corr, mu, res.moments = NULL, Q.sim=NULL){

	n<-length(pi_1)
  	D  = diag(pi_1)   

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

  	Q.Temp = t(res)%*%Z
  	Q = Q.Temp %*% t(Q.Temp)/2

  	Q.res = NULL
  	if(n.Resampling > 0){
  		Q.Temp.res = t(res.out)%*%Z
  		Q.res = rowSums(rbind(Q.Temp.res^2))/2
  	}
  	Z1 = (Z * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1))

	
	if(is.null(Q.sim) && !is.null(res.moments)){

		Q.Temp.res1 = t(cbind(res.moments))%*%Z
  		Q.sim = rowSums(rbind(Q.Temp.res1^2))/2

	} 

	Q.all<-c(Q,Q.res)
	p_all<-mu

	type = "Other"
	if(method =="ECP"){
		type = "OnlySim"
	} else if(method=="QuantileAdj"){
		type = "QuantileAdj"
	}
	
	re<-SKAT_PValue_Logistic_VarMatching(Q.all, Z1 /sqrt(2), p_all, Q.sim, type)

	
	# re$p.value is p-values of aSKAT
	
	p.value.resampling = NULL
	p.value.noadj.resampling = NULL

	
	p.value= re$p.value[1]
	if(length(Q.all) > 1){
		p.value.resampling<-re$p.value[-1]
	}

	p.value.noadj=NULL
	if(!is.null(re$p.value.noadj)){
		p.value.noadj= re$p.value.noadj[1]
		if(length(Q.all) > 1){
	
			p.value.noadj.resampling<-re$p.value.noadj[-1]
		}
	}

	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
, p.value.noadj = p.value.noadj, p.value.noadj.resampling = p.value.noadj.resampling
, Test.Type = method, Q = Q,  Q.resampling = Q.res, param=NULL
, pval.zero.msg=re$pval.msg)  
  
  	return(re)
}




SKAT_PValue_Logistic_VarMatching<-function(Q, Z1, p_all, Q.sim, type="Other"){

	param<-SKAT_Logistic_VarMatching_GetParam1(Z1, p_all, Q.sim, type)
	pval.msg = NULL
	#param1<<-param
	# p-value=1 when varQ==0
	#param1<<-param
	#Q.sim1<<-Q.sim
	#Q1<<-Q
	
	if(param$varQ == 0){
		p.value<-rep(1, length(Q))
	
	} else if(class(param) == "QuantileAdj"){
		p.value<-SKAT_Logistic_PVal_QuantileAdj(Q, param)
	} else {
		Q.Norm<-(Q - param$muQ)/sqrt(param$varQ)
		Q.Norm1<-Q.Norm * sqrt(2*param$df) + param$df
		p.value<- pchisq(Q.Norm1,  df = param$df, ncp=0, lower.tail=FALSE)
		

		df1=param$df
		if(is.null(Q.sim) && param$n.lambda==1){
			re<-Get_Satterthwaite(param$muQ, param$varQ)
			Q.Norm1<-Q / re$a 
			p.value<- pchisq(Q.Norm1,  df = re$df, ncp=0, lower.tail=FALSE)
			df1=re$df
			#cat("df1:", re$df, "\n")
		
		}
		
		if(p.value[1] == 0){
			pval.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q.Norm1, muQ=0, muX=0, sigmaQ=1, sigmaX=1, l=df1, d=0)
	
		}
	}
	

	# SKAT pvalue
	p.value.noadj = NULL
	if(!is.null(param$param.noadj)){
		param.noadj<-param$param.noadj
		
		if(param.noadj$sigmaQ == 0){
			p.value.noadj<-rep(1, length(Q))
		} else {
			Q.Norm<-(Q - param.noadj$muQ)/param.noadj$sigmaQ
			Q.Norm1<-Q.Norm * param.noadj$sigmaX + param.noadj$muX
			p.value.noadj<- pchisq(Q.Norm1,  df = param.noadj$l,ncp=0, lower.tail=FALSE)
		}
	}
	#cat("df2:", param.noadj$l, "\n")
	
	out<-list(p.value=p.value, p.value.noadj=p.value.noadj, param=param, pval.msg=pval.msg)	

	return(out)

}

SKAT_GET_skewness <-  function(x) {
	m3 <- mean((x-mean(x))^3)
	skew <- m3/(sd(x)^3)
	return(skew)
}


SKAT_GET_kurtosis <- function(x) {  

	if(sd(x) == 0){
		return(-100)
	}
	m4 <- mean((x-mean(x))^4) 
	kurt <- m4/(sd(x)^4)-3  
	kurt
}


SKAT_Get_DF_Sim<-function(Q.sim){


	s2.sim<-SKAT_GET_kurtosis(Q.sim)
	df.sim<-12/s2.sim

	if(s2.sim <= 0){
		
		df.sim=100000
	} else if(df.sim < 0.01 ){
		s1.sim<-SKAT_GET_skewness(Q.sim)
		df.sim<-8/s1.sim^2
	}

	return(df.sim)
}

SKAT_Logistic_VarMatching_GetParam1<-function(Z1, p_all, Q.sim, type="Other"){

	type.org = type;
	if(type != "OnlySim"){

		try1<-try(Get_Lambda_U_From_Z(Z1),silent = TRUE)
		if(class(try1) == "try-error"){
			type="OnlySim"
		} else {
			out.svd = try1
			lambda<-out.svd$lambda
			
			if(length(lambda) > 0){
	
    			U<-out.svd$U
				param<-SKAT_Logistic_VarMatching_GetParam(lambda, U, p_all, Q.sim)
			} else {
				type="OnlySim"
			}
		}

	}

	if(type == "OnlySim"){
		param<-SKAT_Logistic_VarMatching_GetParam1_OnlySim(Z1, p_all, Q.sim)
		
	}

	if(type.org == "QuantileAdj"){
		param<-SKAT_Logistic_VarMatching_GetParam1_QuantileAdj(Z1, Q.sim, param)
		return(param)
	}
	
	return(param)

}

####################################################
# in Ver 0.92 for spline adj

SKAT_Logistic_VarMatching_GetParam1_QuantileAdj<-function(Z1,  Q.sim, param){


	#lambda<-Get_Lambda(t(Z1) %*% Z1)
	#muQ = sum(lambda)

	param<-c(param$muQ, param$varQ, param$df)
	
	out.s<-SKAT_Logistic_VarMatching_QuantileAdj_Param(Q.sim, param)
	re<-list(muQ = param$muQ, varQ = param$varQ, df=param$df, out.s=out.s)
	
	class(re)<-"QuantileAdj"
	
	return(re)

}


SKAT_Logistic_VarMatching_QuantileAdj_Param<-function(Q.sim, param){

	#Q.sim<-re.Q1$Q.sim[,1];param<-re.Q1$param;Q<-re.Q1$Q.m[,1]
	
	muQ<-param[1]
	varQ<-param[2]	
	df<-param[3]
	
	if(varQ == 0){
	
		return(NULL)
	}
	
	n1<-length(Q.sim)
	Q.sim<-(Q.sim - muQ)/sqrt(varQ) * sqrt(2 * df) + df
	
	Q.sim.sort<-sort(Q.sim)
	Q.sim.e<-qchisq((1:n1)/(n1 +1), df=df)
	

	n2<-round(n1/20)
	prob<-0:n2/(n2+1) + 10/n1
	Quant1 =  -quantile(-Q.sim.sort, prob)
	Quant2 =  -quantile(-Q.sim.e, prob)
	
	#Q11<<-Quant1
	#Q21<<-Quant2
	#Q.sim.sort1<<-Q.sim.sort
	#Q.sim.e<<-Q.sim.e
	#param1<<-param
	
	fun1<-try(approxfun(Quant1, Quant2), silent=TRUE)
	maxq<-max(Quant1)
	minq<-min(Quant1)
	#ratio<-mean(max(Quant1[1:5]/Quant2[1:5]), Quant[1]/Quant[2])
	ratio<-mean(c(max(Quant1[1:5]/Quant2[1:5]), Quant1[1]/Quant2[1]))
	
	if(class(fun1) == "try-error"){
		out.s = NULL
	} else {
		out.s<-list(fun1=fun1, maxq=maxq, minq=minq, ratio=ratio)
	}
	
	return(out.s)

}

####################################################
# in Ver 0.92 for spline adj

SKAT_Logistic_PVal_QuantileAdj<-function(Q, param.s){

	# PValues
	Q.norm<-(Q - param.s$muQ)/sqrt(param.s$varQ) * sqrt(2*param.s$df) + param.s$df
	out.s<-param.s$out.s
	if(is.null(out.s) == FALSE){
		
		pval<-rep(NA, length(Q.norm))
		
		IDX1<-which(Q.norm <= out.s$minq)
		IDX2<-which(Q.norm >= out.s$maxq)
		IDX3<-setdiff(setdiff(1:length(Q.norm), IDX1), IDX2)
	
		if(length(IDX1) > 0){
			pval[IDX1]<-pchisq(Q.norm[IDX1],df=param.s$df, lower.tail = FALSE) 
		}
		if(length(IDX2) > 0){
			pval[IDX2]<-pchisq(Q.norm[IDX2]/out.s$ratio,df=param.s$df, lower.tail = FALSE) 
		}
		if(length(IDX3) > 0){
			temp<-out.s$fun1(Q.norm[IDX3])
			pval[IDX3]<-pchisq(temp,df=param.s$df, lower.tail = FALSE) 
		}
	
		#Q.norm1<<-Q.norm
		#out1.1<<-out.s
		#param.s1<<-param.s
		#pval1<<-pval
	} else {
		pval<-pchisq(Q.norm,df=param.s$df, lower.tail = FALSE) 
	}

	return(pval)


}




################################
# ver 0.82 changed
SKAT_Logistic_VarMatching_GetParam1_OnlySim<-function(Z1, p_all, Q.sim){

	#out.svd = Get_Lambda_U_From_Z(Z1)
	#lambda<-out.svd$lambda
	
	lambda<-Get_Lambda(t(Z1) %*% Z1)


	muQ = sum(lambda)
	varQ.sim<-var(Q.sim)

	#print(c(varQ, varQ.sim))
	df.sim<-SKAT_Get_DF_Sim(Q.sim)

	# No adjustment
	c1<-rep(0,4)	
	for(i in 1:4){
		c1[i]<-sum(lambda^i)
	}	
	param<-Get_Liu_Params_Mod(c1)

	return(list(muQ = muQ, varQ = varQ.sim, df=df.sim, lambda.new=NULL, param.noadj = param, n.lambda=0))

}


#
#	lambda : eigenvalues
#	U : eigenvectors
#	p_all : 
#	If Q.resample == NULL, it does not estimate kurtosis and df.sim < 0
SKAT_Logistic_VarMatching_GetParam<-function(lambda, U, p_all, Q.sim){

	# Var match
	re<-SKAT_Get_Cov_Param(lambda, p_all, U)

	# New Lambda 
	lambda.new<- re$lambda.new
	
	# new parameters
	muQ<-re$muQ

	# new var
	varQ<-re$varQ

	# df
	s2 = sum(lambda.new^4) / sum(lambda.new^2)^2
	df<-1/s2
		

	if(!is.null(Q.sim)){
		df<-SKAT_Get_DF_Sim(Q.sim)
	}

	
	# No adjustment
	c1<-rep(0,4)	
	for(i in 1:4){
		c1[i]<-sum(lambda^i)
	}	
	param<-Get_Liu_Params_Mod(c1)

	return(list(muQ = muQ, varQ = varQ, df=df, lambda.new=lambda.new, param.noadj = param, n.lambda=length(lambda.new)))

}


SKAT_Get_Var_Elements<-function(m4,p_all,u1,u2){

	temp1<-u1^2 * u2^2


	a1<-sum(m4 * temp1)
	a2<-sum(u1^2) * sum(u2^2) - sum(temp1)
	a3<-sum(u1*u2)^2 - sum(temp1)

	a3<-a3*2
	
	
	a1+a2+a3 
}


SKAT_Get_Cov_Param<-function(lambda,p_all,U){

	#p_all<-obj$mu
	#U<-out$U
	#lambda<-out$lambda

	#cat("1", dim(U), "lambda:", lambda, "\n")
	
	p.m<-length(lambda)
	m4<-p_all*(1-p_all)*(3*p_all^2-3*p_all +1) / (p_all*(1-p_all))^2
		
	zeta<-rep(0,p.m)
	var_i<-rep(0,p.m)
	varQ<-0	

	for(i in 1:p.m){
		temp.M1<-sum(U[,i]^2)^2 - sum(U[,i]^4) 
		zeta[i]<-sum(m4 * U[,i]^4) + 3* temp.M1 # because ( \sum .)^4, not ^2
		var_i[i]<-zeta[i] - 1	
	}

	if(p.m == 1){
		Cov_Mat<-matrix(zeta* lambda^2, ncol=1,nrow=1)
	} else if(p.m > 1){

		Cov_Mat<-diag(zeta* lambda^2)
		for(i in 1:(p.m-1)){
			for(j in (i+1):p.m){
				Cov_Mat[i,j]<-SKAT_Get_Var_Elements(m4,p_all,U[,i],U[,j])
				Cov_Mat[i,j]<-Cov_Mat[i,j]* lambda[i]* lambda[j]
			}
		}
	} else{
		msg<-sprintf("Error SKAT_Get_Cov_Param: p.m=%d \n", p.m)
		stop(msg)
	
	}
	
	Cov_Mat<-Cov_Mat + t(Cov_Mat)
	diag(Cov_Mat)<-diag(Cov_Mat)/2

	varQ<-sum(Cov_Mat) - sum(lambda)^2	
	muQ=sum(lambda)
	lambda.new<-lambda * sqrt(var_i)/sqrt(2)
	
	
	return(list(zeta=zeta, var_i=var_i, varQ = varQ, muQ=muQ, lambda.new=lambda.new, n.lambda=p.m))
	
}




#
#	Function get parameters of optimal test
#
SKAT_Optimal_Param<-function(Z1,r.all){


	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]
	r.n<-length(r.all)

	z_mean<-rowMeans(Z1)
	Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
	cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

	Z.item1<-Z_mean %*% diag(cof1)	
	Z.item2<-Z1 - Z.item1

	# W3.2 Term : mixture chisq
	W3.2.t<-t(Z.item2) %*% Z.item2
	lambda<-Get_Lambda(W3.2.t)
	
	# W3.3 Term : variance of remaining ...
	W3.3.item<-sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
	
	# Mixture Parameters
	MuQ<-sum(lambda)
	VarQ<-sum(lambda^2) *2 + W3.3.item
	KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
	Df<-12/KerQ

	# W3.1 Term : tau1 * chisq_1
	tau<-rep(0,r.n)
	for(i in 1:r.n){
		r.corr<-r.all[i]
		#term1<-p.m*r.corr + cof1^2 * (1-r.corr)
		term1<-p.m^2*r.corr + sum(cof1^2) * (1-r.corr)
		tau[i]<-sum(term1) *  sum(z_mean^2)
	}

	out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau)
	return(out)
}



#
#	Function get SKAT statistics with given rho
#		Q.all is a matrix with n.q x n.r

SKAT_Optimal_Each_Q<-function(param.m, Q.all, r.all, lambda.all, method=NULL){

	n.r<-length(r.all)
	c1<-rep(0,4)
	n.q<-dim(Q.all)[1]

	pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
	param.mat<-NULL

	for(i in 1:n.r){
		Q<-Q.all[,i]
		r.corr<-r.all[i]
		lambda.temp<-lambda.all[[i]] 
		c1[1]<-sum(lambda.temp)
		c1[2]<-sum(lambda.temp^2)
		c1[3]<-sum(lambda.temp^3)
		c1[4]<-sum(lambda.temp^4)
		param.temp<-Get_Liu_Params_Mod(c1)

		muQ<-param.temp$muQ
		varQ<-param.temp$sigmaQ^2
		df<-param.temp$l

		# get pvalue
		Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
		pval[,i]<- pchisq(Q.Norm,  df = df, lower.tail=FALSE)
		# will be changed later
		
		if(!is.null(method)){
			if(method=="optimal.mod" || method=="optimal.adj" || method=="optimal.moment.adj" ){
				pval[,i]<-Get_PValue.Lambda(lambda.temp,Q)$p.value
			}
		}
		
		param.mat<-rbind(param.mat,c(muQ,varQ,df))
	}

	pmin<-apply(pval,1,min)
	for(i in 1:n.r){
	
		muQ<-param.mat[i,1]
		varQ<-param.mat[i,2]
		df<-param.mat[i,3]

		q.org<-qchisq(1-pmin,df=df)
		q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
		pmin.q[,i]<-q.q

	}
	
	out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
	return(out)

}

SKAT_Optimal_Integrate_Func_Davies<-function(x,pmin.q,param.m,r.all){
	
	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	re<-rep(0,length(x))
	for(i in 1:length(x)){
		#a1<<-temp.min[i]
		min1<-temp.min[i]
		if(min1 > sum(param.m$lambda) * 10^4){
			temp<-0
		} else {
			min1.temp<- min1 - param.m$MuQ			
			sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
			min1.st<-min1.temp *sd1 + param.m$MuQ
			
			dav.re<-SKAT_davies(min1.st,param.m$lambda,acc=10^(-6))
			temp<-dav.re$Qq
			if(dav.re$ifault != 0){
				stop("dav.re$ifault is not 0")
			}
		}
		if(temp > 1){
			temp=1
		}
		#lambda.record<<-param.m$lambda
		#print(c(min1,temp,dav.re$ifault,sum(param.m$lambda)))
		re[i]<-(1-temp) * dchisq(x[i],df=1)
	}
	return(re)

}

# add pmin on 02-13-2013 
SKAT_Optimal_PValue_Davies<-function(pmin.q,param.m,r.all, pmin=NULL){

	#re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15), silent = TRUE)

	re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25), silent = TRUE)

	if(class(re) == "try-error"){
		re<-SKAT_Optimal_PValue_Liu(pmin.q,param.m,r.all, pmin)
		return(re)
	} 

	pvalue<-1-re[[1]]
	if(!is.null(pmin)){
		if(pmin *length(r.all) < pvalue){
			pvalue = pmin *length(r.all)
		}
	}
	
	
	return(pvalue)

}


SKAT_Optimal_Integrate_Func_Liu<-function(x,pmin.q,param.m,r.all){
	
	#x<-1
	#print(length(x))
	#print(x)
	#X1<<-x
	#x<-X1

	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	temp.q<-(temp.min - param.m$MuQ)/sqrt(param.m$VarQ)*sqrt(2*param.m$Df) + param.m$Df
	re<-pchisq(temp.q ,df=param.m$Df) * dchisq(x,df=1)
	
	return(re)

}

# add pmin on 02-13-2013 
SKAT_Optimal_PValue_Liu<-function(pmin.q,param.m,r.all, pmin=NULL){

	 re<-integrate(SKAT_Optimal_Integrate_Func_Liu, lower=0, upper=40, subdivisions=2000
	,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25)
	
	pvalue<-1-re[[1]]
	
	if(!is.null(pmin)){
		if(pmin *length(r.all) < pvalue){
			pvalue = pmin *length(r.all)
		}
	}
	
	return(pvalue)

}


SKAT_Optimal_Get_Q<-function(Z1, res, r.all, n.Resampling, res.out, res.moments=NULL, Q.sim=NULL){

	n.r<-length(r.all)
	p.m<-dim(Z1)[2]

	Q.r<-rep(0,n.r)
	Q.r.res<-NULL
	
	temp<-t(res) %*% Z1
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q1<-(1-r.corr) * rowSums(temp^2)
		Q2<-r.corr * p.m^2 * rowMeans(temp)^2
		Q.r[i]<-Q1 + Q2
	}
	Q.r = Q.r /2
  	if(n.Resampling > 0){
	
		temp<-t(res.out) %*% Z1
		Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp^2)
			Q2<-r.corr * p.m^2 * rowMeans(temp)^2
			Q.r.res[,i]<-Q1 + Q2
		}
		Q.r.res = Q.r.res/2
  	}

	if(!is.null(res.moments) && is.null(Q.sim)){

		temp<-t(res.moments) %*% Z1
		n.moments<-dim(res.moments)[2]
		Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp^2)
			Q2<-r.corr * p.m^2 * rowMeans(temp)^2
			Q.sim[,i]<-Q1 + Q2
		}
		Q.sim = Q.sim/2

	}

	re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.sim=Q.sim)
	return(re)


}

SKAT_Optimal_Get_Pvalue<-function(Q.all, Z1, r.all, method){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Z1)[2]

	lambda.all<-list()
	for(i in 1:n.r){
		r.corr<-r.all[i]
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z2<- Z1 %*% t(L)
		K1<-t(Z2) %*% Z2

		lambda.all[[i]]<-Get_Lambda(K1)
		 
	}

	# Get Mixture param 
	param.m<-SKAT_Optimal_Param(Z1,r.all)
	Each_Info<-SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
	pmin.q<-Each_Info$pmin.q
	pmin<-Each_Info$pmin
	pval<-rep(0,n.q)

	if(method == "davies" || method=="optimal" || method=="optimal.mod" || method=="optimal.adj"){

		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])
			
		}


	} else if(method =="liu" || method =="liu.mod" || method=="optimal.moment" || method=="optimal.moment.adj" ){
		
		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all, pmin[i])
		}

	} else {
		stop("Invalid Method!")
	}
	
	# Check the pval 
	# Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
	# To correct conservatively, we use min(p-values) * 3
	
	multi<-3
	if(length(r.all) < 3){
		multi<-2
	}

	for(i in 1:n.q){
		pval.each<-Each_Info$pval[i,]
		IDX<-which(pval.each > 0)
		
		pval1<-min(pval.each) * multi
		if(pval[i] <= 0 || length(IDX) < length(r.all)){
			pval[i]<-pval1
		}
		
		# if pval==0, use nonzero min each.pval as p-value
		if(pval[i] == 0){
			if(length(IDX) > 0){
				pval[i] = min(pval.each[IDX])
			}
		}
	
	}
	
	return(list(p.value=pval,p.val.each=Each_Info$pval))

}

#######################################################33
#	Linear

SKAT_Optimal_Linear = function(res,Z,X1, kernel, weights = NULL, s2,method=NULL
, res.out=NULL, n.Resampling =0, r.all){

	# if r.all >=0.999 ,then r.all = 0.999. It is just for computation.
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999	
	}

	n<-dim(Z)[1]
	p.m<-dim(Z)[2]	
	n.r<-length(r.all)
	
	if (kernel == "linear.weighted") {
		Z = t(t(Z) * (weights))
	}

  	#Get P0Z, wher P0 = diag(n) - X1%*%solve( t(X1)%*%X1)%*%t(X1)
	Z1<-Z - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z)

	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-SKAT_Optimal_Get_Q(Z, res, r.all, n.Resampling, res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2

	##################################################
	# Compute P-values 
	#################################################

	out<-SKAT_Optimal_Get_Pvalue(Q.all, Z1 / sqrt(2), r.all, method)
	
	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each[1,]
	param$q.val.each<-Q.all[1,]
	param$rho<-r.all
	param$minp<-min(param$p.val.each)

	id_temp<-which(param$p.val.each == min(param$p.val.each))
	id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
	if(length(id_temp1) > 0){
		param$rho[id_temp1] = 1
	}
	param$rho_est<-param$rho[id_temp]


	p.value<-out$p.value[1]

	p.value.resampling<-NULL
	if(n.Resampling > 0){
		p.value.resampling<-out$p.value[-1]
		#param$pval.each.resample<-out$p.val.each[-1]
	}

 	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
	, Test.Type = method, Q = NA, param=param )  
  	
	return(re)	

}

#######################################################
#	Logistic 
#		pi_1 = mu *  (1-mu) : variance

SKAT_Optimal_Logistic  = function(res, Z, X1, kernel, weights = NULL, pi_1 , method = NULL
, res.out=NULL, n.Resampling =0, r.all){

	# if r.all >=0.999 ,then r.all = 0.999
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999	
	}

	n<-dim(Z)[1]
	p.m<-dim(Z)[2]
	n.r<-length(r.all)
	
	D  = diag(pi_1)   
	if (kernel == "linear.weighted") {
		Z = t(t(Z) * (weights))
	}

          
  	Z1 = (Z * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1))

	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-SKAT_Optimal_Get_Q(Z, res, r.all, n.Resampling, res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) 

	##################################################
	# Compute P-values 
	#################################################

	out<-SKAT_Optimal_Get_Pvalue(Q.all, Z1 / sqrt(2), r.all, method)

	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each[1,]
	param$q.val.each<-Q.all[1,]
	param$rho<-r.all
	param$minp<-min(param$p.val.each)

	id_temp<-which(param$p.val.each == min(param$p.val.each))
	id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
	if(length(id_temp1) > 0){
		param$rho[id_temp1] = 1
	}

	param$rho_est<-param$rho[id_temp]


	p.value<-out$p.value[1]
	p.value.resampling<-NULL
	if(n.Resampling > 0){
		p.value.resampling<-out$p.value[-1]
		#param$pval.each.resample<-out$p.val.each[-1]
	}

 	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
	, Test.Type = method, Q = NA, param=param )  
  	
	return(re)	

}



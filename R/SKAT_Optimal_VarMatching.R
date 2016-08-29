SKAT_Optimal_Get_Kertosis_Mixture<-function(df1, df2, v1, a1, a2){


	v2<-2*df2

	S4.1<-(12/df1 +3) * v1^2
 	S4.2<-(12/df2 +3) * v2^2

	#v1<-2*df1 + var.add
	

	S4<-a1^4*S4.1 + a2^4*S4.2 + 6 * a1^2 * a2^2 * v1 * v2
	S2<-a1^2*v1 + a2^2*v2
	
	K<-S4/(S2^2) - 3

	if(K < 0){
		K<-0.0001
	}

	
	#print(c(S4.1, S4.2, v1, v2, K, a1, a2, df1,df2))
	return(K)

}


SKAT_Optimal_GetZItem2<-function(Z1){

	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]

	z_mean<-rowMeans(Z1)
	Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
	cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

	Z.item1<-Z_mean %*% diag(cof1)	
	Z.item2<-Z1 - Z.item1

	return(Z.item2)

}

#
#	Function get parameters of optimal test
#
SKAT_Optimal_Param_VarMatching<-function(Z1, r.all, p_all, res.moments, method="Other", Q.sim=NULL){


	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]
	r.n<-length(r.all)

	z_mean<-rowMeans(Z1)
	Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
	cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

	Z.item1<-Z_mean %*% diag(cof1)	
	Z.item2<-Z1 - Z.item1

	# W3.2.t<-t(Z.item2) %*% Z.item2 follows mixture of chisq distribution
	# apply adjustment

	if(!is.null(res.moments) && is.null(Q.sim)){

		Q.Temp.res1 = t(cbind(res.moments))%*%Z.item2
  		Q.sim = rowSums(rbind(Q.Temp.res1^2))/2

	} 

	type = "Other"
	if(method == "ECP"){
		type = "OnlySim"
	} else if(method=="QuantileAdj"){
		type = "QuantileAdj"
	}


	re.param<-SKAT_Logistic_VarMatching_GetParam1(Z.item2, p_all, Q.sim, type)
	
	# W3.3 Term : variance of remaining ...
	W3.3.item<-sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
	
	# W3.1 Term : tau1 * chisq_1
	tau<-rep(0,r.n)
	for(i in 1:r.n){
		r.corr<-r.all[i]
		term1<-p.m*r.corr + cof1^2 * (1-r.corr)
		tau[i]<-sum(term1) *  sum(z_mean^2)
	}

	out<-list(param=re.param, VarRemain=W3.3.item, tau=tau)
	return(out)
}


#
#	Function get SKAT statistics with given rho
#		Q.all is a matrix with n.q x n.r

SKAT_Optimal_Each_Q_VarMatching<-function(param.m, Q.all, r.all, Z2.all, p_all, Q.sim.all, method="Other"){

	type = "Other"
	if(method == "ECP"){
		type = "OnlySim"
	} else if(method=="QuantileAdj"){
		type = "QuantileAdj"
	}

	Is.SIM<-!is.null(Q.sim.all)
	n.r<-length(r.all)
	c1<-rep(0,4)
	n.q<-dim(Q.all)[1]

	pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pval.sim<-matrix(rep(0,n.r*n.q),ncol=n.r)
	
	pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pmin.q.sim<-matrix(rep(0,n.r*n.q),ncol=n.r)

	re.param<-list()
	for(i in 1:n.r){
		Q<-Q.all[,i]
		r.corr<-r.all[i]
		
		out<-SKAT_PValue_Logistic_VarMatching(Q, Z2.all[[i]], p_all, Q.sim.all[,i],type)

		re.param[[i]]<-out$param
		pval[,i]<- out$p.value

	}

	#pval1<<-pval
	#re.param1<<-re.param
	pmin<-apply(pval,1,min)

	# re-adjust the kurtosis of Q using the estimated kurtosis
	for(i in 1:n.r){
		r.corr<-r.all[i]
		muQ<-re.param[[i]]$muQ
		varQ1<-re.param[[i]]$varQ
		varQ<-(1-r.corr)^2*(param.m$param$varQ + param.m$VarRemain) + param.m$tau[i]^2*2
		
		#print(c(varQ, varQ1, param.m$tau[i]^2*2))

		df<-re.param[[i]]$df

		vq1<-param.m$param$varQ + param.m$VarRemain
		ker<-SKAT_Optimal_Get_Kertosis_Mixture(param.m$param$df, 1, vq1 , (1-r.corr), param.m$tau[i])

		
		df<-12/ker

		q.org<-qchisq(1-pmin,df=df)
		q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
		pmin.q[,i]<-q.q

	}
	out<-list(pmin=pmin, pval=pval, pmin.q=pmin.q)
	return(out)

}



SKAT_Optimal_Integrate_Func_VarMatching<-function(x, pmin.q, muQ, varQ, df, tau, r.all){
	

	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	if(varQ > 0){
		temp.q<-(temp.min - muQ)/sqrt(varQ)*sqrt(2*df) + df
		val=pchisq(temp.q ,df=df)
	} else {
		val=rep(0, length(temp.q))
	}
	
	re<-val * dchisq(x,df=1)

	#df.x<-1
	#x.norm<-(x -1)/sqrt(2) * sqrt(2*df.x) + df.x
	#re<-pchisq(temp.q ,df=df) * dchisq(x.norm,df=df.x)
	return(re)

}


SKAT_Optimal_PValue_VarMatching<-function(pmin.q, muQ, varQ, df, tau, r.all, pmin=NULL){
	
	
	re<-integrate(SKAT_Optimal_Integrate_Func_VarMatching, lower=0, upper=40, subdivisions=2000, pmin.q=pmin.q, muQ=muQ, varQ=varQ, df=df
	, tau=tau, r.all=r.all, abs.tol = 10^-25)
	
	pvalue<-1-re[[1]]
	
	# added Jan 31, 2014
	if(!is.null(pmin)){
		if(pmin *length(r.all) < pvalue){
			pvalue = pmin *length(r.all)
		}
	}
	
	return(pvalue)

}



SKAT_Optimal_Get_Pvalue_VarMatching<-function(Q.all, Z1, r.all, p_all, Q.sim.all, res.moments, method=NULL, Q.sim=NULL){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Z1)[2]

	lambda.all<-list()
	Z2.all<-list()
	for(i in 1:n.r){
		r.corr<-r.all[i]
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z2.all[[i]]<- Z1 %*% t(L)
	}

	# Get Mixture param 
	param.m<-SKAT_Optimal_Param_VarMatching(Z1,r.all,p_all, res.moments,method=method, Q.sim=Q.sim)

	Each_Info<-SKAT_Optimal_Each_Q_VarMatching(param.m, Q.all, r.all, Z2.all,p_all, Q.sim.all,method)
	pmin.q<-Each_Info$pmin.q
	pmin.q.sim<-Each_Info$pmin.q.sim
	pmin<-Each_Info$pmin
	
	pval<-rep(0,n.q)
	pval.sim<-rep(0,n.q)

	muQ 	= param.m$param$muQ
	varQ 	= param.m$param$varQ + param.m$VarRemain
	df 	= param.m$param$df
	tau 	= param.m$tau

	#
	# We only calculate one types of p-values to save computing time 
	#
	p.val.each=Each_Info$pval
	for(i in 1:n.q){
		# there was bug in this part, and fixed it
		pval[i]<-SKAT_Optimal_PValue_VarMatching(pmin.q[i,], muQ, varQ, df, tau, r.all, pmin=pmin[i])
	}
	
	# Check the pval 
	# If there is any Each_Info$pval ==0, it does not work properly.
	# Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
	# To correct conservatively, we use min(p-values) * 3 when number(r.all) >= 3
	
	
	multi<-3
	if(length(r.all) < 3){
		multi<-2
	}
	
	for(i in 1:n.q){
		pval.each<-Each_Info$pval[i,]
		IDX<-which(pval.each > 0)
		
		pval1<-min(pval.each) * multi
		if(pval[i] < 0 || length(IDX) < length(r.all)){
			pval[i]<-pval1
		}
		
		# if pval==0, use nonzero min each.pval as p-value
		if(pval[i] == 0){
			if(length(IDX) > 0){
				pval[i] = min(pval.each[IDX])
			}
		}
	}
	
	
	
	return(list(p.value=pval, p.val.each=p.val.each))

}


SKAT_Optimal_Logistic_VarMatching  = function(res, Z, X1, kernel, weights = NULL, pi_1 , method = NULL
, res.out=NULL, n.Resampling =0, r.all, mu, res.moments = NULL, Q.sim=NULL, Q.sim.a=NULL){

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
	out.Q<-SKAT_Optimal_Get_Q(Z, res, r.all, n.Resampling, res.out, res.moments, Q.sim.a)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) 
	Q.sim.all<-out.Q$Q.sim

	##################################################
	# Compute P-values 
	#################################################

	p_all<-mu
	out<-SKAT_Optimal_Get_Pvalue_VarMatching(Q.all, Z1 / sqrt(2), r.all, p_all, Q.sim.all, res.moments, method=method, Q.sim=Q.sim)

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



	p.value.resampling = NULL
	p.value= out$p.value[1]

	if(n.Resampling > 1){
		p.value.resampling<-out$p.value[-1]
	}
	

 	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
	, Test.Type = "moments.matching", Q = NA, param=param )  
  	
	return(re)	

}




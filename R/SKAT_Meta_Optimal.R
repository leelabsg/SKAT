
SKAT_META_Optimal_Get_Q<-function(Score, r.all){

	n.r<-length(r.all)
	Q.r<-rep(0,n.r)
	
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q.r[i]<-(1-r.corr) * sum(Score^2) + r.corr * sum(Score)^2

	}
	Q.r = Q.r /2
 
	re<-list(Q.r=Q.r)
	return(re)


}


SKAT_META_Optimal_Get_Q_Res<-function(Score.res, r.all){

	n.r<-length(r.all)
	p<-dim(Score.res)[1]
	Q.r<-matrix(rep(0,n.r*p), ncol=n.r)
	
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q.r[,i]<-(1-r.corr) * rowSums(Score.res^2) + r.corr * rowSums(Score.res)^2

	}
	Q.r = Q.r /2
 
	re<-list(Q.r=Q.r)
	return(re)


}

SKAT_META_Optimal_Get_Pvalue<-function(Q.all, Phi, r.all, method){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Phi)[2]

	lambda.all<-list()
	for(i in 1:n.r){
		r.corr<-r.all[i]
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Phi_rho<- L %*% (Phi %*% t(L))
		lambda.all[[i]]<-Get_Lambda(Phi_rho)
		 
	}

	# Get Mixture param 
	param.m<-SKAT_META_Optimal_Param(Phi,r.all)
	Each_Info<-SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
	pmin.q<-Each_Info$pmin.q
	pval<-rep(0,n.q)
	
	# added
	pmin<-Each_Info$pmin

	if(method == "davies" || method=="optimal" ||  method=="optimal.adj" || method=="optimal.mod"){

		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])
		}


	} else if(method =="liu" || method =="liu.mod" ){
		
		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all, pmin[i])
		}

	} else {
		
		stop("Invalid Method:", method)
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




#
#	Function get parameters of optimal test
#
SKAT_META_Optimal_Param<-function(Phi,r.all){


	p.m<-dim(Phi)[2]
	r.n<-length(r.all)

	# ZMZ
	Z.item1.1<- Phi %*% rep(1,p.m)
	ZZ<-Phi
	ZMZ<- Z.item1.1 %*% t(Z.item1.1) / sum(ZZ)

	# W3.2 Term : mixture chisq
	W3.2.t<-ZZ - ZMZ
	lambda<-Get_Lambda(W3.2.t)
	
	# W3.3 Term : variance of remaining ...
	W3.3.item<-sum(ZMZ *(ZZ-ZMZ)) * 4
	
	# tau term 
	z_mean_2<- sum(ZZ)/p.m^2
	tau1<- sum(ZZ %*% ZZ) / p.m^2 / z_mean_2



	# Mixture Parameters
	MuQ<-sum(lambda)
	VarQ<-sum(lambda^2) *2 + W3.3.item
	KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
	Df<-12/KerQ

	# W3.1 Term : tau1 * chisq_1
	tau<-rep(0,r.n)
	for(i in 1:r.n){
		r.corr<-r.all[i]
		term1<-p.m^2*r.corr * z_mean_2 + tau1 * (1-r.corr)
		tau[i]<-term1
	}

	out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau,
z_mean_2=z_mean_2, p.m=p.m,
tau.1 = tau1,
tau.2= p.m*z_mean_2 )

	#param2<<-out
	return(out)
}




#######################################################33
#	Linear and Logistic

SKAT_META_Optimal  = function(Score, Phi, r.all, method="davies", Score.Resampling){

	# if r.all >=0.999 ,then r.all = 0.999
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999	
	}

	p.m<-dim(Phi)[2]
	n.r<-length(r.all)
	

	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-SKAT_META_Optimal_Get_Q(Score, r.all)
	Q.res=NULL
	if(!is.null(Score.Resampling)){
		Q.res<-SKAT_META_Optimal_Get_Q_Res(Score.Resampling, r.all)$Q.r
	}
	Q.all<-rbind(out.Q$Q.r, Q.res) 

	##################################################
	# Compute P-values 
	#################################################

	out<-SKAT_META_Optimal_Get_Pvalue(Q.all, Phi/2, r.all, method)

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
	p.value.resampling= NULL
	if(!is.null(Q.res)){
		p.value.resampling=out$p.value[-1]
	}

 	re<-list(p.value = p.value, param=param, p.value.resampling=p.value.resampling)  
  	
	return(re)	

}

##################################################################
#

Met_SKAT_Get_Pvalue<-function(Score, Phi, r.corr, method, Score.Resampling=NULL){

	#Score.Resampling1<<-Score.Resampling
	p.m<-nrow(Phi)
	Q.res = NULL

	# if Phi==0
	if(sum(abs(Phi)) == 0){
		warning("No polymorphic SNPs!",call.=FALSE)
		return(list(p.value=1, p.value.resampling= NULL, pval.zero.msg=NULL))
	}
	
	if(!is.null(Score.Resampling)){
		Score.Resampling<-t(Score.Resampling)
	}
	if(length(Phi) <=1){
		r.corr=0
	} else{
	
		if(ncol(Phi) <=10){
			if(qr(Phi)$rank <= 1){
				r.corr=0
			}
			
		}
	}

	if(length(r.corr) > 1){
		
		re = SKAT_META_Optimal(Score, Phi, r.corr, method=method, Score.Resampling)
		return(re)
	} 
	
	if (r.corr == 0){
		Q<-sum(Score^2)/2
		
		if(!is.null(Score.Resampling)){
			Q.res<-rowSums(Score.Resampling^2)/2
		}

	} else if (r.corr==1){
		Q=SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
		if(!is.null(Score.Resampling)){
			Q.res<-SKAT_META_Optimal_Get_Q_Res(Score.Resampling, r.corr)$Q.r
		}

		a<- as.matrix(sum(Phi))
		re<-Get_Liu_PVal(Q, a, Q.res)
		return(re)
	} else {

	
		Q=SKAT_META_Optimal_Get_Q(Score, r.corr)$Q.r
		if(!is.null(Score.Resampling)){
			Q.res<-SKAT_META_Optimal_Get_Q_Res(Score.Resampling, r.corr)$Q.r
		}
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Phi<- L %*% (Phi %*% t(L))

		
	}

	re<-Get_Davies_PVal(Q, Phi, Q.res)
	if(length(r.corr)==1){
		re$Q = Q
	}
	return(re)
}

#
#	out_type
#		C: continuous, D:binary, V: Kinship
#
#
SKAT_RunFrom_MetaSKAT<-function(res,Z, X1, kernel, weights=NULL, s2=NULL, pi_1=NULL, P0=NULL, out_type="C", method, res.out, n.Resampling, r.corr){
	
	if (kernel == "linear.weighted") {
    	Z = t(t(Z) * (weights))
  	}
  	
  	# Get Score
  	Score = as.vector(t(Z) %*% res)
  	Score.Resampling=NULL
  	if(!is.null(Score.Resampling)){
  		Score.Resampling = t(res.out) %*% Z
  	}
  	
  	# Phi
  	if(out_type=="C"){
  		Score=Score/ sqrt(s2)
  		if(!is.null(Score.Resampling)){
  			Score.Resampling = Score.Resampling / sqrt(s2)
  		}
  		Phi = t(Z) %*% Z - (t(Z) %*%X1)%*%solve(t(X1)%*%X1)%*% (t(X1) %*% Z )
  	} else if(out_type=="D"){
  		Phi = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*%X1)%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z * pi_1))
  	} else if(out_type=="V"){
  		Phi = t(Z) %*% (P0 %*% Z) # t(Z) P0 Z
  	} else {
  		stop("SKAT_RunFrom_MetaSKAT: no-known out_type!")
  	}
  	
	re = Met_SKAT_Get_Pvalue(Score=Score, Phi=Phi, r.corr=r.corr, method=method, Score.Resampling=Score.Resampling)
	re$IsMeta=TRUE
	return(re)
	
}




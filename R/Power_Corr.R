
#
#	W is beta^2, so make it sqrt(w)
Get_Power_Corr_WR_rho<-function(W,Pi,r.corr){
	
	n<-length(W)
	R<-matrix(rep(1,n*n),ncol=n) * r.corr
	diag(R)<-rep(1,n)

	temp<-sqrt(W) * Pi
	
	W_rho<- t(t(R * temp) * temp)
	diag(W_rho)<-W * Pi

	return(W_rho)

}

Get_Power_Corr.GetPower<-function(K,Mu,alpha.ALL,N.Sample){

	n.a<-length(alpha.ALL)
	OUT.Power<-rep(0,n.a)

	c1<-rep(0,4)
	c2<-rep(0,4)
	c_a<-rep(0,4)

	K1<-K
	K2<-K %*% K
	K3<-K2 %*% K
	
	c1[1] = trace.SKAT(K) * N.Sample
	c1[2] = sum(K *t(K)) * N.Sample^2
	c1[3] = sum(K2 * t(K) ) * N.Sample^3
	c1[4] = sum(K2 * t(K2)) * N.Sample^4

	c2[1] = trace.SKAT(  Mu) * N.Sample^2
	c2[2] = 2 * sum(  K * t(Mu)) * N.Sample^3
	c2[3] = 3 * sum(  K2 * t(Mu)) * N.Sample^4
	c2[4] = 4 * sum(  K3 * t(Mu)) * N.Sample^5

	for(i in 1:4){
		c_a[i]<-c1[i] + c2[i]
	}
		

	for(k in 1:n.a){
		alpha<-alpha.ALL[k]

		out<-Get_Critical_Value(c1, alpha)
		param<-Get_Liu_Params(c_a)

		t = (param$sigmaX/param$sigmaQ)*(out$q.crit-param$muQ) + param$muX
		power<-pchisq(t, df = param$l, ncp = param$d, lower.tail = FALSE )
			
		OUT.Power[k]<-power
	}
	return(OUT.Power)

}

Get_Power_Logistic_R<-function(Z,eta,Beta0,ratio,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25),r.corr){

	Dist<-Dist_Case_Control(eta,Beta0,ratio)
	MAF<-colSums(Z * Dist) /2

	A<-Get_A(Z,eta,Dist,ratio)
	B<-Get_B(Z,eta, Dist,ratio)
	
	W<-Get_W_New(MAF,Weight.Param)		
	
	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL

	for(j in 1:length(N.Sample.ALL)){
			
		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_R<-Get_Power_Corr_WR_rho(W,Pi1,r.corr)

		K<-A %*% W_R
		Mu1<-B %*% W_R
		
		OUT.Power[j,]<-Get_Power_Corr.GetPower(K,Mu1,alpha.ALL,N.Sample)
		
	}

	return(OUT.Power)
}

Get_Power_Continuous_R<-function(Z,eta,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25),r.corr){

	A<-Get_A_Q(Z)
	B<-Get_B_Q(Z,eta)

	MAF<-colMeans(Z ) /2	
	W<-Get_W_New(MAF,Weight.Param)	

	#print(MAF)
	#print(Weight.Param)
	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL
	
	for(j in 1:length(N.Sample.ALL)){
		
		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_R<-Get_Power_Corr_WR_rho(W,Pi1,r.corr)

		#print(sum(Pi1))
		#print(sum(W_R))

		K<-A %*% W_R
		Mu1<-B %*% W_R
		
		OUT.Power[j,]<-Get_Power_Corr.GetPower(K,Mu1,alpha.ALL,N.Sample)
	}
	return(OUT.Power)

}


#
#	p1: causal prop
#	p2: negative prop
Get_Optimal_rho_est<-function(p1,p2){
	rho_est<-p1^2* (p2^2 + (1-p2)^2 - 2* p2*(1-p2))
	return(rho_est)
}


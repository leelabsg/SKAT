SKATPower.env <- new.env()

trace.SKAT = function(x){
  sum(diag(x))
}


Get_Critical_Value<-function(c1,alpha){

	param.null<-Get_Liu_Params_Mod(c1)
	#param.null<-Get_Liu_Params(c1)
	q.crit = (qchisq(1-alpha
	, df = param.null$l, ncp =param.null$d)
	-param.null$muX)*(param.null$sigmaQ/param.null$sigmaX) + param.null$muQ

	return(list(q.crit=q.crit,c.all=c1,param.null=param.null))

}



# MaxValue : Maximum Value at MAF = 0.0001
Get_Beta<-function(Type, MAF, MaxValue,Sign=0){

	n<-length(MAF)
	re<-rep(0,n)
	IDX<-which(MAF > 0)
	if(Type == "Log"){
		re[IDX]<-abs(log10(MAF[IDX])) /4 * MaxValue
	} else if (Type == "Fixed") {
		re[IDX]<-MaxValue
	}	

    	if(Sign > 0){
      		#temp.n<-round(n * Sign)
		temp.n<-floor(n * Sign)
		if(temp.n > 0){
      			temp.idx<-sample(1:n, temp.n)
      			re[temp.idx]<--re[temp.idx]
		}
    	} 
	return(re)
  
}
	
#
Beta_Weights<-function(MAF,Weight.Param){
	

	p<-length(MAF)
	w<-rep(0,p)	
	IDX<-which(MAF > 0)
	
	w[IDX]<-dbeta(MAF[IDX],Weight.Param[1],Weight.Param[2])
	return(w)
}




Get_RequiredSampleSize_OLD<-function(Power.Est, Power=0.8){

	N.Sample.ALL<-as.numeric(rownames(Power.Est))
	alpha<-as.numeric(colnames(Power.Est))

	re<-rep(0,length(alpha))
	for(i in 1:length(alpha)){
		temp1<-smooth.spline(Power.Est[,i],N.Sample.ALL)
		yy1<-predict(temp1,0.8)$y
		#print(yy1)
		re[i]<-ceiling(yy1 /10) *10
	}
	names(re)<-alpha
	return(re)
}


Get_CausalSNPs<-function(MAF, Causal.Ratio, Causal.MAF.Cutoff){

	IDX<-which(MAF < Causal.MAF.Cutoff)
	if(length(IDX) == 0){
		msg<-sprintf("No SNPs with MAF < %f",Causal.MAF.Cutoff)
		stop(msg)
	}
	

	N.causal<-round(Causal.Ratio * length(IDX))
	if(N.causal < 1){
		N.causal = 1
	}
	#print(N.causal)
	#print(Causal.Ratio)
	#print(length(IDX))
	re<-sort(sample(IDX,N.causal))
	return(re)
}

Dist_Case<-function(eta, Beta0){
	
	temp<-exp(Beta0 + eta)/(1+exp(Beta0 + eta))
	temp / sum(temp)
	
}

Dist_Control<-function(eta,Beta0){

	temp<-1/(1+exp(Beta0 + eta))
	temp / sum(temp)

}

Dist_Case_Control<-function(eta,Beta0, ratio_case){

	temp1<-Dist_Case(eta, Beta0) *ratio_case  
	temp2<-Dist_Control(eta, Beta0)*(1-ratio_case)
	
	temp1 + temp2 
}


Get_PI_W_ALL<-function(MAF,N.Sample,Weight.Param){

	W<-Beta_Weights(MAF,Weight.Param)^2
	Pi<- 1 - (1 - MAF)^(2 * N.Sample) 
	W_Pi<-W * Pi

	return(W_Pi)

}


Get_W_New<-function(MAF,Weight.Param){

	W<-Beta_Weights(MAF,Weight.Param)^2
	return(W)

}


Get_PI_New<-function(MAF,N.Sample){

	Pi_1<- 1 - (1 - MAF)^(2 * N.Sample) 
	return(Pi_1)

}

Get_A<-function(Z,eta, Dist,ratio){

	alpha<-log(ratio/(1-ratio))
	Z1<- t(t(Z) - colMeans(Z))
	m<-length(eta)
	V<-exp(eta + alpha) / ( 1 + exp(eta + alpha))^2
	A<- t(Z1) %*% (Z1 * V * Dist)
	
	return(A)

}

Get_B<-function(Z, eta, Dist,ratio){

	m<-length(eta)
	alpha<-log(ratio/(1-ratio))
	mu<-( exp(eta + alpha) / (1 + exp(eta + alpha)) - ratio)
	#mu<- exp(eta) / (1 + exp(eta))
	Z1<- t(t(Z) - colMeans(Z))
	
	B1<- (t(Z1) %*% ( mu * Dist))[,1]
	B<- B1 %*% t(B1)


	return(B)
}


Get_A_Q<-function(Z){

	N<-dim(Z)[1]
	Z1<- t(t(Z) - colMeans(Z))
	A<- t(Z1) %*% (Z1 ) / N
	return(A)

}

Get_B_Q<-function(Z, eta){

	N<-dim(Z)[1]
	m<-length(eta)
	mu<-eta

	Z1<- t(t(Z) - colMeans(Z))
	B1<- (t(Z1) %*%  mu )[,1] / N
	B<- B1 %*% t(B1)
	return(B)
}

Get_Power_Logistic<-function(Z,eta,Beta0,ratio,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25)){

	Dist<-Dist_Case_Control(eta,Beta0,ratio)
	MAF<-colSums(Z * Dist) /2

	A<-Get_A(Z,eta,Dist,ratio)
	B<-Get_B(Z,eta, Dist,ratio)
	
	W<-Get_W_New(MAF,Weight.Param)	
	A1<-t( t(A) * W)
	B1<-t( t(B) * W)

	c1<-rep(0,4)
	c2<-rep(0,4)
	c_a<-rep(0,4)
	re.all<-list()
	idx<-1
	
	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL

	for(j in 1:length(N.Sample.ALL)){
			
		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_Pi<-(W*Pi1)

		A2<-t( t(A1 * Pi1) * Pi1)
		diag(A2)<-diag(A2) / Pi1


		K<-t( t(A1) * Pi1)
		K2<-A1 %*% A2

		Mu1<-t( t(B) * W_Pi)
		Mu2<-B1 %*% A2
		
		#K2<<-K2
		#Pi1<<-Pi1
		#A2<<-A2
	
		c1[1] = trace.SKAT(K) * N.Sample
		c1[2] = trace.SKAT(K2) * N.Sample^2
		c1[3] = sum(K2 * t(K) ) * N.Sample^3
		c1[4] = sum(K2 * t(K2)) * N.Sample^4

		c2[1] = trace.SKAT(  Mu1) * N.Sample^2
		c2[2] = 2 * trace.SKAT( Mu2) * N.Sample^3
		c2[3] = 3 * sum(  Mu2 * t(K)) * N.Sample^4
		c2[4] = 4 * sum(  Mu2 * t(K2)) * N.Sample^5
		
		#print(c1)
		#print(c2)
		#print(diag(A2))

		for(i in 1:4){
			c_a[i]<-c1[i] + c2[i]
		}
		for(k in 1:length(alpha.ALL)){
			alpha<-alpha.ALL[k]
			out<-Get_Critical_Value(c1, alpha)

			param<-Get_Liu_Params(c_a)
			t = (param$sigmaX/param$sigmaQ)*(out$q.crit-param$muQ) + param$muX
			power<-pchisq(t, df = param$l, ncp = param$d, lower.tail = FALSE)

			#re<-list(N.Sample=N.Sample, alpha=alpha, power=power
			#,q.crit=out$q.crit,out=out,param=param,param.null=out$param.null
			#, c1=c1, c2=c2, c_a = c_a)
			#re.all[[idx]]<-re
			
			idx= idx + 1
			OUT.Power[j,k]<-power
		}
	}
	#return(re.all)
	return(OUT.Power)
}

Get_Power_Continuous<-function(Z,eta,alpha.ALL,N.Sample.ALL,Weight.Param=c(1,25)){

	A<-Get_A_Q(Z)
	B<-Get_B_Q(Z,eta)
	MAF<-colMeans(Z ) /2

	W<-Get_W_New(MAF,Weight.Param)		
	A1<-t( t(A) * W)
	B1<-t( t(B) * W)

	c1<-rep(0,4)
	c2<-rep(0,4)
	c_a<-rep(0,4)
	re.all<-list()
	idx<-1

	OUT.Power<-matrix(rep(0,length(alpha.ALL)*length(N.Sample.ALL)),ncol=length(alpha.ALL))
	rownames(OUT.Power)<-N.Sample.ALL
	colnames(OUT.Power)<-alpha.ALL
	
	for(j in 1:length(N.Sample.ALL)){

		N.Sample<-N.Sample.ALL[j]
		Pi1<-Get_PI_New(MAF,N.Sample)
		W_Pi<-(W*Pi1)

		A2<-t( t(A1 * Pi1) * Pi1)
		diag(A2)<-diag(A2) / Pi1


		K<-t( t(A1) * Pi1)
		K2<-A1 %*% A2

		Mu1<-t( t(B) * W_Pi)
		Mu2<-B1 %*% A2
	
		c1[1] = trace.SKAT(K) * N.Sample
		c1[2] = trace.SKAT(K2) * N.Sample^2
		c1[3] = sum(K2 * t(K) ) * N.Sample^3
		c1[4] = sum(K2 * t(K2)) * N.Sample^4

		c2[1] = trace.SKAT(  Mu1) * N.Sample^2
		c2[2] = 2 * trace.SKAT( Mu2) * N.Sample^3
		c2[3] = 3 * sum(  Mu2 * t(K)) * N.Sample^4
		c2[4] = 4 * sum(  Mu2 * t(K2)) * N.Sample^5


		#K1<-K
		#K2<-K %*% K
		#K3<-K2 %*% K
	
		#c1[1] = trace.SKAT(K) * N.Sample
		#c1[2] = sum(K *t(K)) * N.Sample^2
		#c1[3] = sum(K2 * t(K) ) * N.Sample^3
		#c1[4] = sum(K2 * t(K2)) * N.Sample^4

		#c2[1] = trace.SKAT(  Mu1) * N.Sample^2
		#c2[2] = 2 * sum(  K * t(Mu1)) * N.Sample^3
		#c2[3] = 3 * sum(  K2 * t(Mu1)) * N.Sample^4
		#c2[4] = 4 * sum(  K3 * t(Mu1)) * N.Sample^5



		for(i in 1:4){
			c_a[i]<-c1[i] + c2[i]
		}
		for(k in 1:length(alpha.ALL)){
			alpha<-alpha.ALL[k]
			out<-Get_Critical_Value(c1, alpha)

			param<-Get_Liu_Params(c_a)
			t = (param$sigmaX/param$sigmaQ)*(out$q.crit-param$muQ) + param$muX
			power<-pchisq(t, df = param$l, ncp = param$d, lower.tail = FALSE)

			#re<-list(N.Sample=N.Sample, alpha=alpha, power=power
			#,q.crit=out$q.crit,out=out,param=param,param.null=out$param.null
			#, c1=c1, c2=c2, c_a = c_a)

			#re.all[[idx]]<-re
			idx= idx + 1
			
			OUT.Power[j,k]<-power

		}
	}
	#return(re.all)
	return(OUT.Power)

}


Get_RandomRegion<-function(SNP.Dist,SubRegion.Length){
	
	if(SubRegion.Length < 0){
		return(1:length(SNP.Dist))
	}
	total.first<-min(SNP.Dist)
	total.last<-max(SNP.Dist)
	total.length<-total.last -total.first
	
	Region.Start<-runif(1) * (total.length - SubRegion.Length)  + total.first
	Region.End<-Region.Start + SubRegion.Length

	#print(c(Region.Start,Region.End,total.first,total.last))

	Marker.Idx1<-which(SNP.Dist >= Region.Start)
	Marker.Idx2<-which(SNP.Dist <= Region.End)
	IDX<-sort(intersect(Marker.Idx1,Marker.Idx2))
	
	return(IDX)

}

Power_Logistic<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Prevalence=0.01, Case.Prop=0.5, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10), Weight.Param=c(1,25),N.Sim=100, OR.Type = "Log", MaxOR=5, Negative.Percent=0){
	
	
	if(is.null(Haplotypes)){
	
		MSG_SKAT_Example()
		data(SKAT.haplotypes, package="SKAT", envir=SKATPower.env)
		SKAT.haplotypes1<-get("SKAT.haplotypes", envir=SKATPower.env)
		Haplotypes<-SKAT.haplotypes1$Haplotype
		SNP.Location<-SKAT.haplotypes1$SNPInfo$CHROM_POS


	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}

	Marker.MAF.ALL<-colMeans(Haplotypes) 

	# Remove monomorphic 
	id.non<-which(Marker.MAF.ALL==0)
	if(length(id.non) > 0){
		Haplotypes<-Haplotypes[,-id.non]
		SNP.Location<-SNP.Location[-id.non]
		Marker.MAF.ALL<-Marker.MAF.ALL[-id.non]
	}

	# approximated Beta0
	Beta0<-log(Prevalence/(1-Prevalence))
	
	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)	
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:n1,5000, replace=TRUE)
			H2<-sample(1:n1,5000, replace=TRUE)
			#H1[1:n1]<-sample(1:n1,replace=FALSE)
			#H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
	
			
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]

		Causal.Idx<-Get_CausalSNPs(Marker.MAF,  Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		Beta = Get_Beta(OR.Type, Marker.Causal.MAF, log(MaxOR),Negative.Percent/100)

		Causal.Idx1<-IDX.Marker[Causal.Idx]
		eta<-(as.matrix(Haplotypes[,Causal.Idx1]) %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		eta1<-eta[H1] + eta[H2]

		#print(Beta)
		#####################################
		#
		#	Power

		OUT<-Get_Power_Logistic(X1,eta1,Beta0, Case.Prop, alpha,N.Sample.ALL,Weight.Param)
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}
		
		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}
	re<-list(Power = OUT.ALL)
	class(re)<-"SKAT_Power"

	return(re)
}


Power_Continuous<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10)
, Weight.Param=c(1,25), N.Sim=100, BetaType = "Log", MaxBeta=1.6, Negative.Percent=0){
	

	if(is.null(Haplotypes)){

		MSG_SKAT_Example()
		data(SKAT.haplotypes, package="SKAT", envir=SKATPower.env)
		SKAT.haplotypes1<-get("SKAT.haplotypes", envir=SKATPower.env)
		Haplotypes<-SKAT.haplotypes1$Haplotype
		SNP.Location<-SKAT.haplotypes1$SNPInfo$CHROM_POS

	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}

	Marker.MAF.ALL<-colMeans(Haplotypes) 
	# Remove monomorphic 
	id.non<-which(Marker.MAF.ALL==0)
	if(length(id.non) > 0){
		Haplotypes<-Haplotypes[,-id.non]
		SNP.Location<-SNP.Location[-id.non]
		Marker.MAF.ALL<-Marker.MAF.ALL[-id.non]
	}
	
	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]
	out.r_2<-rep(0,N.Sim)	

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)
		
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:n1,5000, replace=TRUE)
			H2<-sample(1:n1,5000, replace=TRUE)
			#H1[1:n1]<-sample(1:n1,replace=FALSE)
			#H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
				
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]
	
		Causal.Idx<-Get_CausalSNPs(Marker.MAF, Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		Beta = Get_Beta(BetaType, Marker.Causal.MAF, MaxBeta,Negative.Percent/100)

		Causal.Idx1<-IDX.Marker[Causal.Idx]
		eta<-(as.matrix(Haplotypes[,Causal.Idx1]) %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		eta1<-eta[H1] + eta[H2]
		

		out.r_2[i]<-sum(Beta^2*2*Marker.Causal.MAF*(1-Marker.Causal.MAF))

		
		#print(Causal.Idx)
		#print(Beta)
		#####################################
		#
		#	Power
		
		OUT<-Get_Power_Continuous(X1,eta1,alpha,N.Sample.ALL,Weight.Param)
		
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}

		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}

	r_sq.v<-mean(out.r_2 /(out.r_2 +1))
	re<-list(Power = OUT.ALL, R.sq = r_sq.v)
	class(re)<-"SKAT_Power"

	return(re)
}


MSG_SKAT_Example<-function(){

		cat("SKAT uses a haplotype dataset (200kb, 10,000 haplotypes) in SKAT.haplotypes to compute power. This dataset was generated using the calibrated coalescent model with mimicking LD structure of European ancestry. See SKAT.haplotypes.\n")

}







Power_Logistic_R<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Prevalence=0.01, Case.Prop=0.5, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10), Weight.Param=c(1,25),N.Sim=100, OR.Type = "Log", MaxOR=5, Negative.Percent=0,r.corr=0){
	
	if(r.corr == 2){
		r.corr = Get_Optimal_rho_est(Causal.Percent/100,Negative.Percent/100)
	}
	
	if(r.corr < 0 || r.corr > 1){
		msg<-("Error: r.corr should be either 0<= r.corr <=1 or r.corr=2 (SKAT-O)")
		stop(msg)	
	}
	
	if(is.null(Haplotypes)){
	
		MSG_SKAT_Example()
		data(SKAT.haplotypes, package="SKAT", envir=SKATPower.env)
		SKAT.haplotypes1<-get("SKAT.haplotypes", envir=SKATPower.env)
		Haplotypes<-SKAT.haplotypes1$Haplotype
		SNP.Location<-SKAT.haplotypes1$SNPInfo$CHROM_POS

		

	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}

	Marker.MAF.ALL<-colMeans(Haplotypes) 
	# Remove monomorphic 
	id.non<-which(Marker.MAF.ALL==0)
	if(length(id.non) > 0){
		Haplotypes<-Haplotypes[,-id.non]
		SNP.Location<-SNP.Location[-id.non]
		Marker.MAF.ALL<-Marker.MAF.ALL[-id.non]
	}


	# approximated Beta0
	Beta0<-log(Prevalence/(1-Prevalence))
	
	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)	
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:n1,5000, replace=TRUE)
			H2<-sample(1:n1,5000, replace=TRUE)
			#H1[1:n1]<-sample(1:n1,replace=FALSE)
			#H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
	
			
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]
		#print(Marker.MAF)

		Causal.Idx<-Get_CausalSNPs(Marker.MAF,  Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		#print(length(Causal.Idx))
		#print(Marker.Causal.MAF)
		Beta = Get_Beta(OR.Type, Marker.Causal.MAF, log(MaxOR),Negative.Percent/100)

		#print(Marker.Causal.MAF)
		#print(Beta)

		Causal.Idx1<-IDX.Marker[Causal.Idx]

		# Seunggeun Change
		#eta<-(Haplotypes[,Causal.Idx1] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		#eta1<-eta[H1] + eta[H2]
		eta1<-(as.matrix(X1[,Causal.Idx]) %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]

		#print(Beta)
		#####################################
		#
		#	Power

		OUT<-Get_Power_Logistic_R(X1,eta1,Beta0, Case.Prop, alpha,N.Sample.ALL,Weight.Param,r.corr)
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}
		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}
	re<-list(Power = OUT.ALL, r.corr=r.corr)
	class(re)<-"SKAT_Power"

	return(re)
}


Power_Continuous_R<-function(Haplotypes=NULL, SNP.Location=NULL, SubRegion.Length=-1, Causal.Percent=5, Causal.MAF.Cutoff=0.03, alpha =c(0.01,10^(-3),10^(-6)),N.Sample.ALL = 500 * (1:10)
,Weight.Param=c(1,25),N.Sim=100,BetaType = "Log", MaxBeta=1.6, Negative.Percent=0,r.corr=0){
	
	if(r.corr == 2){
		r.corr = Get_Optimal_rho_est(Causal.Percent/100,Negative.Percent/100)
	}
	
	if(r.corr < 0 || r.corr > 1){
		msg<-("Error: r.corr should be either 0<= r.corr <=1 or r.corr=2 (SKAT-O)")
		stop(msg)	
	}
	
	if(is.null(Haplotypes)){
	
		MSG_SKAT_Example()
		data(SKAT.haplotypes, package="SKAT", envir=SKATPower.env)
		SKAT.haplotypes1<-get("SKAT.haplotypes", envir=SKATPower.env)
		Haplotypes<-SKAT.haplotypes1$Haplotype
		SNP.Location<-SKAT.haplotypes1$SNPInfo$CHROM_POS


	}
	if(is.null(SNP.Location)){
		stop("Error : SNP.Location is NULL")
	}


	Marker.MAF.ALL<-colMeans(Haplotypes) 
	# Remove monomorphic 
	id.non<-which(Marker.MAF.ALL==0)
	if(length(id.non) > 0){
		Haplotypes<-Haplotypes[,-id.non]
		SNP.Location<-SNP.Location[-id.non]
		Marker.MAF.ALL<-Marker.MAF.ALL[-id.non]
	}

	#####################################
	# 	Compute Power
	######################################	
	OUT.ALL<-NULL
	n1<-dim(Haplotypes)[1]
	out.r_2<-rep(0,N.Sim)

	for(i in 1:N.Sim){
		
		IDX.Marker<-Get_RandomRegion(SNP.Location,SubRegion.Length)
		
		if(n1 >= 5000){		
			H1<-sample(1:n1,replace=FALSE)
			H2<-sample(1:n1,replace=FALSE)
		} else {
			H1<-sample(1:n1,5000, replace=TRUE)
			H2<-sample(1:n1,5000, replace=TRUE)
			#H1[1:n1]<-sample(1:n1,replace=FALSE)
			#H2[1:n1]<-sample(1:n1,replace=FALSE)
		}
				
		X1<-Haplotypes[H1,IDX.Marker] + Haplotypes[H2,IDX.Marker]
		Marker.MAF<-Marker.MAF.ALL[IDX.Marker]

		Causal.Idx<-Get_CausalSNPs(Marker.MAF, Causal.Percent/100, Causal.MAF.Cutoff)
		Marker.Causal.MAF<-Marker.MAF[Causal.Idx]
		Beta = Get_Beta(BetaType, Marker.Causal.MAF, MaxBeta,Negative.Percent/100)
		Causal.Idx1<-IDX.Marker[Causal.Idx]

		out.r_2[i]<-sum(Beta^2*2*Marker.Causal.MAF*(1-Marker.Causal.MAF))

		
		# Seunggeun Change
		#eta<-(Haplotypes[,Causal.Idx1] %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]
		#eta1<-eta[H1] + eta[H2]
		eta1<-(as.matrix(X1[,Causal.Idx]) %*% Beta)[,1] - (t(Marker.Causal.MAF *2)  %*% Beta)[1,1]

		#print(Causal.Idx)
		#print(Beta)
		#####################################
		#
		#	Power
		
		OUT<-Get_Power_Continuous_R(X1,eta1,alpha,N.Sample.ALL,Weight.Param,r.corr)
		
		if(i==1){
			OUT.ALL<-OUT/N.Sim
		} else {
			OUT.ALL<-OUT.ALL + OUT/N.Sim
		}

		if(floor(i/10) * 10 == i){
			msg<-sprintf("%d/%d",i,N.Sim)
			print(msg)
		}
	}
	r_sq.v<-mean(out.r_2 /(out.r_2 +1))
	re<-list(Power = OUT.ALL, R.sq = r_sq.v, r.corr=r.corr)
	class(re)<-"SKAT_Power"

	return(re)
}









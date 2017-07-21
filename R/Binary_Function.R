

SKATExactBin.SKATO_GetQParam<-function(Z, res, idx, r.all, p1, prob_k, k, res.out=res.out, N.sim=10000, Is.sim=FALSE){
	
	
	#void ResampleSTAT(double * Z0, double *Z1, double * Z0_C, double * Z1_C, 
    #              double * teststat_Z0, double *teststat_Z1, double *pteststat_Z0_C, double *pteststat_Z1_C,
    #              double * r_corr, int *pn_r, int *pk, int *pm, int * pn,
    #              int * total_k, double * ncase_k, double * p1,
    #              int *buf1, int * buf2, int *buf3, double * teststat_one, /* buffers */
    #              double * Q, int *err)
                
	Z.1<-cbind(Z[idx,])
	Z0_m<-Z.1 * (-p1[idx])
	Z1_m<-Z.1 * (1-p1[idx])
		
	Z0<-as.vector(t(Z0_m))
	Z1<-as.vector(t(Z1_m))
	Z0_C<-rowSums(Z0_m)
	Z1_C<-rowSums(Z1_m)
	teststat_Z0<-colSums(Z0_m)
	teststat_Z1<-colSums(Z1_m)
	teststat_Z0_C<-sum(teststat_Z0)
	teststat_Z1_C<-sum(teststat_Z1)
	
	m<-ncol(Z.1)
	n<-nrow(Z.1)

	prob<-p1[idx]/sum(p1[idx])
	
	prob_k1<-prob_k / sum(prob_k)
	
	N.resample<-rmultinom(1, N.sim, prob_k1)
	N.resample.total<-sum(N.resample)
	

	RE<-.C("ResampleSTAT"
	, as.double(Z0), as.double(Z1), as.double(Z0_C), as.double(Z1_C)
	, as.double(teststat_Z0), as.double(teststat_Z1), as.double(teststat_Z0_C), as.double(teststat_Z1_C)
	, as.double(r.all), as.integer(length(r.all)), as.integer(n), as.integer(m), as.integer(n)
	, as.integer(N.resample), as.integer(0:n), as.double(prob)
	, integer(n), integer(n), integer(n), double(m)
	, double(N.resample.total * length(r.all)), integer(1))

	err<-RE[[22]]
	if (err != 1) {
        stop("Error in ResampleSTAT!")
    }
    Q.sim<-matrix(RE[[21]], byrow=TRUE, ncol=length(r.all))
	
	# Get moments
	p.m<-ncol(Z.1)
	n.r = length(r.all)
	param<-matrix(rep(0, 3* n.r), nrow=n.r)
	#cat("R.all [",r.all,"]\n")
	
	for(i in 1:n.r){
		muQ<-mean(Q.sim[,i])
		varQ<-var(Q.sim[,i])
		df<-SKAT_Get_DF_Sim(Q.sim[,i])
		param[i,]<-c(muQ, varQ, df)
	}
	
	# compute test statistics
	# no /2
	Q.r.res<-NULL
	if(!is.null(res.out)){
		res.out<-cbind(res[idx], rbind(res.out[idx,]))
	} else {
		res.out<-cbind(res[idx])
	}
	
	temp<-t(res.out) %*% Z.1
	n1<-ncol(res.out)
	
	n.r = length(r.all)
	Q.r.res<-matrix(rep(0,n1 *n.r),ncol=n.r)
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q1<-(1-r.corr) * rowSums(temp^2)
		Q2<-r.corr * p.m^2 * rowMeans(temp)^2
		Q.r.res[,i]<-Q1 + Q2
	}
	
	re<-list(param=param, Q.m=Q.r.res, n.Q = nrow(Q.r.res), Q.sim=Q.sim, N.resample=N.resample)
	return(re);
}



SKATExactBin.GetQ<-function(Z, res, idx, res.out=NULL){

	Z.1<-cbind(Z[idx,])
	# Get Test Statistics
	Q1<-sum((t(Z.1) %*% cbind(res[idx]))^2)
	temp1<-NULL
	if(!is.null(res.out)){
		res.out1<-rbind(res.out[idx,]) 
		temp<-t(Z.1) %*% res.out1
		temp1<-colSums(temp^2)
	
	}

	# Get Q	
	# no /2	
	Q<-c(Q1, temp1)
	n.Q<-length(Q)

	re<-list(Q=Q, n.Q=n.Q)
	return(re)
}


SKATExactBin.ComputeProb_Group = function(idx, pi1, n, ncase, type.group=1){


	#cat("test_type:", test_type, "\n")
	k<-length(idx)
	
	# default : adaptively set groups
	ngroup1=10
	if(type.group==1){
		if(k > 400){
			ngroup1=2
		} else if (k > 100){
			ngroup1=5
		} else if (k > 50){
			ngroup1=6
			
		}
	} else {
		if(k > 400){
			ngroup1=2
		} else if(k > 300){
			ngroup1=6
		} else if (k > 100){
			ngroup1=8
		} else if (k > 50){
			ngroup1=10
		}
	}

	
	if(length(idx)==0){
		return(list(pval=1, k=0, is.return=TRUE))
	}	
	
	p1<-pi1[idx]
	p2<-pi1[-idx]
	
	id.temp<-which(p1 >= 1)
	if(length(id.temp) > 0){
		p1[id.temp]<-0.999
	}
	
	
	weight<-NULL
	group<-NULL
	for(i in 1:ngroup1){
		if(i<ngroup1){
			IDX<-intersect(which(p1 >= (i-1)/ngroup1 ), which(p1 < i/ngroup1 ))
		} else {
			IDX<-intersect(which(p1 >= (i-1)/ngroup1 ), which(p1 <= i/ngroup1))
		}
		
		if(length(IDX)> 0){
			p1.temp<-mean(p1[IDX])
			odd.temp<-p1.temp/(1-p1.temp)
			weight<-c(weight, odd.temp)
			group<-c(group, length(IDX))
		}
	}
	odds.p2<-mean(p2)/(1-mean(p2))
	weight<-c(weight, odds.p2) / odds.p2
	group<-c(group, n-k)
	
	
	# compute probability
	prob_k<-rep(0, k+1)
	
	# compute each prob
	# RGetProb(int* k, int* ngroup, int* ncase, int * group, double * weight, double * prob)
	
	#group1<<-group
	#weight1<<-weight
	#prob_k1<<-prob_k
	#idx<<-idx; pi1<<-pi1; n<<-n; ncase<<-ncase
	
	RE<-.C("RGetProb"
	, as.integer(k), as.integer(length(group)), as.integer(ncase), as.integer(group)
	, as.double(weight), as.double(prob_k));
	
	prob_k = RE[[6]]
	#cat("PROB_K", prob_k)
	return(list(prob_k = prob_k, is.return=FALSE))

}

Get_Total_K = function(k){

	n.total.k<-rep(0,k+1)
	for(i in 0:k){
		n.total.k[i+1]<-choose(k,i)
	}
	
	re = list(n.total = sum(n.total.k), n.total.k = n.total.k)
	return(re)	
}


#
#	input : 
#		idx: idx for samples with variant alleles, pi1: estimated pi, n: total # of samples
#		ncase: number of cases, N.Resampling: number of total resampling
#	
#
#
#
#
SKATExactBin.ComputProb_New<-function(idx, pi1, n, ncase, N.Resampling, ExactMax=1000, test_type=1, type.group=1){

	#idx<<-idx; pi1<<-pi1; n<<-n; ncase<<-ncase; N.Resampling<<-N.Resampling; ExactMax<<-ExactMax; test_type<<-test_type
	#stop(1)
	#cat("test_type:", test_type, "\n")
	
	k<-length(idx)	
	p1<-pi1[idx]
		
	obj.prob_k =SKATExactBin.ComputeProb_Group(idx, pi1, n, ncase, type.group=type.group)
	
	#temp1<-system.time({obj.prob_k =SKATExactBin.ComputeProb_Group(idx, pi1, n, ncase, type.group=type.group)})
	#cat(temp1)
	
	
	if(obj.prob_k$is.return){
		return(obj.prob_k)
	}
	prob_k = obj.prob_k$prob_k
	
	# compute probability
	IsExact<-rep(1,k+1)
	
	# compute numbers
	obj.total = Get_Total_K(k)
	n.total.k = obj.total$n.total.k
	Is.ExactP=TRUE
	if(sum(n.total.k) > N.Resampling){
		
		for(i in 0:k){
			if(n.total.k[i+1] > ExactMax){
				n.total.k[i+1]<-ceiling(N.Resampling * prob_k[i+1])
				IsExact[i+1] = 0;
			}
		}
		
		Is.ExactP=FALSE
	}
	n.total<-sum(n.total.k)
	
	re<-list(prob_k=prob_k, k=k, n=n, n.total=n.total, n.total.k=n.total.k, IsExact=IsExact, p1=p1, Is.ExactP=Is.ExactP)
	return(re)


}

#
#	input : 
#		idx: idx for samples with variant alleles, pi1: estimated pi, n: total # of samples
#		ncase: number of cases, N.Resampling: number of total resampling
#	
#
#
#
#
SKATExactBin.ComputProb_Random<-function(obj.prob_k, idx, pi1, n, ncase, N.Resampling, ExactMax=1000, test_type=1){

	#idx<<-idx; pi1<<-pi1; n<<-n; ncase<<-ncase; N.Resampling<<-N.Resampling; ExactMax<<-ExactMax; test_type<<-test_type
	#stop(1)
	#cat("test_type:", test_type, "\n")
	
	k<-length(idx)
	p1<-pi1[idx]
		
	if(obj.prob_k$is.return){
		return(obj.prob_k)
	}
	prob_k = obj.prob_k$prob_k
	
	# compute probability
	n.total.s<-rmultinom(1, N.Resampling, prob=prob_k)
	IsExact<-rep(0,k+1)
	
	# compute numbers
	obj.total = Get_Total_K(k)
	n.total.k = obj.total$n.total.k
	
	if(sum(n.total.k) > N.Resampling){
		
		for(i in 0:k){
		
			if(n.total.s[i+1] > 0){
				if(n.total.s[i+1] <= ExactMax && n.total.k <= ExactMax){
					n.total.s[i+1] = n.total.k[i+1]
					IsExact[i+1] = 1;
				} else if( n.total.s[i+1] <= ExactMax && n.total.k > ExactMax){
					n.total.s[i+1] = ExactMax
				}
			}
		}
	}
	n.total<-sum(n.total.s)
	
	re<-list(prob_k=prob_k, k=k, n=n, n.total=n.total, n.total.k=n.total.s, IsExact=IsExact, p1=p1)
	return(re)


}


Get_Res_Arrays<-function(res, res.out, idx){

	if(!is.null(res.out)){
		res.out<-cbind(res[idx], rbind(res.out[idx,]))
	} else {
		res.out<-cbind(res[idx])
	}
	
	resarray<-NULL
	nres<-ncol(res.out)
	nres_k<-rep(0, nres)
	for(i in 1:nres){
		id1<-sort(which(res.out[,i] > 0)) -1 	# array should start from 0
		nres_k[i]<-length(id1)
		resarray<-c(resarray, id1)
	}
	
	re<-list(resarray=resarray, nres=nres, nres_k = nres_k)
	return(re)
	
}


#
#	Get p-values using res.out
#
Get_Resample_P<-function(Z, res, res.out){

	Q<-sum((t(Z) %*% res)^2)
	temp<-t(Z) %*% res.out
	temp1<-colSums(temp^2)
	
	n1<-length(temp1)
	idx1<-which(Q == temp1)
	idx2<-which(Q <= temp1)
	
	return(c(length(idx2)/n1, length(idx1)/n1))
	
}



#
#	Check genotype matrix
#

SKATBinary.Single.CheckZ<-function(Z, id_include, impute.method, is_check_genotype, is_dosage, estimate_MAF=1){

	#############################################
	# Check parameters

	if(is.matrix(Z)){
		if(ncol(Z) > 1){
			stop("Z should have one column")
		}
		Z = Z[,1]
	} else if(!is.vector(Z)){

		msg<-sprintf("Z should be either a vector or a matrix with one column")
		stop(msg)
	}

 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}


	#####################################################
	# Check Z

	if(!is_check_genotype && !is_dosage){
		Z.test<-Z[id_include]
		return(list(Z.test=Z.test, MAF=0, return=0) )
	}

	if(estimate_MAF==2){
		Z<-Z[id_include]
		id_include<-1:length(id_include)
	}

	##############################################
	# Check Missing 

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	##################################################################
	# doing imputation

	MAF<-mean(Z, na.rm = TRUE)/2
	MAF1<-mean(Z[id_include],na.rm=TRUE)/2
	MAF_Org=MAF
	
	##########################################
	# Missing Imputation
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
	
		warning(msg,call.=FALSE)
		Z<-Impute(cbind(Z),impute.method)[,1]
	} 

	#########################################
	# Check and recal
	
	MAF<-mean(Z, na.rm = TRUE)/2
	MAF1<-mean(Z[id_include],na.rm=TRUE)/2
	
	###########################################
	# Check non-polymorphic
	if(MAF1==0){
		
		msg<-sprintf("No polymorphic SNP. P-value = 1" )
		warning(msg,call.=FALSE)
		re<-list(p.value = 1, MAP=1, K=0, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   
		return(re)
	} 
	
	Z.test=Z[id_include]
	if(MAF1 > 0.5){
		Z.test = 2-Z.test
	}
	
	
	MAC = sum(Z.test)
	idx<-which(Z.test > 0)
	m<-length(idx)

	return(list(Z.test=Z.test, MAF=MAF, MAC=MAC, m=m, return=0))

}


#
#	Check genotype matrix
#
SKATExactBin_Check<-function(Z, obj, kernel = "linear.weighted", weights.beta=c(1,25), weights = NULL, impute.method = "bestguess"
, r.corr=0, is_dosage = FALSE, missing_cutoff, max_maf, SetID, estimate_MAF, Is.Single=FALSE, Is.MakeZ1 = TRUE){

	is_check_genotype=TRUE
	obj.res = SKATExactBin_CheckObj(obj)
	y = round(obj.res$mu + obj.res$res)
	pi1 = obj.res$mu;
	ncase = sum(y)
	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	
	SKAT_Check_RCorr(kernel, r.corr)
	out.z<-SKAT_MAIN_Check_Z(Z, n, obj.res$id_include, SetID, weights, weights.beta, impute.method
	, is_check_genotype, is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
	
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		out.z$param$n.marker.test =0
		return(out.z)
	}
	
	##############################################
	#
	#	Calculate MAC and m
	#


	Z1 = cbind(out.z$Z.test)
	#Z2<<-Z1
	#Z3<<-Z
	#out.z<<-out.z
	#id_include<<-obj.res$id_include
	MAC= sum(Z1)
 	rsum<-rowSums(Z1)
	idx<-which(rsum > 0)
	m<-length(idx)
  	
  	if(!Is.MakeZ1){
  		re=list(Z1=Z1, obj.res=obj.res, pi1=pi1, ncase=ncase, r.corr=r.corr, idx=idx, y=y, MAC=MAC, m=m, return=0)
		return(re)
  	}
  	
	##############################################
	#
	#	Get Z1
	#
		
	if(Is.Single){
	
		re=list(Z1=Z1, obj.res=obj.res, pi1=pi1, ncase=ncase, r.corr=0, idx=idx, y=y, MAC=MAC, m=m, return=0)
		return(re)
	
	}

	#################################
	# 	if there is only one variant or rank of Z1=1
  	  	
	if(length(r.corr) > 1 && ncol(Z1) <= 1){
		r.corr=0
	} else if(length(r.corr) > 1 && sum(abs(Z1 - Z1[,1])) == 0){
		r.corr=0
		
		msg<-sprintf("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!" )
		warning(msg,call.=FALSE)
	}
	
	#################################
	# 	weighting and etc
	
  	if (kernel == "linear.weighted") {
    	Z1 = t(t(Z1) * (out.z$weights))
  	}
	
	n.marker.test = ncol(Z1)
 	if(length(r.corr) == 1){
  		
  		if(r.corr == 1){
  			Z1<-cbind(rowSums(Z1))
  		} else if(r.corr > 0){

   			p.m<-dim(Z1)[2]	
			R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
			L<-chol(R.M,pivot=TRUE)
			Z1<- Z1 %*% t(L) 
  		}
  	}

	re=list(Z1=Z1, obj.res=obj.res, pi1=pi1, ncase=ncase, r.corr=r.corr, idx=idx, y=y, MAC=MAC, m=m, n.marker.test=n.marker.test, return=0)
	return(re)

}



#
#
#	Z: weighted 
#

SKATExactBin_CheckAccuracy<-function(pval, total, r.critical=25 ){
	
	# if total
	
	if(total * pval >= r.critical){
			return(TRUE)
	}
	return (FALSE)

}


SKATExactBin_CheckObj<-function(obj.res){

	if(class(obj.res) == "SKAT_NULL_Model_ADJ"){
		obj.res = obj.res$re1
	} else if(class(obj.res) == "SKAT_NULL_Model"){
		if(obj.res$out_type !="D"){
			stop("out_type in SKAT_Null_Model should be D!")
		}
	} else {
		stop("Error: obj.res!")
	}
	
	return(obj.res)

}


SKATExactBin.Adaptive<-function(Z, obj, kernel = "linear.weighted", weights.beta=c(1,25), weights = NULL, impute.method = "bestguess"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE 
, N.Iter = 50000, N.Resampling=2 * 10^6, missing_cutoff=0.15, max_maf=1, SetID = NULL, estimate_MAF=1, epsilon){


	obj = SKATExactBin_Check(Z=Z, obj=obj, kernel = kernel, weights.beta=weights.beta
	, weights = weights, impute.method = impute.method, r.corr=r.corr, is_dosage = is_dosage, estimate_MAF=estimate_MAF, 
	missing_cutoff=missing_cutoff, max_maf=max_maf, SetID = SetID)
	
	if(obj$return ==1){
		return(obj)
	}
		
	res.out=obj$obj.res$res.out
  	#p.value.resampling=NULL
  	#n.resample.test=0
  	#if(!is.null(res.out)){
  	#	n.resample.test<-ncol(res.out) 
  	#	p.value.resampling=rep(NA, n.resample.test)
  	#}
  	
  	
  	re = SKATExactBin.Work_Adaptive(obj$Z1, obj$obj.res$res, obj$pi1, obj$ncase, obj$idx, obj$r.corr, res.out=res.out, 
  	Is.testdata=FALSE, File=NULL, N.Iter = N.Iter, N.Resampling=N.Resampling, test_type=3, epsilon=epsilon)

	re$is.accurate =TRUE
	if(!SKATExactBin_CheckAccuracy(re$p.value, re$n.total, r.critical=25 )){
  		re$is.accurate = FALSE
  	}
  	
  	#if(n.resample.test > 0){
  	#	for(i in 1:n.resample.test){
  	#	
  	#		re1 = SKATExactBin.Work_Adaptive(obj$Z1, res.out[,i], obj$pi1, obj$ncase, obj$idx, obj$r.corr, 
  	#Is.testdata=FALSE, File=NULL, N.Iter = N.Iter, N.Resampling=N.Resampling, test_type=3, epsilon=epsilon)
	#		p.value.resampling[i] = re1$p.value
	#	}
  	#}
  			
	param = list(n.marker=ncol(Z), n.marker.test=obj$n.marker.test)

	re$param=param
	re$m = length(obj$idx)
	re$MAC= obj$MAC
	re$MAP = -1
	

	return(re)


}


SKATExactBin<-function(Z, obj, kernel = "linear.weighted", method.bin="ER"
, weights.beta=c(1,25), weights = NULL, impute.method = "bestguess"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, SetID = NULL, estimate_MAF=1,
N.Resampling=10^7, ExactMax=10000, test_type=1, Is.testdata=FALSE, File=NULL, Is.Single=FALSE, epsilon){

	
	obj = SKATExactBin_Check(Z=Z, obj=obj, kernel = kernel, weights.beta=weights.beta
	, weights = weights, impute.method = impute.method, r.corr=r.corr, is_dosage = is_dosage, max_maf=max_maf, estimate_MAF=estimate_MAF,
	missing_cutoff=missing_cutoff, SetID = SetID, Is.Single=Is.Single)

	
	if(obj$return ==1){
		return(obj)
	}

	
	# Run ER.R if n* > 100
	if(method.bin=="ER" && length(obj$idx) >= 100){
		method.bin="ER.R"
	}
  	
  	if(method.bin=="QA" || method.bin=="MA"){
	
		re=SKATExactBin.Moment(obj$Z1, obj$obj.res$res, obj$pi1, obj$ncase, obj$idx, obj$r.corr, 
			res.out=obj$obj.res$res.out, X1=obj$obj.res$X1, pi_1=obj$obj.res$pi_1, method.bin=method.bin)
		re$MAP=-1
		
	} else if(method.bin=="ER.R"){
		re=SKATExactBin.ER_R(obj$Z1, obj$obj.res$res, obj$pi1, obj$ncase, obj$idx, obj$r.corr, 
			res.out=obj$obj.res$res.out, X1=obj$obj.res$X1, pi_1=obj$obj.res$pi_1, N.Resampling=N.Resampling)

  			re$is.accurate =TRUE
  			if(!SKATExactBin_CheckAccuracy(re$p.value, re$n.total, r.critical=25 )){
  				re$is.accurate = FALSE
  			}
  				
	} else if(method.bin=="ER") {

  			re = SKATExactBin.Work(obj$Z1, obj$obj.res$res, obj$pi1, obj$ncase, obj$idx, obj$r.corr, res.out=obj$obj.res$res.out, 
  			Is.testdata=Is.testdata, File=File, N.Resampling=N.Resampling, ExactMax=ExactMax, epsilon=epsilon, test_type=test_type)
  			
  			re$is.accurate =TRUE
  			if(!SKATExactBin_CheckAccuracy(re$p.value, re$n.total, r.critical=25 ) && !re$Is.ExactP){
  				re$is.accurate = FALSE
  			}
  	
	} else {
		msg = sprintf("Error! %s is not correct method type !", method.bin)
		stop(msg)
	}
	
	if(is.null(re$param)){
		param = list(n.marker=ncol(Z), n.marker.test=obj$n.marker.test)
		re$param=param
	} else {
		re$param$n.marker=ncol(Z)
		re$param$n.marker.test=obj$n.marker.test
	}
	re$m = length(obj$idx)
	re$MAC= obj$MAC
		
	return(re)
}


SKATExactBin.Firth<-function(Z, obj, kernel = "linear.weighted", method.bin="ER"
, weights.beta=c(1,25), weights = NULL, impute.method = "bestguess"
, r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, SetID = NULL, estimate_MAF=1, Is.Single=FALSE){

	if(ncol(Z) != 1 && r.corr !=1){
		stop("Firth method only can be used for single variant or burden test!")
	} 
	
	obj = SKATExactBin_Check(Z=Z, obj=obj, kernel = kernel, weights.beta=weights.beta
	, weights = weights, impute.method = impute.method, r.corr=r.corr, is_dosage = is_dosage, estimate_MAF=estimate_MAF,
	missing_cutoff=missing_cutoff, max_maf=max_maf, SetID = SetID, Is.Single=Is.Single)
	
	if(obj$return ==1){
		return(obj)
	}
	
	obj.res = obj$obj.res
	z1 = obj$Z1
	if(ncol(z1) > 1){
		stop("Firth method only can be used for single variant or burden test!")
	}
	
	k = ncol(obj.res$X1)
	COV.all = cbind(z1, obj.res$X1) 
	
	#COV.all1<<-COV.all;y1<<-obj$y;k<<-k
	fit.full <- fast.logistf.fit(x = COV.all, y=obj$y)
	fit.i <- fast.logistf.fit(x = COV.all, col.fit = (1:(k+1))[-1], y=obj$y)
    pval = pchisq(2 * (fit.full$loglik - fit.i$loglik), 1, lower.tail = FALSE)
    
    
    p.value.resampling=NULL
    prob.tie.resampling=NULL
	if(obj.res$n.Resampling > 0){
		p.value.resampling = rep(NA, obj.res$n.Resampling)
		prob.tie.resampling= rep(0, obj.res$n.Resampling)
		y1 = round(obj.res$mu + obj.res$res.out)
		
		for(i in 1:obj.res$n.Resampling){
			fit.full <- fast.logistf.fit(x = COV.all, y=y1[,i])
			fit.i <- fast.logistf.fit(x = COV.all, col.fit = (1:(k+1))[-1], y=y1[,i])
    		p.value.resampling[i] = pchisq(2 * (fit.full$loglik - fit.i$loglik), 1, lower.tail = FALSE)
		}
	}
	
    re=list(p.value=pval, MAP=0, p.value.resampling=NULL, n.total=NULL)
    re$m = length(obj$idx)
	re$MAC= obj$MAC
	
    return(re)
    
}

SKATExactBin.Moment<-function(Z, res, pi1, ncase, idx, r.corr, res.out=NULL, X1, pi_1, method.bin){

	if(method.bin=="QA"){
		N.sim = 500000
	} else {
		#method.bin="ECP"
		N.sim = 10000
	}
	
	n<-length(res)
	p1<-pi1[idx]
	p2<-pi1[-idx] 
	
	Z.1<-cbind(Z[idx,])
	Z0<-as.vector(t(Z.1 * (-p1)))
	Z1<-as.vector(t(Z.1 * (1-p1)))
	test_Z0<-colSums(cbind(Z.1 * (-p1)))
	m<-ncol(Z.1)
	

	n.Resampling<-0
	if(!is.null(res.out)){
		n.Resampling<-ncol(res.out)
	}
	pr<-SKATExactBin.ComputProb_New(idx, pi1, n, ncase, N.Resampling=10^7, ExactMax=1, test_type=1, type.group=1)
	
	if(length(r.corr)==1){
		
		re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, 0, pi1, pr$prob_k, pr$k, res.out=res.out, N.sim=N.sim, Is.sim=TRUE)
		#re.Q1<<-re.Q
		out<-KMTest.logistic.Linear.VarMatching(res, Z, X1=X1, kernel="linear"
		, weights = NULL,pi_1=pi_1 , method=method.bin, res.out,n.Resampling, r.corr=0, mu=pi1, res.moments = NULL, Q.sim=re.Q$Q.sim[,1]/2)
	} else {

		Z.item2<-SKAT_Optimal_GetZItem2(Z)

		re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, r.corr, pi1, pr$prob_k, pr$k, res.out=res.out, N.sim=N.sim, Is.sim=TRUE)
		re.Q2<-SKATExactBin.SKATO_GetQParam(cbind(Z.item2), res, idx, 0, pi1, pr$prob_k, pr$k, res.out=NULL, N.sim=N.sim, Is.sim=TRUE)
	

		out<-SKAT_Optimal_Logistic_VarMatching(res, Z, X1, kernel="linear", weights = NULL, pi_1=pi_1, method = method.bin
		, res.out=res.out, n.Resampling =n.Resampling, r.all=r.corr, mu=pi1, res.moments = NULL, Q.sim=re.Q2$Q.sim[,1]/2, Q.sim.a=re.Q$Q.sim/2)
	
	}
	
	return(out)


}


SKATExactBin.ER_R<-function(Z, res, pi1, ncase, idx, r.corr, res.out=NULL, X1, pi_1, N.Resampling=N.Resampling){

	N.sim = N.Resampling
	n<-length(res)
	p1<-pi1[idx]
	p2<-pi1[-idx] 
	
	Z.1<-cbind(Z[idx,])
	Z0<-as.vector(t(Z.1 * (-p1)))
	Z1<-as.vector(t(Z.1 * (1-p1)))
	test_Z0<-colSums(cbind(Z.1 * (-p1)))
	m<-ncol(Z.1)
	

	n.Resampling<-0
	if(!is.null(res.out)){
		n.Resampling<-ncol(res.out)
	}
	pr<-SKATExactBin.ComputProb_New(idx, pi1, n, ncase, N.Resampling=10^7, ExactMax=1, test_type=1, type.group=1)
	
	
	if(length(r.corr)==1){
		
		re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, 0, pi1, pr$prob_k, pr$k, res.out=res.out, N.sim=N.sim, Is.sim=TRUE)
		#re.Q1<<-re.Q
		
		# calculate p-values
		pval.a<-rep(0,re.Q$n.Q)
		for(i in 1:re.Q$n.Q){
			pval.a[i]<-(length(which(re.Q$Q.sim >= re.Q$Q.m[i])) +1)/(N.sim+1)
		}	
		
		
	} else {

		Z.item2<-SKAT_Optimal_GetZItem2(Z)

		re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, r.corr, pi1, pr$prob_k, pr$k, res.out=res.out, N.sim=N.sim, Is.sim=TRUE)
		#re.Q<<-re.Q
		param<-re.Q$param
		Q.m<-(t(re.Q$Q.m)  - param[,1]) /sqrt(param[,2]) * sqrt(2 * param[,3]) + param[,3]
		Q.sim1<-(t(re.Q$Q.sim)  - param[,1]) /sqrt(param[,2]) * sqrt(2 * param[,3]) + param[,3]
		
		Test.stat<-apply(Q.m, 2, max)
		Test.stat.Sim<-apply(Q.sim1, 2, max)

		pval.a<-rep(0,re.Q$n.Q)
		for(i in 1:re.Q$n.Q){
			pval.a[i]<-(length(which(Test.stat.Sim >= Test.stat[i])) +1)/(N.sim+1)
		}	

		
		#re.Q2<-SKATExactBin.SKATO_GetQParam(cbind(Z.item2), res, idx, 0, pi1, pr$prob_k, pr$k, res.out=NULL, N.sim=N.sim, Is.sim=TRUE)
		#re.Q2<<-re.Q2

		#out<-SKAT_Optimal_Logistic_VarMatching(res, Z, X1, kernel="linear", weights = NULL, pi_1=pi_1, method = method.bin
		#, res.out=res.out, n.Resampling =n.Resampling, r.all=r.corr, mu=pi1, res.moments = NULL, Q.sim=re.Q2$Q.sim[,1]/2, Q.sim.a=re.Q$Q.sim/2)
	
	}
	
	p.value.resampling=NULL
	if(n.Resampling> 0){
		p.value.resampling=pval.a[-1]
	}
	p.value=pval.a[1]
	
	re=list(p.value=p.value, MAP=-1, p.value.resampling=p.value.resampling, n.total=N.sim)
	return(re)

}

SKATExactBin.Work<-function(Z, res, pi1, ncase, idx, r.corr, res.out=NULL, Is.testdata=FALSE, File=NULL, N.Resampling=2 *10^6
, ExactMax,epsilon, test_type=1){

	n<-length(res)
	p1<-pi1[idx]
	p2<-pi1[-idx] 
	
	Z.1<-cbind(Z[idx,])
	Z0<-as.vector(t(Z.1 * (-p1)))
	Z1<-as.vector(t(Z.1 * (1-p1)))
	test_Z0<-colSums(cbind(Z.1 * (-p1)))
	m<-ncol(Z.1)
	
	pr<-SKATExactBin.ComputProb_New(idx, pi1, n, ncase, N.Resampling, ExactMax, test_type, type.group=2)
	p1_adj<-as.double(pr$p1/mean(pr$p1))
	odds<-pr$p1/(1-pr$p1)
	
	#
	#	test_type=1 : org
	#	test_type=3 : use resampling_random when count < 5000
	#	
	test_type=1

	re.arr<-Get_Res_Arrays(res, res.out, idx)
	n.Q<-re.arr$nres
	pval<-rep(0, n.Q)
	pval1<-rep(0,n.Q)
	minP<-100
	
	#odds1<<-odds
	if(length(r.corr) == 1){
			
		if(Is.testdata){
			Generate_TestData(File, re.arr$resarray, re.arr$nres, re.arr$nres_k, Z0, Z1, test_Z0, pr$k, 
			m, pr$n.total, pr$n.total.k, pr$prob_k, odds, p1_adj, pr$IsExact);
			return(1);
		}
	
		RE<-.C("RSKATExact"
		, as.integer(re.arr$resarray), as.integer(re.arr$nres), as.integer(re.arr$nres_k), as.double(Z0), as.double(Z1)
		, as.integer(pr$k), as.integer(m), as.integer(pr$n.total), as.integer(pr$n.total.k), as.double(pr$prob_k)
		, odds, p1_adj, as.integer(pr$IsExact), as.double(pval) , as.double(pval1), as.double(minP),as.integer(test_type), as.double(epsilon))


		pval.re<-cbind(RE[[14]], RE[[15]])
		minP<-RE[[16]][1] /2 	# MidP-value MAP
		prob_k_out<-RE[[10]]

		pval1<-cbind(pval.re[,1], pval.re[,1] - pval.re[,2] /2, pval.re[,2])
		
	} else {
		re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, r.corr, pi1, pr$prob_k, pr$k, res.out=res.out)
		param = as.vector(t(re.Q$param))
	
		RE<-.C("RSKATOExact"
			, as.integer(re.arr$resarray), as.integer(re.arr$nres), as.integer(re.arr$nres_k), as.double(Z0), as.double(Z1)
			, as.double(r.corr), as.integer(length(r.corr)), as.double(param)
			, as.integer(pr$k), as.integer(m), as.integer(pr$n.total), as.integer(pr$n.total.k), as.double(pr$prob_k)
			, odds, p1_adj, as.integer(pr$IsExact), as.double(pval) , as.double(pval1), as.double(minP), as.integer(test_type), as.double(epsilon))

	
		pval.re<-cbind(RE[[17]], RE[[18]])
		minP<-RE[[19]][1] /2 	# MidP-value MAP
		prob_k_out<-RE[[10]]
		
		pval1<-cbind(pval.re[,1], pval.re[,1] - pval.re[,2] /2, pval.re[,2])

	} 
	if(!pr$Is.ExactP){
		minP=-1
	}

	p.value.resampling=NULL
	p.value.standard.resampling=NULL
	Q.resampling=NULL
	if(n.Q > 1){
		p.value.resampling=pval1[-1,2]
		p.value.standard.resampling=pval1[-1,1]
	}

	
	re=list(p.value=pval1[1,2], p.value.standard=pval1[1,1], MAP=minP
	, p.value.resampling=p.value.resampling, p.value.standard.resampling=p.value.standard.resampling
	,m=pr$k, n.total=pr$n.total, Is.ExactP=pr$Is.ExactP)
	#, n.total.k=pr$n.total.k, prob_k = prob_k_out)
	
	return(re)

}





SKATExactBin.Work_Adaptive<-function(Z, res, pi1, ncase, idx, r.corr, res.out=NULL,  Is.testdata=FALSE, File=NULL, N.Iter = 50000
, N.Resampling=2 * 10^6, epsilon, test_type=3){

	n<-length(res)
	p1<-pi1[idx]
	p2<-pi1[-idx] 
	
	Z.1<-cbind(Z[idx,])
	Z0<-as.vector(t(Z.1 * (-p1)))
	Z1<-as.vector(t(Z.1 * (1-p1)))
	test_Z0<-colSums(cbind(Z.1 * (-p1)))
	m<-ncol(Z.1)
	
	#pr<-SKATExactBin.ComputProb_New(idx, pi1, n, ncase, N.Resampling, ExactMax, test_type)
	
	total.iter = ceiling(N.Resampling / N.Iter)
	
	r.critical=25
	
	#
	#	test_type=1 : org
	#	test_type=3 : use resampling_random when count < 5000
	#	
	test_type=3
	obj.prob_k =SKATExactBin.ComputeProb_Group(idx, pi1, n, ncase)	

  	n.resample.test=0
	if(!is.null(res.out)){
		n.resample.test<-ncol(res.out) 
	}
	
	pval.all<-matrix(rep(NA, 3*(1+n.resample.test)), ncol=3)
	
	for(k in 1:(n.resample.test+1)){
	
		pval.a<-NULL
		if(k==1){
			re.arr<-Get_Res_Arrays(res, NULL, idx)
		} else {
			re.arr<-Get_Res_Arrays(res.out[,k-1], NULL, idx)
		}
		
		for(l in 1:total.iter){ 

			pr<-SKATExactBin.ComputProb_Random(obj.prob_k, idx, pi1, n, ncase, N.Iter, ExactMax=0, test_type=1)
	
			p1_adj<-as.double(pr$p1/mean(pr$p1))
			odds<-pr$p1/(1-pr$p1)
			
			n.Q<-re.arr$nres
			pval<-rep(0, n.Q)
			pval1<-rep(0,n.Q)
			minP<-100
	
			if(length(r.corr) ==1){
			
	
				RE<-.C("RSKATExact"
				, as.integer(re.arr$resarray), as.integer(re.arr$nres), as.integer(re.arr$nres_k), as.double(Z0), as.double(Z1)
				, as.integer(pr$k), as.integer(m), as.integer(pr$n.total), as.integer(pr$n.total.k), as.double(pr$prob_k)
				, odds, p1_adj, as.integer(pr$IsExact), as.double(pval) , as.double(pval1), as.double(minP),as.integer(test_type), as.double(epsilon))

				pval.re<-cbind(RE[[14]], RE[[15]])
				prob_k_out<-RE[[10]]

				pval1<-cbind(pval.re[,1], pval.re[,1] - pval.re[,2] /2, pval.re[,2])
		
			} else {
				if(l==1){
					re.Q<-SKATExactBin.SKATO_GetQParam(Z, res, idx, r.corr, pi1, pr$prob_k, pr$k, res.out=res.out)
					param = as.vector(t(re.Q$param))
				}
	
				RE<-.C("RSKATOExact"
					, as.integer(re.arr$resarray), as.integer(re.arr$nres), as.integer(re.arr$nres_k), as.double(Z0), as.double(Z1)
					, as.double(r.corr), as.integer(length(r.corr)), as.double(param)
					, as.integer(pr$k), as.integer(m), as.integer(pr$n.total), as.integer(pr$n.total.k), as.double(pr$prob_k)
					, odds, p1_adj, as.integer(pr$IsExact), as.double(pval) , as.double(pval1), as.double(minP), as.integer(test_type), as.double(epsilon))

	
				pval.re<-cbind(RE[[17]], RE[[18]])
				prob_k_out<-RE[[10]]
			
				pval1<-cbind(pval.re[,1], pval.re[,1] - pval.re[,2] /2, pval.re[,2])

			}	 

			if(l==1){
				pval.a = pval1
			} else {
				pval.a = pval.a *(l-1)/l + pval1/l
			}
		
			n.Total = N.Iter * l
			if(SKATExactBin_CheckAccuracy(pval.a[1], n.Total, r.critical=25 )){
				break
			}
		}
		pval.all[k,] = pval.a
		
	}
	p.value.resampling=NULL
	p.value.standard.resampling=NULL
	
	if(	n.resample.test > 0){
	
		p.value.resampling=pval.all[-1,2]
		p.value.standard.resampling=pval.all[-1,1]
	
	}
	
	re=list(p.value=pval.all[1,2], p.value.standard=pval.all[1,1], MAP=-1, p.value.resampling=p.value.resampling, p.value.standard.resampling=p.value.standard.resampling
	, m=pr$k, n.total=n.Total)
	return(re)

}



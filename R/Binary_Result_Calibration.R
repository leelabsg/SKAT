

SKATExactBin_SingleMAP<-function(z, obj.res){

	# Z should be a vector
	if(!is.vector(z)){
		stop("z should be a vector!")
	}


	if(class(obj.res) == "SKAT_NULL_Model_ADJ"){
		obj.res = obj.res$re1
	} else if(class(obj.res) == "SKAT_NULL_Model"){
		if(obj.res$out_type !="D"){
			stop("out_type in SKAT_Null_Model should be D!")
		}
	} else {
		stop("Error: obj.res!")
	}


	pi1 = obj.res$mu;
	if(mean(z,  na.rm = TRUE) <= 0.5){
		idx<-which(z > 0)
	} else {
		idx<-which(z < 2)
	}
	
	re<-prod(pi1[idx])
	return(re)
	
}



Get_EffectiveNumberTest_Item<-function(ntest, MAP, alpha){

	alpha1<-alpha/ntest
	n1<-length(which(MAP < alpha1))
	re<-n1 - ntest

	return(re)
}


Get_EffectiveNumberTest_Item_MidP<-function(ntest, MAP, alpha){

	#ntest<-n
	alpha1<-alpha/ntest
	
	n1<-length(which(MAP*2 < alpha1))
	idx<-intersect(which(MAP >= alpha1), which(MAP*2 < alpha1))
	
	#k<-n1 * alpha/(alpha - sum(MAP[idx])*2)
	b1<-(alpha - 2 *sum(MAP[idx]))/(alpha * n1)
	k<-1/b1
	re=ntest-k
	return(re)
}



#
#	
# Treat MAP=NA as MAP=0
#
Get_EffectiveNumberTest<-function(MAP, alpha=0.05, Is.MidP=TRUE){
	
	#alpha=0.05
	n<-length(MAP)
	#n1<-ceiling((1:k)/k *n)
	#alpha<-0.05 /n1
	
	id.miss<-which(is.na(MAP))
	if(length(id.miss) > 0){
		MAP[id.miss] = 0
	}	
	
	if(length(which(MAP < alpha)) == 0){
		return(1)
	}

	if(!Is.MidP){
		re<-uniroot(Get_EffectiveNumberTest_Item, interval=c(1, n), MAP=MAP, alpha=alpha)
	} else {
		re<-uniroot(Get_EffectiveNumberTest_Item_MidP, interval=c(1, n), MAP=MAP, alpha=alpha)
	}
	n_test<-ceiling(re$root)
	if(n_test > n){
		n_test = n
	}
	
	return(n_test)
	
}



#
#	
# QQ plot
#
QQPlot_Adj_Simple<-function(Pval, MAP, ntry=50){

	p<-length(Pval)
	id.miss<-which(is.na(MAP))
	if(length(id.miss) > 0){
		MAP[id.miss] = 0
	}	

	expected<-matrix(rep(0, p*ntry), ncol=ntry)
	for(i in 1:p){
		b1<-rbinom(ntry,1, MAP[i])
		expected[i,]<-runif(ntry, MAP[i], 1)*(1-b1) + b1*MAP[i]
	}

	expected<-apply(expected, 2, sort)
	expected1<-apply(expected, 1, median)

	max1<-max(-log10(expected1), -log10(Pval))
	qqplot(-log10(expected1), -log10(Pval), xlim=c(0,max1), ylim=c(0, max1))
	abline(0,1)

}




QQPlot_Adj<-function(Pval, MAP, main="QQ plot",ntry=500, confidence=0.95, Is.unadjsted=TRUE, Is.legend=TRUE, 
xlab="Expected Quantiles (-log10 P-values)", ylab="Observed Quantiles (-log10 P-values)"){

	#  QQPlot_Adj(out.tbl$P.value, out.tbl$MAP)
	#ntry<-500;confidence=0.95;main="QQ plot";confidence=0.95;Is.unadjsted=TRUE;Is.legend=TRUE;Pval=out1$results$P.value; MAP=out1$results$MAP
	alpha=1-confidence
	p<-length(Pval)
	Pval<-sort(Pval)
	id.miss<-union(which(is.na(MAP)), which(MAP < 0))
	if(length(id.miss) > 0){
		MAP[id.miss] = 0
	}	

	if(p > 500000){
		p1<-ceiling(p * 0.01)	# to save computation time...
		p2<-ceiling(p * 0.05)
		p3<-ceiling(p * 0.1)
		range<-c(1:p1, seq(p1+1, p2, 5), seq(p2+1, p3, 20), seq(p3+1, p, 50))
	} else if(p > 100000){
		p1<-ceiling(p * 0.01)	# to save computation time...
		p2<-ceiling(p * 0.05)
		p3<-ceiling(p * 0.1)
		range<-c(1:p1, seq(p1+1, p2, 2), seq(p2+1, p3, 10), seq(p3+1, p, 20))
	} else if (p > 10000){
		p1<-ceiling(p * 0.01)	# to save computation time...
		p2<-ceiling(p * 0.1)
		range<-c(1:p1, seq(p1+1, p2, 2), seq(p2+1, p, 10))
	}	else {
	
		range<-1:p
	}
	
	p2<-length(range)
	
	Pval1<-Pval[range]
	expected<-matrix(rep(0, p2*ntry), ncol=ntry)
	
	if(Is.unadjsted){
		expected_noadj<-matrix(rep(0, p2*ntry), ncol=ntry)
	}
	for(i in 1:ntry){
		b1<-rbinom(p,1, MAP)
		b2<-runif(p)
		temp<-(b2*(1-MAP) + MAP)*(1-b1) + b1*MAP
		expected[,i]<-sort(temp)[range]
		
		if(Is.unadjsted){
			expected_noadj[,i]<-sort(b2)[range]
		}
	}



	expected1<-apply(expected, 1, median)
	expected1.u<-apply(expected, 1, quantile,probs=(1-alpha/2))
	expected1.l<-apply(expected, 1, quantile,probs=alpha/2)

	if(Is.unadjsted){
		unadj.m<-apply(expected_noadj, 1, median)
	}
	

	max1<-max(-log10(expected1), -log10(Pval))
	if(Is.unadjsted){
		max1<-max(c(max1, -log10(unadj.m)))
	}
	plot(-log10(expected1), -log10(Pval1), xlim=c(0,max1), ylim=c(0, max1), col="black", main=main, pch="*"
	, xlab=xlab, ylab=ylab)
	points(-log10(expected1), -log10(expected1.u), type="l", lty="dashed", col="grey")
	points(-log10(expected1), -log10(expected1.l), type="l", lty="dashed", col="grey")
	
	if(Is.unadjsted){
		points(-log10(unadj.m), -log10(Pval1), xlim=c(0,max1), ylim=c(0, max1), col="grey", pch="*")
		if(Is.legend){
			legend("topleft", legend=c("MAP-adjusted", "Unadjusted"), col=c("black", "grey"), pch="*", bty="n")
		}
	}
	abline(0,1)

}



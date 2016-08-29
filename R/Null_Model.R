Get_SKAT_Residuals.Get_X1 = function(X1){
	
	qr1<-qr(X1)
	q1<-ncol(X1)
	if(qr1$rank < q1){
		
		X1.svd<-svd(X1)
		X1 = X1.svd$u	
	} 

	return(X1)

}


Get_SKAT_Residuals.linear = function(formula, data, n.Resampling, type.Resampling, id_include ){

	
 	mod = lm(formula, data=data)
	X1<-model.matrix(formula,data=data)
	X1<-Get_SKAT_Residuals.Get_X1(X1)
	
  	s2 = summary(mod)$sigma**2
  	res = mod$resid
	n1<-length(res)
	res.out<-NULL
	
	if(n.Resampling > 0){

		if(type.Resampling=="permutation"){
			res.out<-res %x% t(rep(1,n.Resampling))
			res.out<-apply(res.out,2,sample)
		} else if(type.Resampling=="bootstrap"){
			res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=sqrt(s2)),ncol=n.Resampling)
			X1_inv<-solve(t(X1) %*% X1)
			res.out<- res.out - (X1 %*% X1_inv) %*% (t(X1) %*% res.out)
		} else if(type.Resampling=="perturbation"){
			res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=1),ncol=n.Resampling)
			res.out<-res.out * res
			stop("Error: Perturbation is no more provided!")
		} else {
			stop("Error: Wrong resampling method!")
		}
	}

  	return(list(res=res, X1=X1,res.out=res.out,out_type="C", 
	n.Resampling=n.Resampling, type.Resampling=type.Resampling,
	id_include=id_include, s2=s2))
}


Get_SKAT_Residuals.logistic = function(formula, data, n.Resampling, type.Resampling,id_include){


 	mod = lm(formula, data)
	X1<-model.matrix(formula,data=data)
	X1<-Get_SKAT_Residuals.Get_X1(X1)
	
	glmfit= glm(formula, data=data, family = "binomial")
 	betas = glmfit$coef
  	mu    = glmfit$fitted.values
  	eta   = glmfit$linear.predictors
	n.case = sum(glmfit$y)

	pi_1 = mu*(1-mu)
  	res = glmfit$y- mu
	n1<-length(res)
	res.out<-NULL

	if(n.Resampling > 0){
		if(type.Resampling=="bootstrap.fast"){
		
			res.out<-Get_Resampling_Bin(n.case, mu, n.Resampling)
			if(is.null(res.out)){
				type.Resampling="bootstrap"
			}
			res.out<-res.out - mu
		} 
		
		if(type.Resampling=="permutation"){
			res.out1<-res %x% t(rep(1,n.Resampling))
			res.out<-apply(res.out1,2,sample)
		} else if(type.Resampling=="bootstrap"){
			mu1<-mu/sum(mu)	# all prob
			res.out<-matrix(rep(0,n.Resampling*n1),ncol=n.Resampling)
			for(i in 1:n.Resampling){
				#id_case<-sample(1:n1,n.case,prob=mu1)
				#res.out[id_case,i]<-1
				#res.out[,i]<-rbinom(n1,1,mu)

				res.out1<-rbinom(n1,1,mu)
				res.out2<-rbinom(n1,1,mu)

				id_case1<-which(res.out1 ==1)
				id_case2<-which(res.out2 ==1)

				id_c1<-intersect(id_case1,id_case2)
				id_c2<-union(setdiff(id_case1,id_case2),setdiff(id_case2,id_case1))
				if(n.case <= length(id_c1)){
					id_case<-sample(id_c1,n.case)
				}else if (n.case > length(id_c1) && n.case <= length(id_c1)+length(id_c2)){
					id_c3<-sample(id_c2,n.case - length(id_c1))
					id_case<-c(id_c1,id_c3)
				}else {					
					id_case3<-union(id_c1,id_c2)
					id_c4<-setdiff(1:n1,id_case3)
					n.needed<-n.case - length(id_case3)
					
					id_c5<-sample(id_c4,n.needed,prob=mu[id_c4])
					id_case<-union(id_case3,id_c5)
				}
				
				res.out[id_case,i]<-1
			}	
			res.out<-res.out - mu
		} else if(type.Resampling=="perturbation"){
			res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=1),ncol=n.Resampling)
			res.out<-res.out * res
			stop("Error: Perturbation is no more provided!")
		} else {
			if(is.null(res.out)){
				stop("Error: Wrong resampling method!")
			}
		}
		
	}

  	return(list(res=res, X1=X1,res.out=res.out,out_type="D", 
	n.Resampling=n.Resampling, type.Resampling=type.Resampling,
	id_include=id_include, mu=mu,pi_1=pi_1))

}

#
#	type	:  	
#		: permu - permutation
#		: bootstrap - bootstrap
#		: 
#
SKAT_Null_Model = function(formula, data=NULL, out_type="C", n.Resampling=0, type.Resampling="bootstrap", Adjustment=TRUE){
	
	SKAT_MAIN_Check_OutType(out_type)

	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)
	n1<-length(id_include)
	
	# Check whether n < 2000 and out_type="D", apply the adjustment 
	# if No_Adjustment = FALSE
	if(n1< 2000 && out_type=="D" && Adjustment){
		MSG<-sprintf("Sample size (non-missing y and X) = %d, which is < 2000. The small sample adjustment is applied!\n",n )
		cat(MSG)
		n.Resampling.kurtosis=10000
		#if(n > 1000){
		#	n.Resampling.kurtosis = floor(10000 - (n-1000) * 5)	
		#} 
		#if(n.Resampling.kurtosis < 5000){
		#	n.Resampling.kurtosis = 5000
		#}

		
		re<-SKAT_Null_Model_MomentAdjust(formula, data, n.Resampling, type.Resampling=type.Resampling, is_kurtosis_adj=TRUE, n.Resampling.kurtosis=n.Resampling.kurtosis)
		return(re)
	}


	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	if(out_type=="C"){
		re<-Get_SKAT_Residuals.linear(formula, data, n.Resampling, type.Resampling, id_include )
	} else {
		re<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
	}


	class(re)<-"SKAT_NULL_Model"
	re$n.all<-n
	return(re)
	
}


SKAT_Null_Model_MomentAdjust = function(formula, data=NULL, n.Resampling=0, type.Resampling="bootstrap", is_kurtosis_adj=TRUE, n.Resampling.kurtosis=10000){
	
	
	# check missing 
	obj1<-model.frame(formula,na.action = na.omit,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n<-dim(obj2)[1]
	n1<-dim(obj1)[1]
	id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)


	if(n - n1 > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
		warning(MSG,call.=FALSE)
	}

	re1<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
	re1$n.all<-n
	re2<-NULL

	if(is_kurtosis_adj == TRUE){
		re2<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling.kurtosis, type.Resampling, id_include )
	}

	class(re1)<-"SKAT_NULL_Model"
	re<-list(re1=re1, re2=re2, is_kurtosis_adj= is_kurtosis_adj, type = "binary")

	
	class(re)<-"SKAT_NULL_Model_ADJ"
	return(re)
	
}


SKAT_Null_Model_Get_Includes<-function( obj_omit, obj_pass){

	ID1<-rownames(obj_omit)
	ID2<-rownames(obj_pass)

	d1<-data.frame(ID=ID1)
	d2<-data.frame(ID=ID2, idx=1:length(ID2))

	d3<-merge(d1, d2,by.x="ID", by.y="ID")
	id_include = sort(d3$idx)

	return(id_include)
}




SKAT = function(Z,obj, kernel = "linear.weighted", method="davies", weights.beta=c(1,25)
, weights = NULL, impute.method = "fixed", r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf = 1, estimate_MAF=1){


	if(kernel != "linear" && kernel != "linear.weighted"){

		if(class(obj) == "SKAT_NULL_Model_ADJ"){
			msg<-sprintf("The small sample adjustment only can be applied for linear and linear.weighted kernel in the current version of SKAT! No adjustment is applied")
			warning(msg,call.=FALSE)
			obj<-obj$re1
		}

	}

	if(class(obj) == "SKAT_NULL_Model_EMMAX"){

		re = SKAT_emmaX(Z, obj, kernel = kernel, method=method, weights.beta=weights.beta, weights = weights, impute.method = impute.method,  r.corr=r.corr, is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)

	} else if(class(obj) == "SKAT_NULL_Model_ADJ"){

		re<-SKAT_With_NullModel_ADJ(Z, obj, kernel = kernel, method=method, weights.beta=weights.beta, weights = weights, impute.method = impute.method,  r.corr=r.corr, is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)

	} else if(class(obj) == "SKAT_NULL_Model"){

		re<-SKAT_With_NullModel(Z,obj, kernel = kernel, method=method, weights.beta=weights.beta, weights = weights, impute.method = impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype, is_dosage = is_dosage, missing_cutoff=missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)

	} else {
		#re<-SKAT_MAIN(Z,obj, ...)
		stop("The old interface is defunct! Please run SKAT_NULL_Model first!")
	}	
	class(re)<-"SKAT_OUT"
	return(re)
}


SKAT_1 = function(Z,obj, ...){

	if(class(obj) == "SKAT_NULL_Model_EMMAX"){
		re<-SKAT_emmaX( Z, obj, ...)
	} else if(class(obj) == "SKAT_NULL_Model_ADJ"){
		re<-SKAT_With_NullModel_ADJ(Z,obj, ...)
	} else if(class(obj) == "SKAT_NULL_Model"){
		re<-SKAT_With_NullModel(Z,obj, ...)
	} else {
		#re<-SKAT_MAIN(Z,obj, ...)
		stop("The old interface is defunct! Please run SKAT_NULL_Model first!")
	}	
	class(re)<-"SKAT_OUT"
	return(re)
}


#
#	Check the out_type
#
SKAT_MAIN_Check_OutType<-function(out_type){
 	
	if(out_type != "C" && out_type != "D"){
		stop("Invalid out_type!. Please use either \"C\" for the continous outcome or \"D\" for the dichotomous outcome.")
	}

}


SKAT_MAIN_Check_Z_Flip<-function(Z, id_include, Is.chrX=FALSE, SexVar=NULL){

	MAF<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){
		#msg<-sprintf("Genotypes of some variants are not the number of minor allele! It is fixed!")
		msg<-sprintf("Genotypes of some variants are not the number of minor alleles! These genotypes are flipped!")
		warning(msg,call.=FALSE)
		
		
		if(!Is.chrX){
			Z[,IDX.Err]<-2-Z[,IDX.Err]
		} else {
			id.male<-which(SexVar==1)
			id.female<-which(SexVar==2)
		
			if(length(id.male) > 0){
				Z[id.male,IDX.Err]<-1-Z[id.male,IDX.Err]
			} 
			if(length(id.female) > 0){
				Z[id.female,IDX.Err]<-2-Z[id.female,IDX.Err]
			}
		}
	}

	return(Z)
}

SKAT_MAIN_Check_Z_Impute<-function(Z, id_include,impute.method, SetID, Is.chrX=FALSE, SexVar=NULL){

	##################################################################
	# doing imputation


	MAF<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)
	MAF1<-SKAT_Get_MAF(Z, id_include=id_include, Is.chrX=Is.chrX, SexVar=SexVar)
	MAF_Org=MAF

	
	##########################################
	# Missing Imputation
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		if(is.null(SetID)){
			msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
		} else {
			msg<-sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
		}

		warning(msg,call.=FALSE)
		if(!Is.chrX){
			Z<-Impute(Z,impute.method)
		} else {
			Z<-Impute_XChr(Z, impute.method, SexVar)
		}
	} 

	#########################################
	# Check and recal
	
	Z<-SKAT_MAIN_Check_Z_Flip(Z, id_include, Is.chrX=Is.chrX, SexVar=SexVar)

	MAF<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)
	MAF1<-SKAT_Get_MAF(Z, id_include=id_include, Is.chrX=Is.chrX, SexVar=SexVar)
	IDX.Err<-which(MAF > 0.5)	
	if(length(IDX.Err) > 0){
		
		msg<-sprintf("ERROR! genotype flipping")
		stop(msg)
		
	}

	return(Z)

}


#
#	Check the Z, and do imputation
#
#


SKAT_MAIN_Check_Z<-function(Z, n, id_include, SetID, weights, weights.beta, impute.method, is_check_genotype
, is_dosage, missing_cutoff, max_maf, estimate_MAF=1, Is.chrX=FALSE, SexVar=NULL){

	#############################################
	# Check parameters

	if (class(Z)!= "matrix") stop("Z is not a matrix")
	if (nrow(Z)!=n) stop("Dimensions of y and Z do not match")
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}


	#####################################################
	# Check Z

	if(!is_check_genotype && !is_dosage){
		Z.test<-Z[id_include,]
		if(!is.matrix(Z.test)){
			Z.test<-as.matrix(Z.test)
		}
		return(list(Z.test=Z.test,weights=weights, MAF=rep(0, ncol(Z)), return=0) )
	}

	if(estimate_MAF==2){
		Z<-cbind(Z[id_include,])
		id_include<-1:length(id_include)
	}

	##############################################
	# Check Missing 

	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	m = ncol(Z)
	ID_INCLUDE_SNP<-NULL
	MAF_toCutoff<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)

	# Changed by SLEE, 07/21/2017
	for(i in 1:m){
		missing.ratio<-length(which(is.na(Z[,i])))/n
		sd1<-sd(Z[,i], na.rm=TRUE)
		if(missing.ratio < missing_cutoff && sd1 > 0){
			if(MAF_toCutoff[i] < max_maf){
				ID_INCLUDE_SNP<-c(ID_INCLUDE_SNP,i)
			}
		}
		
	}
	
	if(length(ID_INCLUDE_SNP) == 0){

		if(is.null(SetID)){
			msg<-sprintf("ALL SNPs have either high missing rates or no-variation. P-value=1")
		} else {
			msg<-sprintf("In %s, ALL SNPs have either high missing rates or no-variation. P-value=1",SetID )
		}
		warning(msg,call.=FALSE)
		
		re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 ) 
		return(re)		  

	} else if(m - length(ID_INCLUDE_SNP) > 0 ){

		if(is.null(SetID)){
			msg<-sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
		} else {
			msg<-sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
		}

		warning(msg,call.=FALSE)	
		Z<-as.matrix(Z[,ID_INCLUDE_SNP])
	}


	##################################################################
	# doing imputation

	MAF_Org<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)
	Z<-SKAT_MAIN_Check_Z_Impute(Z, id_include,impute.method, SetID, Is.chrX, SexVar)
	

	MAF<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=Is.chrX, SexVar=SexVar)
	MAF1<-SKAT_Get_MAF(Z, id_include=id_include, Is.chrX=Is.chrX, SexVar=SexVar)
	
	###########################################
	# Check non-polymorphic
	if(length(which(MAF1 > 0)) == 0){
		
		if(is.null(SetID)){
			msg<-sprintf("No polymorphic SNP. P-value = 1" )
		} else {
			msg<-sprintf("In %s, No polymorphic SNP. P-value = 1",SetID )
		}
		warning(msg,call.=FALSE)
		re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   
		return(re)
	}
	
	##########################################
	# Get Weights

	if(is.null(weights)){
		weights<-Beta.Weights(MAF,weights.beta)
	} else {
		weights = weights[ID_INCLUDE_SNP]
	}
	
	###########################################
	# Check missing of y and X

	if(n - length(id_include)  > 0){
	
		id_Z<-which(MAF1 > 0)

		if(length(id_Z) == 0){

			if(is.null(SetID)){
				msg<-sprintf("No polymorphic SNP. P-value = 1" )
			} else {
				msg<-sprintf("In %s, No polymorphic SNP. P-value = 1",SetID )
			}
			warning(msg,call.=FALSE)
			re<-list(p.value = 1, p.value.resampling =NA, Test.Type = NA, Q = NA, param=list(n.marker=0, n.marker.test=0), return=1 )   

		} else if (length(id_Z) == 1){
			Z<-cbind(Z[,id_Z])
		} else {
			Z<-Z[,id_Z]
		}

		if(!is.null(weights)){
			weights<-weights[id_Z]
		}

	}	
	
	if( dim(Z)[2] == 1){

		if(is.null(SetID)){
			msg<-sprintf("Only one SNP in the SNP set!" )
		} else {
			msg<-sprintf("In %s, Only one SNP in the SNP set!"
			,SetID )
		}
		# Suppress this warning
		#warning(msg,call.=FALSE)

		Z.test<-as.matrix(Z[id_include,])

	} else {

		Z.test<-Z[id_include,]

	}

	return(list(Z.test=Z.test,weights=weights, MAF=MAF_Org, id_include.test=id_include, return=0) )

}

SKAT_Check_RCorr<-function(kernel, r.corr){

	if(length(r.corr) == 1 && r.corr[1] == 0){
		return(1)
	}
	if(kernel != "linear" && kernel != "linear.weighted"){
		stop("Error: non-zero r.corr only can be used with linear or linear.weighted kernels")
	}

	for(i in 1:length(r.corr)){
		if(r.corr[i] < 0 || r.corr[i] > 1){
			stop("Error: r.corr should be >= 0 and <= 1")
		}
	}



}

SKAT_Check_Method<-function(method,r.corr, n=NULL, m=NULL){


	IsMeta=FALSE
	
	if(method != "liu"  && method != "davies" && method != "liu.mod" && method != "optimal" && method != "optimal.moment" 
	&& method != "optimal.mod" && method != "adjust" && method != "optimal.adj" && method != "optimal.moment.adj"
	 && method != "SKAT" && method != "SKATO" && method != "Burden" && method != "SKATO.m" 
	 && method !="davies.M" && method != "optimal.adj.M"  ){
		stop("Invalid method!")
	}
	
	# Run Meta-code
	if(method=="davies.M"){
		IsMeta=TRUE
		method="davies"	
	} else if(method=="optimal.adj.M"){
		IsMeta=TRUE
		method="optimal.adj"
	}

	if(method=="SKAT"){
		method="davies"	
		r.corr=0
	} else if(method == "SKATO"){
		method="optimal.adj"
	} else if(method =="SKATO.m"){
		method="optimal.moment.adj"
	} else if(method == "Burden"){
		method="davies"
		r.corr=1
	}
	
	
	if((method == "optimal" || method == "optimal.moment" ) && length(r.corr) == 1){
		r.corr = (0:10)/10
		#r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
	} else if( (method == "optimal.mod" || method == "optimal.adj" || method == "optimal.moment.adj" ) && length(r.corr)==1){
		r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
	}
	if(method =="optimal" ){
		method="davies"
	} else if (method =="optimal.moment") {
		method="liu.mod"
	}
	
	# if # of r.corr > 1, m/n < 1 and n > 10000, use Meta
	if(!is.null(n) && !is.null(m) && length(r.corr) > 1){
		if(m/n < 1 && n > 5000){
			IsMeta=TRUE
		}
	}
	

	re<-list(method=method,r.corr=r.corr, IsMeta=IsMeta)
	return(re)

}




SKAT_With_NullModel = function(Z, obj.res, kernel = "linear.weighted", method="davies", weights.beta=c(1,25), weights = NULL
, impute.method = "fixed", r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1
, SetID = NULL, out.z=NULL){

	

	n<-dim(Z)[1]
	m<-dim(Z)[2]
	out.method<-SKAT_Check_Method(method, r.corr, m=m, n=n)

	method=out.method$method
	r.corr=out.method$r.corr
	IsMeta=out.method$IsMeta

	SKAT_Check_RCorr(kernel, r.corr)
	
	# for old version
	if(is.null(obj.res$n.all)){
		obj.res$n.all=n
	}

	if(is.null(out.z)){
		out.z<-SKAT_MAIN_Check_Z(Z, obj.res$n.all, obj.res$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
	}
	
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}

	if(length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1){
		r.corr=0
		method="davies"
	} else if(length(r.corr) > 1 && sum(abs(out.z$Z.test - out.z$Z.test[,1])) == 0){
		r.corr=0
		method="davies"	
		
		msg<-sprintf("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!" )
		warning(msg,call.=FALSE)
		
	}


	if(obj.res$out_type == "C"){
		  #if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		  if( kernel =="linear" || kernel == "linear.weighted"){
		    re = SKAT.linear.Linear(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr, IsMeta=IsMeta)
		  } else {  
		    re = SKAT.linear.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights,obj.res$s2,method
			,obj.res$res.out, obj.res$n.Resampling)  
		  }
	} else if (obj.res$out_type == "D"){

		#if( (kernel =="linear" || kernel == "linear.weighted") && n > m){
		if( kernel =="linear" || kernel == "linear.weighted"){
			re = SKAT.logistic.Linear(obj.res$res, out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1,method
			,obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr, IsMeta=IsMeta)
		} else {  
			re = SKAT.logistic.Other(obj.res$res,out.z$Z.test
			,obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			,obj.res$res.out, obj.res$n.Resampling)  
		}
	}

	re$param$n.marker<-m
	re$param$n.marker.test<-dim(out.z$Z.test)[2]
	return(re)

}
 
#
#	Adjustment methods only use liu.mod, so it doesn't need method the "method" parameter
#	I use this field for outcome.type for subfunctions
#	
# 	estimate_MAF=1 using all samples, 2 only non-missing samples
SKAT_With_NullModel_ADJ = function(Z, obj.res.a, kernel = "linear.weighted", method="adjust", weights.beta=c(1,25), weights = NULL,
impute.method = "fixed", r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1
, SetID = NULL, out.z=NULL){

	
	n<-dim(Z)[1]
	m<-dim(Z)[2]
	obj.res<-obj.res.a$re1

	out.method<-SKAT_Check_Method(method,r.corr)
	method=out.method$method
	r.corr=out.method$r.corr

	SKAT_Check_RCorr(kernel, r.corr)
	# Use method field for the type of outcome
	method = obj.res.a$type

	
	# for old version
	if(is.null(obj.res$n.all)){
		obj.res$n.all=n
	}
	if(is.null(out.z)){
		out.z<-SKAT_MAIN_Check_Z(Z, obj.res$n.all, obj.res$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
	}
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}

	res2<-NULL
	if(obj.res.a$is_kurtosis_adj){
		res2<-obj.res.a$re2$res.out
	}	

	if(length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1){
		r.corr=0
		method="davies"
	} else if(length(r.corr) > 1 && sum(abs(out.z$Z.test - out.z$Z.test[,1])) == 0){
		r.corr=0
		method="davies"	
		msg<-sprintf("Rank of the genotype matrix is one! SKAT is conducted instead of SKAT-O!" )
		warning(msg,call.=FALSE)
	}



	if(length(r.corr) == 1){

		re = KMTest.logistic.Linear.VarMatching (obj.res$res,out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling,r.corr=r.corr
			, obj.res$mu, res.moments=res2)

	} else {

		re = SKAT_Optimal_Logistic_VarMatching(obj.res$res, out.z$Z.test
			, obj.res$X1, kernel, out.z$weights, obj.res$pi_1, method
			, obj.res$res.out, obj.res$n.Resampling, r.corr, obj.res$mu
			, res.moments=res2)

	}

	re$param$n.marker<-m
	re$param$n.marker.test<-dim(out.z$Z.test)[2]
	return(re)


}
 
 

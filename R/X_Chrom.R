
Impute_XChr<-function(Z, impute.method, SexVar){
	
	p<-dim(Z)[2]
	MAF<-SKAT_Get_MAF(Z, id_include=NULL, Is.chrX=TRUE, SexVar=SexVar)
	
	id.male<-which(SexVar==1)
	id.female<-which(SexVar==2)
	
	for(i in 1:p){
		IDX<-which(is.na(Z[,i]))
		IDX.m<-intersect(IDX, id.male)
		IDX.f<-intersect(IDX, id.female)
		
		maf1<-MAF[i]
		
		if(impute.method =="fixed"){
			if(length(IDX.m) > 0){
				Z[IDX.m,i]<- maf1
			}
			if(length(IDX.f) > 0){
				Z[IDX.f,i]<-2 * maf1
			}
		} else if(impute.method =="bestguess"){
			if(length(IDX.m) > 0){
				Z[IDX.m,i]<- round(maf1)
			}
			if(length(IDX.f) > 0){
				Z[IDX.f,i]<-round(2 * maf1)
			}
		
		} else {
			stop("Error: Imputation method shoud be \"fixed\" or \"bestguess\" ")
		}
		
	}

	return(Z)
}

#
#	X.inact : x inactivation
#
SKAT_Code_XChr<-function(Z, is_X.inact, SexVar){
	
	id.male<-which(SexVar==1)
	if(is_X.inact && length(id.male) > 0){
		Z[id.male,] = 2* Z[id.male,]
	}
	return(cbind(Z))

}

#
# Return Sex variable
#

Check_Sex_Var = function(formula, SexVar, data=NULL){

	# Check SexVar (1=male, 2=female)
	temp.dat<-get_all_vars(formula, data)
	idx.var.sex = which(colnames(temp.dat) == SexVar)
	if(length(idx.var.sex) == 0){
		msg<-sprintf("No variable %s exists in the formula!\n", SexVar);
		stop(msg)
	}
	
	# check missing
	Sex.Var<-temp.dat[,idx.var.sex]
	Sex.Var.Org<-Sex.Var
	id.missing<-which(is.na(Sex.Var))
	if(length(id.missing) == length(Sex.Var)){
		msg<-sprintf("All sex values are NA\n");
		stop(msg);
	} else if(length(id.missing) > 0){
		msg<-sprintf("There are %d missing values in the Sex variable\n", length(id.missing));
		cat(msg);
		Sex.Var<-Sex.Var[id.missing]
	}
	
	# check whether SexVar is 1, 2 variable
	
	id<-intersect(which(Sex.Var!=1), which(Sex.Var!=2))
	if(length(id) >0){
		msg<-sprintf("%d subjects have sex values neither 1 nor 2. They should be either 1 (=male) or 2(=female)!", length(id))
		stop(msg);
	}
	
	return(Sex.Var.Org)
}



SKAT_Null_Model_ChrX = function(formula, SexVar, data=NULL, out_type="C", n.Resampling=0, type.Resampling="bootstrap", Adjustment=TRUE){
	
	SKAT_MAIN_Check_OutType(out_type)
	SexVar=Check_Sex_Var(formula, SexVar, data=NULL)

	
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
		re<-SKAT_Null_Model_MomentAdjust(formula, data, n.Resampling, type.Resampling=type.Resampling, is_kurtosis_adj=TRUE, n.Resampling.kurtosis=n.Resampling.kurtosis)
		
		class(re)<-"SKAT_NULL_Model_Adj_ChrX"
		re$SexVar=SexVar
		re$n.all<-n
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


	class(re)<-"SKAT_NULL_Model_ChrX"
	re$SexVar=SexVar
	re$n.all<-n
	return(re)
	
}



SKAT_ChrX<-function(Z, obj, is_X.inact =TRUE, kernel = "linear.weighted", method="davies", weights.beta=c(1,25)
, weights = NULL, impute.method = "fixed", r.corr=0, is_check_genotype=TRUE
, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID=NULL){

	
	if(kernel != "linear" && kernel != "linear.weighted"){

		if(class(obj) == "SKAT_NULL_Model_Adj_ChrX"){
			msg<-sprintf("The small sample adjustment only can be applied for linear and linear.weighted kernel in the current version of SKAT! No adjustment is applied")
			warning(msg,call.=FALSE)
			obj<-obj$re1
		}

	}
	
	if(class(obj) != "SKAT_NULL_Model_Adj_ChrX" && class(obj) != "SKAT_NULL_Model_ChrX"){
		msg<-sprintf("The obj parameter should be a returned object from SKAT_Null_Model_ChrX!")
		stop(msg)
	}
	

	SexVar = obj$SexVar
	n.all=obj$n.all

 	if(class(obj) == "SKAT_NULL_Model_Adj_ChrX"){
		
		id_include = obj$re1$id_include
	

	} else if(class(obj) == "SKAT_NULL_Model_ChrX"){

		id_include = obj$id_include

	} else {
		#re<-SKAT_MAIN(Z,obj, ...)
		stop("The old interface is defunct! Please run SKAT_NULL_Model first!")
	}	

	out.z<-SKAT_MAIN_Check_Z(Z, n.all, id_include, SetID
		, weights, weights.beta, impute.method, is_check_genotype
		, is_dosage, missing_cutoff,max_maf=max_maf, estimate_MAF=estimate_MAF, Is.chrX=TRUE, SexVar=SexVar)

	out.z$Z.test<-SKAT_Code_XChr(out.z$Z.test, is_X.inact=is_X.inact, SexVar=SexVar[out.z$id_include.test])
	
	
 	if(class(obj) == "SKAT_NULL_Model_Adj_ChrX"){
		
	
		re<-SKAT_With_NullModel_ADJ(Z, obj, kernel = kernel, method=method, weights.beta=weights.beta
		, weights = weights, impute.method = impute.method,  r.corr=r.corr, is_check_genotype=is_check_genotype
		, is_dosage = is_dosage, missing_cutoff=missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF, out.z=out.z)

	} else if(class(obj) == "SKAT_NULL_Model_ChrX"){

		re<-SKAT_With_NullModel(Z,obj, kernel = kernel, method=method, weights.beta=weights.beta
		, weights = weights, impute.method = impute.method, r.corr=r.corr, is_check_genotype=is_check_genotype
		, is_dosage = is_dosage, missing_cutoff=missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF, out.z=out.z)

	} 
		
	class(re)<-"SKAT_OUT"
	return(re)

}

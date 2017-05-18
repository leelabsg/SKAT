
#
#	Read Fam File
#
Read_Plink_FAM<-function(Filename, Is.binary=TRUE, flag1=0){

	Check_File_Exists(Filename)
	Fam.Obj<-read.table(Filename, header=FALSE)
	colnames(Fam.Obj)<-c("FID","IID","PID","MID", "Sex", "Phenotype")

	id.missing<-which(Fam.Obj$Phenotype == -9)
	if(length(id.missing) > 0){
		Fam.Obj$Phenotype[id.missing]<-NA
	}

	if(Is.binary && flag1 == 0){
		
		id_0<-which(Fam.Obj$Phenotype == 0)
		id_1<-which(Fam.Obj$Phenotype == 1)
		id_2<-which(Fam.Obj$Phenotype == 2)

		if(length(id_0)> 0){
			Fam.Obj$Phenotype[id_0] <-NA
		}
		if(length(id_1)> 0){
			Fam.Obj$Phenotype[id_1] <-0
		}
		if(length(id_2)> 0){
			Fam.Obj$Phenotype[id_2] <-1
		}

	}

	return(Fam.Obj)

}

Read_Plink_FAM_Cov<-function(Filename, File_Cov, Is.binary=TRUE, flag1=0, cov_header=TRUE){

	Check_File_Exists(Filename)
	Check_File_Exists(File_Cov)
	Fam.Obj<-read.table(Filename, header=FALSE)
	
	colnames(Fam.Obj)<-c("FID","IID","PID","MID", "Sex", "Phenotype")
	
	Cov.Obj<-read.table(File_Cov, header=cov_header)
	ncov<-dim(Cov.Obj)[2]
	
	if(ncov <= 2){
		msg<-sprintf("Error: Cov file only has <= 2 columns!\n")
		stop(msg)
	}
		
	if(cov_header == FALSE){
		
		colnames(Cov.Obj)<-c("FID", "IID", sprintf("COV%d",1:(ncov-2)))
	}
	
	if(colnames(Cov.Obj)[1] != "FID" || colnames(Cov.Obj)[2] != "IID"){
		msg<-sprintf("Error: the first two columns of Cov file must be labelled as FID and IID. If the Cov file does not have a header row, please set cov_header=FALSE!\n")
		
		stop(msg) 
	}
	
	for(i in 3:ncov){
		# check numeric coding
		if(is.numeric(Cov.Obj[,i]) == FALSE){
			msg<-sprintf("Error: column [%d] in Cov file isn't numerically coded! If the Cov file has a header row, please set cov_header=TRUE!", i)
			stop(msg)
		}
		id.missing<-which(Cov.Obj[,i] == -9)
		if(length(id.missing) > 0){
			Cov.Obj[id.missing,i]<-NA
		}
	}

	id.missing<-which(Fam.Obj$Phenotype == -9)
	if(length(id.missing) > 0){
		Fam.Obj$Phenotype[id.missing]<-NA
	}
	if(Is.binary && flag1 == 0){
		
		id_0<-which(Fam.Obj$Phenotype == 0)
		id_1<-which(Fam.Obj$Phenotype == 1)
		id_2<-which(Fam.Obj$Phenotype == 2)

		if(length(id_0)> 0){
			Fam.Obj$Phenotype[id_0] <-NA
		}
		if(length(id_1)> 0){
			Fam.Obj$Phenotype[id_1] <-0
		}
		if(length(id_2)> 0){
			Fam.Obj$Phenotype[id_2] <-1
		}

	}

	# Combine
	Fam.Obj$idx_fam<-1:(dim(Fam.Obj)[1])
	ID1<-paste(as.character(Fam.Obj$FID), as.character(Fam.Obj$IID))
	ID2<-paste(as.character(Cov.Obj$FID), as.character(Cov.Obj$IID))
	
	Fam.Obj$MergeID_SKAT=ID1
	Cov.Obj$MergeID_SKAT=ID2
	
	Cov.Obj$FID=NULL
	Cov.Obj$IID=NULL
	
	Out<-merge(Fam.Obj, Cov.Obj, by.x="MergeID_SKAT", by.y="MergeID_SKAT", all.x=TRUE, all.y=FALSE, sort=FALSE)
	Out<-Out[order(Out$idx_fam),]
	Out$idx_fam<-NULL
	Out$MergeID_SKAT<-NULL
	
	return(Out)

}


#
# x is either y or SKAT_NULL_Model 
#
SKAT.SSD.OneSet = function(SSD.INFO, SetID, obj, ..., obj.SNPWeight=NULL){
	
	id1<-which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	SetIndex<-SSD.INFO$SetInfo$SetIndex[id1]
	re = SKAT.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=obj.SNPWeight)

	return(re)
}

SKAT.SSD.GetSNP_Weight<-function(SSD.INFO, SetIndex, obj.SNPWeight=NULL){


	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]
	
	is_ID = FALSE
	if(!is.null(obj.SNPWeight)){
		is_ID = TRUE
	}

	try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=is_ID),silent = TRUE)
	if(class(try1) != "try-error"){
		Z<-try1
		Is.Error<-FALSE	
	} else {
		err.msg<-geterrmessage()
		msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
		stop(msg)
	}

	if(!is_ID){
		re=list(Z=Z, Is.weights=FALSE)
		return(re)
	}


	SNP_ID<-colnames(Z)
	p<-ncol(Z)
	weights<-rep(0, p)
	for(i in 1:p){
		val1<-SNP_ID[i]			
		val2<-obj.SNPWeight$hashset[[val1]]
			
		if(is.null(val2)){
			msg<-sprintf("SNP %s is not found in obj.SNPWeight!", val1)
			stop(msg)
		}

		weights[i]<-val2
	}
	
	# Change FALSE to TRUE
	re=list(Z=Z, Is.weights=TRUE, weights=weights)
	return(re)
		

}


SKAT.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){
	
	re1 = SKAT.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.SNPWeight=obj.SNPWeight)

	
	if(!re1$Is.weights){
		re<-SKAT(re1$Z, obj, ...)
	} else {
	
		re<-SKAT(re1$Z, obj, weights=re1$weights, ...)
	}
	
	return(re)
}



#
# x is either y or SKAT_NULL_Model 
#
SKAT.SSD.OneSet_SetIndex_OLD = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){
	
	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	#id1 = SetIndex
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]
	
	is_ID = FALSE
	if(!is.null(obj.SNPWeight)){
		is_ID = TRUE
	}
	try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=is_ID),silent = TRUE)
	if(class(try1) != "try-error"){
		Z<-try1
		Is.Error<-FALSE	
	} else {
		err.msg<-geterrmessage()
		msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
		stop(msg)
	}
		
	
	if(is.null(obj.SNPWeight)){

		re<-SKAT(Z, obj, ...)
	} else {
	
		SNP_ID<-colnames(Z)
		p<-ncol(Z)
		weights<-rep(0, p)
		for(i in 1:p){
			val1<-SNP_ID[i]			
			val2<-obj.SNPWeight$hashset[[val1]]
			
			if(is.null(val2)){
				msg<-sprintf("SNP %s is not found in obj.SNPWeight!", val1)
				stop(msg)
			}

			weights[i]<-val2
		}
		re<-SKAT(Z, obj, weights=weights, ...)
	}
	
	return(re)
}


#
# Only SKAT_Null_Model obj can be used
#
SKAT.SSD.All = function(SSD.INFO, obj, ..., obj.SNPWeight=NULL){
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	Is.Resampling = FALSE
	n.Resampling = 0
	
	if(class(obj) == "SKAT_NULL_Model"){
		if(obj$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$n.Resampling

			OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	} else if(class(obj) == "SKAT_NULL_Model_ADJ"){
		if(obj$re1$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$re1$n.Resampling

			OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	}

	pb <- txtProgressBar(min=0, max=N.Set, style=3)
	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1 = try(SKAT.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, obj=obj, ..., obj.SNPWeight=obj.SNPWeight))

		
		if(class(try1) != "try-error"){
			re<-try1
			Is.Error<-FALSE
		} else {

			err.msg<-geterrmessage()
			msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
			
		}
		
		if(!Is.Error){

			OUT.Pvalue[i]<-re$p.value
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
			if(Is.Resampling){
				OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
			}
		}
		
		#if(floor(i/100)*100 == i){
		#	cat("\r", i, "/", N.Set, "were done");
		#}
		setTxtProgressBar(pb, i)

	}
	close(pb)

	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
	class(re)<-"SKAT_SSD_ALL"

	return(re)	
}


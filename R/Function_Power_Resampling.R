
###########################################
#	Function related to resampling 


Get_Resampling_Pvalue<-function(obj){

	if(class(obj) != "SKAT_OUT"){
		stop("obj is not a SKAT output object")
	}

	if(is.null(obj$p.value.resampling)){
		stop("No resampling was applied!")
	}

	n<-length(obj$p.value.resampling)
	n1<-length(which(obj$p.value >= obj$p.value.resampling))
	pval1<-(n1+1)/(n+1)
	
	re<-list(p.value=pval1, is_smaller=FALSE)
	if(n1==0){
		re$is_smaller=TRUE
	}
	
	return(re)
}

Get_Resampling_Pvalue_1<-function(p.value,p.value.resampling){

	if(is.null(p.value.resampling)){
		stop("No resampling was applied!")
	}

	n<-length(p.value.resampling)
	n1<-length(which(p.value >= p.value.resampling))
	pval1<-(n1+1)/(n+1)

	re<-list(p.value=pval1, is_smaller=FALSE)
	if(n1==0){
		re$is_smaller=TRUE
	}

	return(re)
}


###########################################
#	Function related to FWER 



Resampling_FWER<-function(obj,FWER=0.05){

	if(class(obj) != "SKAT_SSD_ALL" && class(obj) !="SKATBinary_SSD_ALL"){
		stop("obj is not a SKAT.SSD.All output object")
	}
	p.min<-apply(obj$P.value.Resampling,2,min,na.rm=TRUE)
	P.cut<-quantile(p.min,FWER)
	ID<-which(obj$results$P.value < P.cut)

	if(length(ID) == 0){
		re<-list(result=NULL, n=length(ID) ,ID=NULL)
	} else {
		re<-list(result=obj$results[ID,], n=length(ID) ,ID=ID)
	}
	return(re) 
	
}

Resampling_FWER_1<-function(P.value, P.value.Resampling, FWER=0.05){

	if(is.matrix(P.value.Resampling) == FALSE){
		stop("P.value.Resampling should be a matrix")
	}
	p.min<-apply(P.value.Resampling,2,min,na.rm=TRUE)
	P.cut<-quantile(p.min,FWER)
	ID<-which(P.value < P.cut)

	if(length(ID) == 0){
		re<-list(result=NULL, n=length(ID) ,ID=NULL)
	} else {
		re<-list(result=P.value[ID], n=length(ID) ,ID=ID)
	}
	return(re) 	

}


############################################################
#
#	Power Estimation

Get_RequiredSampleSize<-function (obj, Power=0.8){



	if(class(obj) == "SKAT_Power"){
		re<-Get_RequiredSampleSize.SKAT_Power(obj, Power)
	} else {
		re<-Get_RequiredSampleSize.numeric(obj, Power)
	}


	return(re)
}


Get_RequiredSampleSize.SKAT_Power<-function(obj, Power=0.8){

	Get_RequiredSampleSize.numeric(obj$Power, Power)

}


Get_RequiredSampleSize.numeric<-function(Power.Est, Power=0.8){

	N.Sample.ALL<-as.numeric(rownames(Power.Est))
	alpha<-as.numeric(colnames(Power.Est))

	re<-list()
	for(i in 1:length(alpha)){

		
		temp<-which(Power.Est[,i] > Power)
		if(length(temp) == 0){
			temp<-sprintf("> %d",max(N.Sample.ALL))
			#print(temp)
			re[[i]]<-temp
		} else if( min(temp) ==1 ){
			re[[i]]<-min(N.Sample.ALL)
		} 
		else {
			id1<-min(temp)
			re[[i]]<-(N.Sample.ALL[id1] - N.Sample.ALL[id1-1])/(Power.Est[id1,i] - Power.Est[id1-1,i]) * (Power - Power.Est[id1-1,i]) + N.Sample.ALL[id1-1]

		}
	}
	names(re)<-sprintf("alpha = %.2e",alpha)
	return(re)
}







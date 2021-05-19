

SKAT_CommonRare_Robust.SSD.OneSet = function(SSD.INFO, SetID, obj, ..., obj.SNPWeight=NULL){
  
  id1<-which(SSD.INFO$SetInfo$SetID == SetID)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
    stop(MSG)
  }	
  SetIndex<-SSD.INFO$SetInfo$SetIndex[id1]
  
  re = SKAT_CommonRare_Robust.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=obj.SNPWeight)
  
  return(re)
}

SKAT_CommonRare_Robust.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){
  
  re1 = SKAT.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.SNPWeight=obj.SNPWeight)
  SetID = SSD.INFO$SetInfo$SetID[SetIndex]
  if(!re1$Is.weights){
    re<-SKAT_CommonRare_Robust(re1$Z, obj, ...)
  } else {
    
    re<-SKAT_CommonRare_Robust(re1$Z, obj, weights=re1$weights, ...)
  }
  
  return(re)
  
}




#
# Only SKAT_Null_Model obj can be used
#
SKAT_CommonRare_Robust.SSD.All = function(SSD.INFO, obj, ...,obj.SNPWeight=NULL){
  
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.Pvalue.Resampling<-NULL
  OUT.Q<-rep(NA,N.Set)
  OUT.snp.mac<-list()
  
  OUT.nRare<-rep(NA,N.Set)
  OUT.nCommon<-rep(NA,N.Set)
  
  Is.Resampling = FALSE
  n.Resampling = 0
  
  if(Check_Class(obj, "SKAT_NULL_Model")){
    if(obj$n.Resampling > 0){
      Is.Resampling = TRUE
      n.Resampling = obj$n.Resampling
      
      OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
    }
  }
  pb <- txtProgressBar(min=0, max=N.Set, style=3)
  for(i in 1:N.Set){
    Is.Error<-TRUE
    try1<-try(SKAT_CommonRare_Robust.SSD.OneSet_SetIndex(SSD.INFO, i, obj, ..., obj.SNPWeight=obj.SNPWeight) ,silent = TRUE)
    if(Is_TryError(try1)){
      
      err.msg<-geterrmessage()
      msg<-sprintf("Error to run SKATBinary for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
      
    } else {
	  ##Changed by Zhangchen 05/19/2021, remove the situation where there is no variation in the gene. 
      if (max(try1$mac.rare+try1$mac.common,0,na.rm=T)>=1){
     	 re<-try1
     	 OUT.Pvalue[i]<-re$p.value
     	 OUT.Marker[i]<-re$param$n.marker
     	 OUT.Marker.Test[i]<-re$param$n.marker.test
      	 OUT.nRare[i]<-re$n.rare
     	 OUT.nCommon[i]<-re$n.common
     	 OUT.Q[i]<-re$Q
     # if(Is.Resampling){
     #   OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
     # }
      	SetID<-SSD.INFO$SetInfo$SetID[i]
      	OUT.snp.mac[[SetID]]<-re$test.snp.mac
      }
    }
    #if(floor(i/100)*100 == i){
    #	cat("\r", i, "/", N.Set, "were done");
    #}
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)	
  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, Q=OUT.Q
                      , N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test, N.Marker.Rare=OUT.nRare, N.Marker.Common=OUT.nCommon)
  re<-list(results=out.tbl,OUT.snp.mac=OUT.snp.mac)
  class(re)<-"SKAT_SSD_ALL"
  
  return(re)	
}



colMax <- function(data) apply(data,2, max, na.rm = TRUE)

##Changed by Zhangchen 09/06/2020, the default method is changed to SKAT instead of SKATO. 
SKAT_CommonRare_Robust<-function(Z, obj, kernel = "linear.weighted", method="SKAT"
                            , r.corr=NULL, weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weights = NULL
                            , CommonRare_Cutoff=NULL, impute.method = "bestguess",is_check_genotype=TRUE
                            ,is_dosage = FALSE, missing_cutoff=0.15, max_maf=1
                            , estimate_MAF=1){
  
  # Added by SLEE 12/23/2019. Currently only joint is used
  test.type="Joint"
  
  
  SetID1=NULL
  # This function only can be used for SNPs
  #is_check_genotype=TRUE
  
  if(Check_Class(obj, "SKAT_NULL_Model_ADJ")){
    
    obj.res=obj$re1
    
  } else if(Check_Class(obj, "SKAT_NULL_Model")){
    
    obj.res=obj
    
  } else {
    stop("Wrong obj!")
  }
  
  # changed by SLEE 12/23/2019
  if(!Check_Class(Z,  c("matrix", "dgCMatrix", "dgeMatrix"))){
    stop("Z should be a matrix")
  }
  
  # Compute common and rare
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  m.org<-m
  
  
  if(is.null(CommonRare_Cutoff)){
    CommonRare_Cutoff<-1/sqrt(n * 2)
  }
  # Check Cutoff
  if(CommonRare_Cutoff < 0 && CommonRare_Cutoff > 0.5){
    stop("Error in CommonRare_Cutoff! It should be NULL or a numeric value between 0 to 0.5")
  }
  
  # for old version
  if(is.null(obj.res$n.all)){
    obj.res$n.all=n
  }
  # no max_maf cutoff
  
  
  if (kernel=="linear"){weights=rep(1,ncol(Z)) } else{if (kernel != "linear.weighted"){stop("Wrong kernel!")}}
  if (!is.null(weights)){if (length(weights)!=ncol(Z)) {stop("Incorrect length of weights!")}} else {
    if (length(weights.beta.rare)!=2){stop("Incorrect length of weights.beta.rare!")}
    if (length(weights.beta.common)!=2){stop("Incorrect length of weights.beta.common!")}
    
  }

  # changed by SLEE 12/23/2019, removed SKAT::: 
  out<-SKAT_MAIN_Check_Z(Z, obj.res$n.all, id_include=obj.res$id_include, SetID=SetID1, weights=weights, weights.beta=c(1,1), 
                         impute.method="fixed", is_check_genotype=is_check_genotype, is_dosage=is_dosage, missing_cutoff, max_maf= max_maf, estimate_MAF=estimate_MAF)
  if(out$return ==1){
    out$param$n.marker<-m
    out$n.rare = 0
    out$n.common = 0
    out$test.type= test.type
    out$Cutoff = CommonRare_Cutoff
    
    return(out)
  }
  
  Z.org<-Z
  Z<-out$Z.test
  ##weights.org<-weights
  ##weights<-out$weights
  for (jj in 1:ncol(Z)) {
    n.g <- sum(Z[, jj])
    if (n.g/(2 * length(Z[, jj])) > 0.5) {
      Z[, jj] <- 2 - Z[, jj]
      n.g <- sum(Z[, jj])
    }
  }  
  # Since I already used ID include.
  obj.res$n.all =nrow(Z) 
  obj.res$id_include = 1:nrow(Z)	

  mafcutoff=CommonRare_Cutoff

  # changed by SLEE 12/23/2019, removed SKAT::: 
  MAF<-Get_MAF(Z)
  maf_temp = which(MAF <= mafcutoff)
  
  if (is.null(weights)){
    weight=rep(0,length(MAF))
    if (length(maf_temp) > 0 & length(maf_temp) < length(MAF)) {
      weight[maf_temp] = Beta_Weight(MAF[maf_temp], weights.beta.rare)
      weight[-maf_temp] =Beta_Weight(MAF[-maf_temp], weights.beta.common)
      Z1=Z[,maf_temp]
      Z2=Z[,-maf_temp]
      pi_1=obj$pi_1
      X1= obj$X1
      Z1.1<- (Z1 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z1 * pi_1))
		  Z2.1<- (Z2 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z2 * pi_1))
  	  temp1<-t(Z1.1) %*% Z1.1
	    temp2<-t(Z2.1) %*% Z2.1    
      z1.var<-sum(temp1 * temp1)
	    z2.var<-sum(temp2 * temp2)
      weight[maf_temp] <- weight[maf_temp]/sqrt(sqrt(z1.var))
	    weight[-maf_temp] <- weight[-maf_temp]/sqrt(sqrt(z2.var))
    }    else {
      if (length(maf_temp) == 0) {
        weight = Beta_Weight(MAF,weights.beta.common)
      }
      if (length(maf_temp) == length(MAF)) {
        weight = Beta_Weight(MAF,  weights.beta.rare)
      }
    }
    weights=weight
  }
  

  
  if(Check_Class(obj, "SKAT_NULL_Model_ADJ")){
    obj$re1$id_include = obj.res$id_include
    obj$re1$n.all = obj.res$n.all
  } else {
    obj$id_include = obj.res$id_include
    obj$n.all = obj.res$n.all
  }

  colmax_Z=colMax(Z)
  list_tiny=which(colmax_Z<=0.2)
  if (length(list_tiny)>=1 ){
    if (length(list_tiny)<dim(Z)[2]){
      # changed by SLEE 12/23/2019, weigths changed to weights
      Z=Z[,-list_tiny];weights=weights[-list_tiny];
      }else { stop("all genotypes are close to 0!")
    }
  }
  
  if (max(Z)>2 | min(Z)<0) {stop("Z is out of bounds[0,2]!")}
  
  m.test<-ncol(Z)
  
  # changed by SLEE 12/23/2019, removed SKAT::: 
  MAF<-Get_MAF(Z)  
    
  id.rare<-intersect(which(MAF < CommonRare_Cutoff), which(MAF > 0))
  id.common<-intersect(which(MAF >= CommonRare_Cutoff), which(MAF > 0))
  
  n.rare= length(id.rare)
  n.common= length(id.common)
  if (length(id.rare)==0){mac.rare=0}else {mac.rare=sum(Z[,id.rare])}
  if (length(id.common)==0){mac.common=0}else {mac.common=sum(Z[,id.common])}
  
  
  is.run<-FALSE
  
  if ( method=="SKAT"){r.corr=0}
  if ( method=="Burden"){r.corr=1}
  re<-SKATBinary_spa(G=Z,obj=obj.res,weights = weights, method=method, r.corr=r.corr)
  
  if (length(r.corr)==0){r.corr= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)}
  is.run=TRUE
  
  re$param$n.marker<-m.org
  re$param$n.marker.test<-m.test
  re$test.snp.mac<-SingleSNP_INFO(out$Z.test)
  
  re$param$rho=r.corr
  if (method=="SKATO"){
    re$param$minp=min( re$p.value_each)
    re$param$rho_est=r.corr[which.min(re$p.value_each)]
  }
  re$n.rare = n.rare
  re$mac.rare=mac.rare
  re$n.common = n.common
  re$mac.common=mac.common
  re$Cutoff = CommonRare_Cutoff
  
  return(re)
  
}

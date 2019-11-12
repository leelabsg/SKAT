Beta_Weight<-function(MAF,weights.beta){
  
  n<-length(MAF)
  weights<-rep(0,n)	
  IDX_0<-which(MAF == 0)
  if(length(IDX_0) == n){
    stop("No polymorphic SNPs")
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  
  
  #print(length(IDX_0))
  #print(weights[-IDX_0])
  return(weights)
  
}

SPA_ER_kernel<-function(G, obj,  u, Cutoff, variancematrix, weight){
  
  res_all=obj$res
  MAFsum=colSums(G)
  qtemp=matrix(rep(0,dim(G)[2]))
  res_time=1
  
  zscore.all_0 <- matrix(rep(0, (ncol(G) * res_time)), ncol = ncol(G))
  zscore.all_1 <- matrix(rep(0, (ncol(G) * res_time)), ncol = ncol(G))
  VarS = c()
  p.old = c()
  p.new = c()
  g.sum = 0
  q.sum = 0
  for (jj in 1:ncol(G)) {
    n.g <- sum(G[, jj])
    if (n.g/(2 * length(G[, jj])) > 0.5) {
      G[, jj] <- 2 - G[, jj]
      n.g <- sum(G[, jj])
    }
    NAset <- which(G[, jj] == 0)
    G1 <- G[, jj] - obj$XXVX_inv %*% (obj$XV %*% G[, jj])
    q <- (t(G1) %*% (res_all + u))/sqrt(n.g)
    g = G1/sqrt(n.g)
    mu.qtemp = u
    g.qtemp = g
    mu1 <- sum(mu.qtemp * g.qtemp)
    var1 <- sum(mu.qtemp * (1 - mu.qtemp) * g.qtemp^2)
    stat.qtemp <- (q - mu1)^2/var1
    p_temp1 <- pchisq(stat.qtemp, lower.tail = FALSE, df = 1)
    p.old[jj] = p_temp1
    
    zscore.all_0[, jj] = (q - mu1) * sqrt(n.g)
    id1 <- which(stat.qtemp > Cutoff^2)
    qtemp[jj, ] = q
    if (MAFsum[jj] < 10) {
      if (length(id1) > 0) {
        G_temp = G[, jj]
        G_temp[which(G_temp <= 0.2)] = 0
        temp_binary = SKATBinary(as.matrix(G_temp), obj,
                                 method.bin = "ER")
        p_temp_binary = c(temp_binary$p.value, temp_binary$p.value.resampling)
        p_temp1[id1] = p_temp_binary[id1]
      }
    }         else {
      if (length(id1) > 0) {
        qtemp_id1 = qtemp[jj, id1]
        p_temp1_id1 = c()
        for (jjk in 1:length(qtemp_id1)) {
          p_temp1_id1[jjk] = SPAtest:::Saddle_Prob(qtemp_id1[jjk],
                                                   mu = mu.qtemp, g = g.qtemp, Cutoff = Cutoff,
                                                   alpha = 5 * 10^-8)$p.value
          if (p_temp1_id1[jjk] != 0) {
            p_temp1[id1[jjk]] = p_temp1_id1[jjk]
          }
        }
      }
    }
    p.new[jj] = p_temp1
    if (variancematrix[jj, jj] <= 0) {
      zscore.all_1[, jj] = 0
    }         else {
      zscore.all_1[, jj] = qnorm(p_temp1/2, mean = 0, sd = sqrt(variancematrix[jj,
                                                                               jj]), lower.tail = FALSE, log.p = FALSE) * weight[jj] *
        sign(q - mu1)
    }
    if (p_temp1 > 0) {
      VarS[jj] = zscore.all_0[, jj]^2/qchisq(p_temp1, 1,
                                             ncp = 0, lower.tail = FALSE, log.p = FALSE)
    }         else {
      VarS[jj] = zscore.all_0[, jj]^2/500
    }
    
    if (p_temp1<1){        
      g.sum = g.sum + g.qtemp * weight[jj] * sqrt(n.g)
      q.sum = q.sum + q * weight[jj] * sqrt(n.g)
    }
  }
  outlist = list()
  outlist$p.old = p.old
  outlist$p.new = p.new
  outlist$zscore.all_0 = zscore.all_0
  outlist$VarS = VarS
  outlist$mu = mu.qtemp
  outlist$g.sum = g.sum
  outlist$q.sum = q.sum
  outlist$zscore.all_1 = zscore.all_1
  return(outlist)
}


SKATBinary_spa<-function (G, obj, weights, method="SKATO",r.corr=NULL){
  Cutoff=2
  if (length(G)==0) {stop("WARNING: no-variantion in the whole genotype matrix!\n")}
  X = obj$X1
  u = obj$mu
  w = obj$pi_1
  obj$XV = t(X * w)
  temp1 = solve(t(X) %*% (X * w))
  obj$XXVX_inv = X %*% temp1

  
  MAF_0 = which(colSums(G)==0)
  if (length(MAF_0)>0){
    cat("The following columns are removed due to no-variation: ", MAF_0,"\n")
    G=Matrix(G[,-MAF_0],sparse=TRUE)
  }
  if (length(G)==0) {stop("WARNING: no-variantion in the whole genotype matrix!\n")}
  
  
  MAF = colMeans(G)/2
  MAFsum = colSums(G)
  mafcutoff = 0.01 
  
  weight = weights
  variancematrix = t(G) %*% (w * G) - (t(G) %*% (w * X)) %*%
    temp1 %*% (t(w * X) %*% G)
  
  
  out_kernel=SPA_ER_kernel(G, obj,  u, Cutoff, variancematrix, weight)
  
  zscore.all_1 = out_kernel$zscore.all_0 * weight
  VarS = out_kernel$VarS * weight^2
  pi_1 = obj$pi_1
  G_w = Matrix(t(t(G) * weight), sparse = TRUE)
  G2_adj_n = as.matrix(t(G_w) %*% (w * G_w) - (t(G_w) %*% (w *
                                                             X)) %*% temp1 %*% (t(w * X) %*% G_w))
  rm(G_w)
  gc()
  
  if (method=="SKAT"){r.all=0;r.corr=0;} else {
    if (method=="burden"){r.all=1;r.corr=1;} else{
      if (method=="SKATO"){if (length(r.corr)==0) {    r.all = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);    r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1) } else {r.all=r.corr}
      } else {stop("WARNING: wrong method!\n")}
    }
  }
  
  IDX <- which(r.all >= 0.999)
  if (length(IDX) > 0) {
    r.all[IDX] <- 0.999
  }
  list_myfun = list()
  
  VarS_org = diag(G2_adj_n)
  vars_inf = which(VarS == Inf)
  if (length(vars_inf) > 0) {
    VarS[vars_inf] = 0
    zscore.all_1[vars_inf]=0
    G2_adj_n[vars_inf,]=0
    G2_adj_n[,vars_inf]=0
  }
  if (length(VarS)==1){G2_adj_n = G2_adj_n *VarS/VarS_org} else{
    G2_adj_n = G2_adj_n %*% diag(VarS/VarS_org)
  }
  mu = out_kernel$mu
  g.sum = out_kernel$g.sum
  q.sum = out_kernel$q.sum
  p.value_burden <- SPAtest:::Saddle_Prob(q.sum, mu = mu, g = g.sum,
                                          Cutoff = 2, alpha = 2.5 * 10^-6)$p.value
  v1 = rep(1, dim(G2_adj_n)[1])
  VarQ = t(v1) %*% G2_adj_n %*% v1
  p.m <- dim(G)[2]
  Q_b = p.m^2 * rowMeans(zscore.all_1)^2
  VarQ_2 = Q_b/qchisq(p.value_burden, df = 1, ncp = 0, lower.tail = FALSE,
                      log.p = FALSE)
  if (VarQ_2 == 0) {
    r = 1
  }     else {
    r = VarQ/VarQ_2
  }
  r = min(r, 1)
  #list_myfun$r = r
  
  
  if (dim(G2_adj_n)[2]==1){Phi_temp=as.matrix(G2_adj_n *1/r)} else {Phi_temp=as.matrix(G2_adj_n %*% diag(rep(1/r, dim(G2_adj_n)[2])))}
  
  out = try(SKAT:::Met_SKAT_Get_Pvalue(Score = zscore.all_1,
                                       Phi = Phi_temp,
                                       r.corr = r.all, method = "optimal.adj", Score.Resampling = NULL),
            silent = TRUE)
  if (class(out) != "try-error") {
    list_myfun$p.value = out$p.value
    list_myfun$p.value_each = out$param$p.val.each
  }     else {
    
    
    list_myfun$p.value_each=c(p.value_burden)
    list_myfun$p.value=p.value_burden  ##2*min( list_myfun$p_each)
    cat("Try-error message from SKAT. Only report robust burden.\n")

    
  }
  list_myfun$Q=SKAT:::SKAT_META_Optimal_Get_Q(zscore.all_1, r.corr)$Q.r
  ##list_myfun$p.old = out_kernel$p.old
  list_myfun$p.value_singlevariant = out_kernel$p.new
  
  return(list_myfun)
}



#
# x is either y or SKAT_NULL_Model 
#
SKATBinary_Robust.SSD.OneSet = function(SSD.INFO, SetID, obj, ...){
  
  id1<-which(SSD.INFO$SetInfo$SetID == SetID)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
    stop(MSG)
  }	
  Set_Index<-SSD.INFO$SetInfo$SetIndex[id1]
  
  Z<-Get_Genotypes_SSD(SSD.INFO, Set_Index)
  re<-SKATBinary_Robust(Z, obj, ...)
  
  return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKATBinary_Robust.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ...){
  
  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }	
  SetID<-SSD.INFO$SetInfo$SetID[id1]
  
  
  Z<-Get_Genotypes_SSD(SSD.INFO, SetIndex)
  re<-SKATBinary_Robust(Z, obj, ...)
  return(re)
}




#
# Only SKAT_Null_Model obj can be used
#
SKAT_CommonRare_Robust.SSD.All = function(SSD.INFO, obj, ...){
  
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.Pvalue.Resampling<-NULL
  OUT.Q<-rep(NA,N.Set)
  
  OUT.nRare<-rep(NA,N.Set)
  OUT.nCommon<-rep(NA,N.Set)
  
  Is.Resampling = FALSE
  n.Resampling = 0
  
  if(class(obj) == "SKAT_NULL_Model"){
    if(obj$n.Resampling > 0){
      Is.Resampling = TRUE
      n.Resampling = obj$n.Resampling
      
      OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
    }
  }
  pb <- txtProgressBar(min=0, max=N.Set, style=3)
  for(i in 1:N.Set){
    Is.Error<-TRUE
    try1<-try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
    if(class(try1) != "try-error"){
      Z<-try1
      Is.Error<-FALSE
      
      
    } else {
      err.msg<-geterrmessage()
      msg<-sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
    }
    
    if(!Is.Error){
      Is.Error<-TRUE
      try2<-try(SKATBinary_Robust(Z, obj, ...),silent = TRUE)
      
      if(class(try2) != "try-error"){
        re<-try2
        Is.Error<-FALSE
      } else {
        
        err.msg<-geterrmessage()
        msg<-sprintf("Error to run SKATBinary_Robust for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
        warning(msg,call.=FALSE)
      }
    }
    
    if(!Is.Error){
      
      OUT.Pvalue[i]<-re$p.value
      OUT.Marker[i]<-re$param$n.marker
      OUT.Marker.Test[i]<-re$param$n.marker.test
      OUT.nRare[i]<-re$n.rare
      OUT.nCommon[i]<-re$n.common
      OUT.Q[i]<-re$Q
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
  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, Q=OUT.Q
                      , N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test, N.Marker.Rare=OUT.nRare, N.Marker.Common=OUT.nCommon)
  re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
  class(re)<-"SKAT_SSD_ALL"
  
  return(re)	
}



colMax <- function(data) apply(data,2, max, na.rm = TRUE)


SKATBinary_Robust<-function(Z, obj, kernel = "linear.weighted", method="SKATO"
                            , r.corr=NULL, weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weights = NULL
                            , CommonRare_Cutoff=NULL, impute.method = "bestguess"
                            ,is_dosage = FALSE, missing_cutoff=0.15, max_maf=1
                            , estimate_MAF=1){
  
  
  SetID1=NULL
  # This function only can be used for SNPs
  is_check_genotype=TRUE
  
  if(class(obj) == "SKAT_NULL_Model_ADJ"){
    
    obj.res=obj$re1
    
  } else if(class(obj) == "SKAT_NULL_Model"){
    
    obj.res=obj
    
  } else {
    stop("Wrong obj!")
  }
  
  
  
  if(is.matrix(Z) != TRUE  && class(Z)!="dgCMatrix" && class(Z)!="dgeMatrix"){
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
  
  out<-SKAT:::SKAT_MAIN_Check_Z(Z, obj.res$n.all, id_include=obj.res$id_include, SetID=SetID1, weights=weights, weights.beta=c(1,1), 
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
  MAF<-SKAT:::Get_MAF(Z)
  maf_temp = which(MAF <= mafcutoff)
  
  if (is.null(weights)){
    weight=rep(0,length(MAF))
    if (length(maf_temp) > 0 & length(maf_temp) < length(MAF)) {
      weight[maf_temp] = Beta_Weight(MAF[maf_temp], weights.beta.rare)
      weight[-maf_temp] =Beta_Weight(MAF[-maf_temp], weights.beta.common)
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
  

  
  if(class(obj) == "SKAT_NULL_Model_ADJ"){
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
      Z=Z[,-list_tiny];weights=weigths[-list_tiny];
      }else { stop("all genotypes are close to 0!")
    }
  }
  
  if (max(Z)>2 | min(Z)<0) {stop("Z is out of bounds[0,2]!")}
  
  m.test<-ncol(Z)
  MAF<-SKAT:::Get_MAF(Z)  
    
  id.rare<-intersect(which(MAF < CommonRare_Cutoff), which(MAF > 0))
  id.common<-intersect(which(MAF >= CommonRare_Cutoff), which(MAF > 0))
  
  n.rare= length(id.rare)
  n.common= length(id.common)
  if (length(id.rare)==0){mac.rare=0}else {mac.rare=sum(Z[,id.rare])}
  if (length(id.common)==0){mac.common=0}else {mac.common=sum(Z[,id.common])}
  
  
  is.run<-FALSE
  

  
  re<-SKATBinary_spa(G=Z,obj=obj.res,weights = weights, method=method, r.corr=r.corr)
  
  
  is.run=TRUE
  
  re$param$n.marker<-m.org
  re$param$n.marker.test<-m.test
  re$param$n.marker.name<-colnames(Z)
  re$param$rho=r.corr
  re$param$minp=min( re$p.value_each)
  re$param$rho_est=r.corr[which.min(re$p.value_each)]
  re$n.rare = n.rare
  re$mac.rare=mac.rare
  re$n.common = n.common
  re$mac.common=mac.common
  re$Cutoff = CommonRare_Cutoff
  
  return(re)
  
}





  

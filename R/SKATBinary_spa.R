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
	zscore.all_0<-matrix(rep(0, (ncol(G)*res_time)), ncol=ncol(G))
	zscore.all_1<-matrix(rep(0, (ncol(G)*res_time)), ncol=ncol(G))
	VarS=c()
	#p_all=matrix(rep(0, (ncol(G)*res_time)), ncol=ncol(G))
	p.old=c()
	p.new=c()

	g.sum=0
	q.sum=0

		for (jj in 1:ncol(G)){


			n.g<-sum(G[,jj])
			if(n.g/(2*length(G[,jj]))>0.5)
			{
				G[,jj]<-2-G[,jj]
				n.g<-sum(G[,jj])
			}
			NAset<-which(G[,jj]==0)
			G1<-G[,jj]  -  obj$XXVX_inv %*%  (obj$XV %*% G[,jj])
			#q<-(t(G1) %*% (res_all+u)) /sqrt(n.g)
			#g=G1 /sqrt(n.g)
			q<-(t(G1) %*% (res_all+u)) 
			g=G1

			mu.qtemp=u; g.qtemp=g   ###G_tilde*w^(-1/2) 
			mu1 <- sum(mu.qtemp * g.qtemp)
 			var1 <- sum(mu.qtemp * (1 - mu.qtemp) * g.qtemp^2)


    			stat.qtemp<-(q - mu1)^2/var1
    			p_temp1<-pchisq(stat.qtemp, lower.tail = FALSE, df = 1)  
			p.old[jj]=p_temp1


			g.sum=g.sum+g.qtemp*weight[jj]*sqrt(n.g)
			q.sum=q.sum+q*weight[jj]*sqrt(n.g)

			##zscore.all_0[,jj]=qnorm(p_temp1/2, mean = 0, sd =sqrt( variancematrix[jj,jj]),lower.tail = FALSE, log.p = FALSE)*weight[jj]*sign(q-mu1)
			zscore.all_0[,jj]=(q-mu1)*sqrt(n.g)
	
			id1<-which(stat.qtemp > Cutoff^2) 

			qtemp[jj,]=q
	
			if (MAFsum[jj]<10){

				if (length(id1)>0 ){
					G_temp=G[,jj]
					G_temp[which(G_temp<=0.2)]=0

					temp_binary=SKATBinary(as.matrix(G_temp),obj, method.bin="Hybrid")
					p_temp_binary=c(temp_binary$p.value,temp_binary$p.value.resampling)

					p_temp1[id1]=p_temp_binary[id1]						
				}
	
			}	else {
	
				if(length(id1) > 0){
					qtemp_id1=qtemp[jj,id1]
					#id1_order=order(qtemp_id1)
					#qtemp_id1_order=qtemp_id1[id1_order]

					p_temp1_id1=c()
					for (jjk in 1:length(qtemp_id1)){
						p_temp1_id1[jjk]=SPAtest:::Saddle_Prob(qtemp_id1[jjk], mu=mu.qtemp, g=g.qtemp, Cutoff=Cutoff,alpha=5*10^-8)$p.value
						if (p_temp1_id1[jjk]!=0){
							p_temp1[id1[jjk]]=p_temp1_id1[jjk]
						}
					}				


				}
			}
			p.new[jj]=p_temp1
			if (variancematrix[jj,jj]<=0){zscore.all_1[,jj]=0} else{
				zscore.all_1[,jj]=qnorm(p_temp1/2, mean = 0, sd =sqrt( variancematrix[jj,jj]),lower.tail = FALSE, log.p = FALSE)*weight[jj]*sign(q-mu1)

			}
			if (p_temp1>0){
				VarS[jj]= zscore.all_0[,jj]^2/qchisq(p_temp1, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
			} else {
				VarS[jj]= zscore.all_0[,jj]^2/10^6
 			}
		}##for every col of G
		outlist=list();
		outlist$p.old=p.old
		outlist$p.new=p.new
		outlist$zscore.all_0=zscore.all_0
		outlist$VarS=VarS
		outlist$mu=mu.qtemp
		outlist$g.sum=g.sum
		outlist$q.sum=q.sum
		outlist$zscore.all_1=zscore.all_1
		return(outlist) ;
}





SKATBinary_spa<-function(G, obj, Cutoff ){
		X=obj$X1
		u=obj$mu
		w=obj$pi_1
	
		obj$XV = t(X * w)
		temp1= solve(t(X)%*%(X * w))
		obj$XXVX_inv= X %*% temp1


		for (jj in 1:ncol(G)){

			n.g<-sum(G[,jj])
			if(n.g/(2*length(G[,jj]))>0.5)
			{
				G[,jj]<-2-G[,jj]
				n.g<-sum(G[,jj])
			}
		}

		MAF=colMeans(G)/2
		MAFsum=colSums(G)

		mafcutoff=0.01  ###########1/sqrt(nrow(G_o) * 2)
		weight=rep(0,length(MAF))
		maf_temp=which(MAF<mafcutoff)##rare variants
		if (length(maf_temp)>0 & length(maf_temp)<length(MAF)){
			weight[maf_temp]=Beta_Weight(MAF[maf_temp],c(1,25))
			weight[-maf_temp]=Beta_Weight(MAF[-maf_temp],c(0.5,0.5))
			flag=1
		} else {
			if (length(maf_temp)==0){weight=Beta_Weight(MAF,c(0.5,0.5));flag=2} ###flag 1 means both; 2 only common; 3 only rare;
			if (length(maf_temp)==length(MAF)){weight=Beta_Weight(MAF,c(1,25));flag=3}

		}

		variancematrix=t(G)%*%(w*G)-(t(G)%*%(w*X))%*%temp1%*%(t(w*X)%*%G)

		res_all=cbind(obj$res,obj$res.out) 
		res_temp2=(res_all+u)*w^(-1/2)
		qtemp=t(G)%*%res_temp2-t(G)%*%  (w^(0.5)*X )%*%temp1%*% (t(X) %*% (res_temp2* w^(0.5)))


		out_kernel=SPA_ER_kernel(G, obj, res_all, u, Cutoff, MAFsum,qtemp, resampling_time=0,variancematrix,weight);
		zscore.all_1=out_kernel$zscore.all_0* weight
		VarS=out_kernel$VarS*weight^2

		pi_1=obj$pi_1
		#G_adj=as.matrix(G*pi_1^(1/2)- X %*%temp1%*% (t(X) %*% (G *pi_1^(1/2) * pi_1)) )
		#G2_adj=t(G_adj)%*%G_adj
	
		G_w=Matrix(t(t(G)*weight),sparse=TRUE)
		G2_adj_n=as.matrix(t(G_w)%*%(w*G_w)-(t(G_w)%*%(w*X))%*%temp1%*%(t(w*X)%*%G_w))
		rm(G_w)
		gc()

		r.all = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
		r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)


		IDX<-which(r.all >= 0.999)
		if(length(IDX) > 0){
			r.all[IDX]<-0.999	
		}

		out=SKAT::: Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(G2_adj_n), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL)
		list_myfun=list();
		list_myfun$p_skato_old=out$p.value
		list_myfun$p_each_old=out$param$p.val.each



		VarS_org=diag(G2_adj_n)
		
		vars_inf=which(VarS==Inf)
		if (length(vars_inf)>0){
			VarS[vars_inf]=VarS_org[vars_inf]
		}
		
		G2_adj_n=G2_adj_n%*%diag(VarS/VarS_org)
		
		
		mu =out_kernel$mu
		g.sum =out_kernel$g.sum
		q.sum=out_kernel$q.sum
		p.value_burden<-SPAtest:::Saddle_Prob(q.sum , mu=mu, g=g.sum, Cutoff=2,alpha=2.5*10^-6)$p.value


		v1=rep(1,dim(G2_adj_n)[1])
		VarQ=t(v1)%*%G2_adj_n %*%v1


		p.m<-dim(G)[2]

		Q_b=p.m^2 * rowMeans(zscore.all_1)^2

		VarQ_2=Q_b/qchisq(p.value_burden, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

		if (VarQ_2== 0) {r=1} else {r=VarQ/VarQ_2}
		r=min(r,1)
		
		list_myfun$r=r;
		out=try(SKAT::: Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(G2_adj_n), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL),silent=TRUE)		
		if (class(out)!="try-error") {
		
			list_myfun$p_skato=out$p.value
			list_myfun$p_each=out$param$p.val.each
		} else {
			list_myfun$p_skato=NA
			list_myfun$p_each=rep(NA,7)
			
		}

		out=try(SKAT::: Met_SKAT_Get_Pvalue(Score=zscore.all_1, Phi=as.matrix(G2_adj_n%*%diag(rep(1/r,dim(G2_adj_n)[2]))), r.corr=r.all, method="optimal.adj",Score.Resampling=NULL),silent=TRUE)	
		if (class(out)!="try-error") {
			list_myfun$p_skato_2=out$p.value
			list_myfun$p_each_2=out$param$p.val.each
		} else {
			list_myfun$p_skato_2=NA
			list_myfun$p_each_2=rep(NA,7)

		}
		

		
		if (flag==2) {list_myfun$rare_n=0; list_myfun$common_n=length(MAF); list_myfun$rare_mac=0;list_myfun$common_mac=sum(G);}
		if (flag==1){list_myfun$rare_n=length(maf_temp); list_myfun$common_n=length(MAF)-length(maf_temp); list_myfun$rare_mac=sum(G[,maf_temp]);list_myfun$common_mac=sum(G[,-maf_temp]);}
		if (flag==3) {list_myfun$rare_n=length(MAF); list_myfun$common_n=0; list_myfun$rare_mac=sum(G);list_myfun$common_mac=0;}
		list_myfun$p.old=out_kernel$p.old
		list_myfun$p.new=out_kernel$p.new
		return (list_myfun);
}

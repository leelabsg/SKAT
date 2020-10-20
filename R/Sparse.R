############################################################
##This function calculates the variance matrix of weighted score statistics
##G is a spase genotype matrix
## obj is the null model from SKAT_Null_Model
## weight is the weight vector of snps 
## temp1 is  
###############################################################
Sparse_Phi<-function(G,obj, weight,temp1=NULL){
  if (temp1==NULL){
    temp1=solve(t(obj$X1)%*%(obj$pi_1*obj$X1))
  }
  G_w=Matrix(t(t(G)*weight),sparse=TRUE)
  G2_adj_n=as.matrix(t(G_w)%*%(obj$pi_1*G_w)-(t(G_w)%*%(obj$pi_1*obj$X1))%*%temp1%*%(t(obj$pi_1*obj$X1)%*%G_w))
  return( G2_adj_n)
}


#######################################################
##Get Sparse Matrix From SSD file
##Can save much memory compared with Get_Genotype_SSD
##################################################################

#Changed by SLEE 12/23, change the name from Get_Genotypes_SSD_New to Get_Genotypes_SSD_Sparse
Get_Genotypes_SSD_Sparse<-function(SSD_INFO, Set_Index){

	SNP_ID_SIZE=1024 # it should be the same as SNP_ID_SIZE_MAX in error_messages.h 
	
	Is_MakeFile=0
  	if(get("SSD_FILE_OPEN.isOpen", envir=SSD.env) == 0){
		stop("SSD file is not opened. Please open it first!")
	}

	id1<-which(SSD_INFO$SetInfo$SetIndex == Set_Index)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!\n", Set_Index)
		stop(MSG)
	}	
	Set_Index<-SSD_INFO$SetInfo$SetIndex[id1]

	err_code<-0
	N.SNP_total<-SSD_INFO$SetInfo$SetSize[id1]
	N.Sample<-SSD_INFO$nSample
	
	N.SNP_left=N.SNP_total 
	Z.out.t=NULL
	Z.out=NULL
	Pos=SSD_INFO$SetInfo$Offset[id1]
	
	flag=FALSE
	i=1
	while (flag==FALSE){
		if (N.SNP_left>10 ){
			N.SNP=10
			N.SNP_left=N.SNP_left-10
		} else {
			flag=TRUE
			N.SNP=N.SNP_left	
		}
		size<-N.SNP * N.Sample
		Z<-rep(9,size)

		
		SNPID=raw(N.SNP* SNP_ID_SIZE)
			
		

		temp<-.C("R_Get_Genotypes_withID_new",as.integer(Set_Index),as.integer(Z), SNPID, as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), as.integer(Pos),as.integer(N.SNP),PACKAGE="SKAT")

	
			
		error_code<-temp[[6]]
		Print_Error_SSD(error_code)
		
		SNPID.m<-matrix(temp[[3]], byrow=TRUE, nrow=N.SNP)
		SNPID.c<-apply(SNPID.m, 1, rawToChar)

		Z.out.t<-Matrix(temp[[2]],byrow=TRUE, nrow=N.SNP,sparse=TRUE)
		rownames(Z.out.t)<-SNPID.c
			

		Pos=temp[[7]]
		rm(temp)
		gc()
		if (i==1){Z.out<-Z.out.t} else {Z.out=Matrix(rbind(Z.out,Z.out.t), sparse=TRUE)}
		i=i+1;		
			
	}	
	
	
  	Z.out<-t(Z.out)
	return(Z.out)
}


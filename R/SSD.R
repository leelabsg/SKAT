SSD.env <- new.env()


Print_Error_SSD<-function(code){

	if(code == 0 ){
		return(0);
	} else if(code == 1){
		stop("Error Can't open BIM file")
	} else if(code == 2){
		stop("Error Can't open FAM file")
	} else if(code == 3){
		stop("Error Can't open BED file")
	} else if(code == 4){
		stop("Error Can't open SETID file")
	} else if(code == 5){
		stop("Error Can't write SSD file")
	} else if(code == 6){
		stop("Error Can't read SSD file")
	} else if(code == 7){
		stop("Error Can't write INFO file")
	} else if(code == 8){
		stop("Error Can't read INFO file")
	} else if(code == 9){
		stop("Error Can't write INFO file")
	} else if(code == 13){
		stop("Error Wrong SNP or Individual sizes")
	} else if(code == 14){
		stop("Error SetID not found")
	} else {
		MSG<-sprintf("Error [%d]\n",code)
		stop(MSG)
	}
	return(1)
}

Read_SNP_WeightFile<-function(FileName){

	#FileName<-"./Example1_W1.txt"
	Check_File_Exists(FileName)
	dat<-read.table(FileName, header=FALSE, stringsAsFactors=FALSE)
	hashset<-new.env()
	
	n.SNP<-length(dat[,1])
	for(i in 1:n.SNP){
		val1<-dat[i,1]
		val2<-dat[i,2]
		hashset[[val1]]<-val2
	}
	
	obj<-list(hashset=hashset, nSNP=n.SNP)
	class(obj)<-"SNPWeight"
	return(obj)
}


Check_File_Exists<-function(FileName){
	
	if(!file.exists(FileName)){
		Msg<-sprintf("File %s does not exist\n",FileName)
		stop(Msg)
	}

}

Check_Duplicate<-function(SSD.Info){

	cat("\nCheck duplicated SNPs in each SNP set\n")
	colnames(SSD.Info)<-c("SetID","SNPID")

	ID.u<-unique(SSD.Info[,1])
	ID.n<-length(ID.u)
	ID.d<-data.frame(SetID=ID.u, idx=1:ID.n)
	
	temp<-merge(ID.d, SSD.Info, by.x="SetID", by.y="SetID")
	temp1<-temp[order(temp$idx),]
	
	idx.int.end<-findInterval(1:ID.n,temp1$idx)
	idx.int.start<-c(1,idx.int.end[-ID.n]+1)

	for(i in 1:ID.n){
	
		idx.start<-idx.int.start[i]
		idx.end<-idx.int.end[i]
		
		SetID1<-temp1$SetID[idx.start]
		SetID2<-temp1$SetID[idx.end]
		
		if(SetID1 != SetID2){
			msg<-sprintf("Check Duplicate, SetIDs don't match [%s] [%s]\n", SetID1, SetID2)
			stop(msg)
		}
		if(length(unique(temp1$SNPID[idx.start:idx.end])) != length(temp1$SNPID[idx.start:idx.end])){
			msg<-sprintf("%s has duplicates! Check following variants : ", SetID1);
			tbl<-table(temp1$SNPID[idx.start:idx.end])
			idx<-which(tbl > 1)
			
			msg<-paste(msg, names(tbl[idx]))
			stop(msg)
		}
	}
	
	cat("No duplicate\n")
}


Check_ID_Length<-function(FileName){
	

	SSD.Info<-try(read.table(FileName, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SSD.Info)=="try-error"){
		stop("Error in SetID file!") 
	}

	# increase the limit ot 50
	n1<-length(which(nchar(SSD.Info[,1]) > 1024))
	n2<-length(which(nchar(SSD.Info[,2]) > 1024))

	if(n1 > 0){
		stop("Some SetIDs have more than 1024 characters!") 
	}
	if(n2 > 0){
		stop("Some SNP_IDs have more than 1024 characters!") 
	}	

	nSets<-length(unique(SSD.Info[,1]))
	
	Check_Duplicate(SSD.Info)
	
	return(list(nSets=nSets))
}





Print_File_Info<-function(INFO){
	
	MSG<-sprintf("%d Samples, %d Sets, %d Total SNPs\n",INFO$nSample, INFO$nSets, INFO$nSNPs)
	cat(MSG)

}

Read_File_Info<-function(File.Info){

	Check_File_Exists(File.Info)
	
	info1<-read.table(File.Info, header=FALSE, nrows= 6, sep='\t')
	info2<-read.table(File.Info, header=FALSE, skip=7, sep='\t', stringsAsFactors=FALSE,comment.char = "")

	INFO<-list()
	INFO$WindowSize<-info1[1,1]
	INFO$MAFConvert<-info1[2,1]
	INFO$nSNPs<-info1[3,1]	
	INFO$nSample<-info1[4,1]
	INFO$nDecodeSize<-info1[5,1]
	INFO$nSets<-info1[6,1]
	INFO$SetInfo<-data.frame(SetIndex= info2[,1], SetID = as.character(info2[,3])
	, SetSize = info2[,4], Offset = info2[,2],stringsAsFactors=FALSE)

	Print_File_Info(INFO)
	
	return(INFO)


}


Read_File_Info_Head<-function(File.Info, Is.Print=TRUE){

	Check_File_Exists(File.Info)
	
	info1<-read.table(File.Info, header=FALSE, nrows= 6, sep='\t')


	INFO<-list()
	INFO$WindowSize<-info1[1,1]
	INFO$MAFConvert<-info1[2,1]
	INFO$nSNPs<-info1[3,1]	
	INFO$nSample<-info1[4,1]
	INFO$nDecodeSize<-info1[5,1]
	INFO$nSets<-info1[6,1]
	
	if(Is.Print == TRUE){
		Print_File_Info(INFO)
	} 

	return(list(nSets=INFO$nSets))

}

Check_SetID_LOG<-function(File.SSD){

	File.LOG<-sprintf("%s_LOG.txt",File.SSD)
	if(!file.exists(File.LOG)){
		MSG<-sprintf("Error %s file did not generated!\n",File.LOG)
		stop(MSG)
	}
	
	temp<-read.delim(File.LOG, header = FALSE, sep = "\n",comment.char = "")
	N.Miss<-dim(temp)[1] -1

	if(N.Miss > 0){
		MSG1<-sprintf("Warning: %d SNPs in the SetID file were not found in Bim file!\n",N.Miss)
		MSG2<-sprintf("Please check %s file!\n",File.LOG)
		
		cat(MSG1)
		cat(MSG2)
	}

}


Create_Temporaly_SetID<-function(FileName){
	

	FileName.out<-sprintf("%s.temp",FileName)
	SSD.Info<-try(read.table(FileName, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SSD.Info)=="try-error"){
		stop("Error in SepID file!") 
	}

	ord1<-order(SSD.Info[,1])
	SSD.Info1<-SSD.Info[ord1,]

	write.table(SSD.Info1, file=FileName.out, quote=FALSE, row.names=FALSE, col.names=FALSE)
	
	return(FileName.out)
}



##################################################################
#
#	Generate SSD Files

Generate_SSD_SetID_Work<-function(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype){


	File.Bed<-normalizePath(File.Bed ,mustWork =FALSE)
	File.Bim<-normalizePath(File.Bim ,mustWork =FALSE)
	File.Fam<-normalizePath(File.Fam ,mustWork =FALSE)
	File.SetID<-normalizePath(File.SetID ,mustWork =FALSE)
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)


	Check_File_Exists(File.Bed)
	Check_File_Exists(File.Bim)
	Check_File_Exists(File.Fam)
	Check_File_Exists(File.SetID)

	nSet1 = Check_ID_Length(File.SetID)$nSet
	
	err_code<-0
	MAFConvert=0
	if(Is.FlipGenotype){
		MAFConvert=1
	}


	temp<-.C("R_Generate_MWA_SetID_File"
	, as.character(File.Bed), as.character(File.Bim), as.character(File.Fam)
	, as.character(File.SetID), as.character(File.SSD), as.character(File.Info)
	, as.integer(MAFConvert), as.integer(err_code))
	
	
	error_code<-temp[[8]]
	Print_Error_SSD(error_code)

        # Check SetID_LOG file 
	Kill_SSD_SetID()	
	Check_SetID_LOG(File.SSD)

	nSet2 = Read_File_Info_Head(File.Info, Is.Print=FALSE)$nSet

	re=1;
	if(nSet1 < nSet2){
		re=-1;	# re=-1 not matching numbers
	}

	return(re)	
}


Generate_SSD_SetID<-function(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE){

	

	re = Generate_SSD_SetID_Work(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype)

	if(re == -1){
		MSG1<-sprintf("Warning: SSD file has more SNP sets then SetID file.It happens when SNPs in sets are not contiguous!\n")
		MSG2<-sprintf("SKAT generates a temporary SetID file with contiguous SNP sets. However the order of SNP sets will be different from the order of SNP sets in the original SetID file! \n\n")
		cat(MSG1)
		cat(MSG2)

		File.SetID.temp<-Create_Temporaly_SetID(File.SetID)
		re = Generate_SSD_SetID_Work(File.Bed, File.Bim, File.Fam, File.SetID.temp, File.SSD, File.Info, Is.FlipGenotype)
		if(re == -1){
			stop("Keep having more SetIDs in the SSD file!")
		}
	}	

	Read_File_Info_Head(File.Info, Is.Print=TRUE)
	print("SSD and Info files are created!")

}

Kill_SSD_SetID<-function(){
	temp<-.C("R_Kill_MWA_SetID_File")
}


##################################################################
#
#	Open and Close the SSD Files


assign("SSD_FILE_OPEN.isOpen", 0, envir=SSD.env)
assign("SSD_FILE_OPEN.FileName","", envir=SSD.env)

Close_SSD<-function(){

	if(get("SSD_FILE_OPEN.isOpen", envir=SSD.env) == 1){
		temp<-.C("R_Close_MWA")
		Msg<-sprintf("Close the opened SSD file: %s\n"
		,get("SSD_FILE_OPEN.FileName", envir=SSD.env));
		cat(Msg)
		assign("SSD_FILE_OPEN.isOpen", 0, envir=SSD.env);
	} else{
		Msg<-sprintf("No opened SSD file!\n");
		cat(Msg)		
	}
}

Open_SSD<-function(File.SSD, File.Info){

	err_code<-0
	File.SSD<-normalizePath(File.SSD ,mustWork =FALSE)
	File.Info<-normalizePath(File.Info ,mustWork =FALSE)

	Check_File_Exists(File.SSD)
	Check_File_Exists(File.Info)

	if(get("SSD_FILE_OPEN.isOpen", envir=SSD.env) == 1){
		Close_SSD();
	}

	# Read Info File
	INFO<-Read_File_Info(File.Info)

	# Read SSD File
	temp<-.C("R_Open_MWA", as.character(File.SSD), as.character(File.Info)
	, as.integer(err_code))

	error_code<-temp[[3]]
	Print_Error_SSD(error_code)


	Msg<-sprintf("Open the SSD file\n");
	cat(Msg)

	#SSD_FILE_OPEN.isOpen<<-1
	#SSD_FILE_OPEN.FileName<<-File.SSD

	assign("SSD_FILE_OPEN.isOpen", 1, envir=SSD.env)
	assign("SSD_FILE_OPEN.FileName",File.SSD, envir=SSD.env)


	return(INFO)
	
}

#######################################################



##################################################################
#
#	Get Genotype Matrix

Get_Genotypes_SSD<-function(SSD_INFO, Set_Index, is_ID = FALSE){

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
	N.SNP<-SSD_INFO$SetInfo$SetSize[id1]
	N.Sample<-SSD_INFO$nSample
	size<-N.SNP * N.Sample

	Z<-rep(9,size)

	if(!is_ID){
		temp<-.C("R_Get_Genotypes",as.integer(Set_Index),as.integer(Z),as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="SKAT")

		error_code<-temp[[5]]
		Print_Error_SSD(error_code)
		
		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		
	} else {
		SNPID=raw(N.SNP* SNP_ID_SIZE)
	
		temp<-.C("R_Get_Genotypes_withID",as.integer(Set_Index),as.integer(Z), SNPID, as.integer(size)
		,as.integer(Is_MakeFile), as.integer(err_code), PACKAGE="SKAT")

		error_code<-temp[[6]]
		Print_Error_SSD(error_code)
		
		SNPID.m<-matrix(temp[[3]], byrow=TRUE, nrow=N.SNP)
		SNPID.c<-apply(SNPID.m, 1, rawToChar)

		Z.out.t<-matrix(temp[[2]],byrow=TRUE, nrow=N.SNP)
		rownames(Z.out.t)<-SNPID.c
		
		#SNPID.c1<<-SNPID.c
		#SNPID.1<<-SNPID	
	}
	
	return(t(Z.out.t))
}


##################################################################
#
#	Read one SNP from plink BED file 

Open_Plink<-function(File.Bed, File.Bim, File.Fam, Is_ReadBim=TRUE){
  
  File.Bed<-normalizePath(File.Bed ,mustWork =FALSE)
  File.Bim<-normalizePath(File.Bim ,mustWork =FALSE)
  File.Fam<-normalizePath(File.Fam ,mustWork =FALSE)
 
  Check_File_Exists(File.Bed)
  Check_File_Exists(File.Bim)
  Check_File_Exists(File.Fam)
 
  err_code<-0

  temp<-.C("R_Open_Plink_BED"
           , as.character(File.Bed), as.character(File.Bim), as.character(File.Fam),
           as.integer(err_code))
  
  
  error_code<-temp[[4]]
  Print_Error_SSD(error_code)

}







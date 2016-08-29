 \name{Read_Plink_FAM}
 \alias{Read_Plink_FAM}
 \alias{Read_Plink_FAM_Cov}
 \title{Read Plink FAM and covariates files}
 \description{
     Read Plink FAM and covariates files.
 }
 \usage{
	Read_Plink_FAM(Filename, Is.binary=TRUE, flag1=0)
	Read_Plink_FAM_Cov(Filename, File_Cov, Is.binary=TRUE, flag1=0, cov_header=TRUE)

 }
\arguments{
      \item{Filename}{input file name of plink FAM file}
      \item{Is.binary}{if TRUE, the phenotype is binary. If phenotype is continuous, it should be FALSE}
      \item{flag1}{0 represents the default coding of unaffected/affected (1/2) (default=0), and 1 represents 0/1 coding. 
      flag1=1 is the same as --1 flag in plink. Please see the plink manual. }      
	  \item{File_Cov}{an input file name of plink covariates file. The first two columns of this file should be FID and IID.}      
	  \item{cov_header}{a logical value indicating whether the covariate file contains a header row (default=TRUE)}      
	
}
\value{
	A dataframe of Family ID (FID), Individual ID (IID), Paternal ID (PID), Maternal ID(MID), Sex, and Phenotype. 
  	If Read_Plink_FAM_Cov is used with a covariate file, the dataframe has covariates from the 7th column.
}



\author{Seunggeun Lee}






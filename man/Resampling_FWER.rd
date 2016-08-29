 \name{Resampling_FWER}
 \alias{Resampling_FWER}
 \alias{Resampling_FWER_1}
 \title{Obtain significant SNP sets after controlling family wise error rate (FWER)}
 \description{

      Obtain significant SNP sets after controlling for family wise error rate (FWER) 
      using resampled residuals. To use it, SKAT_Null_Model or SKAT_Null_Model_MomentAdjust should have n.Resampling > 0. 
 }
 \usage{

	Resampling_FWER(obj,FWER=0.05)

	Resampling_FWER_1(P.value, P.value.Resampling, FWER=0.05)

 }
\arguments{
      \item{obj}{object returned from SKAT.SSD.All function.}
      \item{P.value}{a vector of SKAT p-values. If 100 genes were tested, this vector should have 100 p-values.}
      \item{P.value.Resampling}{a matrix of p-values of the resampled residuals. 
      Each row represents each gene/snp set, and each column represents resampling set. 
      For example, if you have 100 genes, and conducted resampling 1000 times ( ex.n.Resampling=1000 in SKAT_Null_Model), then it should be a 100 x 1000 matrix.}  
      \item{FWER}{a numeric value of FWER rate to control (default=0.05)}
}
\value{
	\item{results}{If you use the returned object from SKAT.SSD.all function, 
	it is a sub-table of significant snp sets of the result table in the obj. 
	If you use P.value and P.value.Resampling, it is a vector of significant p-values. 
	If there is no significant snp set, it is NULL. }
	\item{n}{a numeric value of the number of significant snp sets.}
  	\item{ID}{a vector of indexes of significant snp sets.}
}

\author{Seunggeun Lee}


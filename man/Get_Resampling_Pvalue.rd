 \name{Get_Resampling_Pvalue}
 \alias{Get_Resampling_Pvalue}
 \alias{Get_Resampling_Pvalue_1}
 \title{Compute a resampling p-value}
 \description{
     Compute a resampling p-value using resampled residuals. 
     To use it, SKAT_Null_Model or SKAT_Null_Model_MomentAdjust should have n.Resampling > 0. 
 }
 \usage{
	Get_Resampling_Pvalue(obj)

	Get_Resampling_Pvalue_1(p.value, p.value.resampling)


 }
\arguments{
      \item{obj}{SKAT outcome object.}
      \item{p.value}{a numeric value of the SKAT p-value.}
      \item{p.value.resampling}{a vector of p-values from the resampled residuals.}  
}
\value{
	\item{p.value}{the resampling p-value. It is computed as (n1 +1)/(n+1), where n is the number of resampling (n.Resampling in SKAT_Null_Model or SKAT_Null_Model_MomentAdjust), 
	and n1 is the number of resampled residual p-values smaller than the original sample p-value. }
  	\item{is_smaller}{a logical value indicates whether the resampling p-value should be smaller. 
  	If n1=0, then it is TRUE, otherwise it is FALSE. }
  	
}
\details{
	See SKAT_Null_Model
}

\author{Seunggeun Lee}


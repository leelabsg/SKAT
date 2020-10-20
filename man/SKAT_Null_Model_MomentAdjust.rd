 \name{SKAT_Null_Model_MomentAdjust}
 \alias{SKAT_Null_Model_MomentAdjust}
 \title{Get parameters and residuals from the NULL model for small sample adjustment}
 \description{
     Compute model parameters and residuals for SKAT with adjusting small sample moments when the trait is binary. 
     You can also obtain resampled residuals that can be used to compute resampling p-value or to control family-wise error rate.
 }
 \usage{

SKAT_Null_Model_MomentAdjust(formula, data=NULL, n.Resampling=0,
type.Resampling="bootstrap", is_kurtosis_adj=TRUE, n.Resampling.kurtosis=10000)

 }
\arguments{
      \item{formula}{object of class ``formula'': a symbolic description of the NULL model to be fitted.}
      \item{data}{optional data frame containing the variables in the model (default=NULL).  If it is NULL, the variables are taken from 'environment(formula)'}
      \item{n.Resampling}{a numeric value of the number of resampling (default=0). If you don't want resampling, please set n.Resampling=0. }
      \item{type.Resampling}{ resampling methods (default="bootstrap"). see details.}
      \item{is_kurtosis_adj}{ If TRUE, the kurtosis adjustment will be applied. The small sample kurtosis will be estimated using the resampled phenotypes.}   
      \item{n.Resampling.kurtosis}{ a numeric value of the number of resampling for kurtosis estimation (default=10000). If is_kurtosis_ad=FALSE, it will be ignored. }     
}
\value{
	This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variants and outcome phenotypes. 
	After obtaining it, please use SKAT function to conduct association tests.

}
\details{

When the trait is binary, the SKAT can produce conservative results when the sample size is small. 
To address this, we developed a small sample adjustment method, which adjust asymptotic null distribution by estimating small sample variance and kurtosis. 
The small smaple variance is estimated analytically, and the small sample kurtosis is estimated using the resampling approach.

There are 2 different methods to get resampled residuals.
"bootstrap" conducts the parametric bootstrap to resample residuals under the NULL model with considering covariates. 
"bootstrap.fast" (only for binary traits) is a fast implementation of "bootstrap".
If there is no covariate, "bootstrap" is equivalent to the permutation method.

Since the kurtosis is estimated using random samples, SKAT with the kurtosis-based small sample adjustment 
can yield slightly different p-values for each run. If you want to reproduce p-values, please set a seed number 
using set.seed function in R.  
                                                                  
}


\author{Seunggeun Lee}

\examples{


data(SKAT.example)
attach(SKAT.example)

#############################################################
#	Compute the P-value of SKAT 

IDX<-c(1:100,1001:1100)

# binary trait
obj<-SKAT_Null_Model_MomentAdjust(y.b[IDX] ~ X[IDX,])
SKAT(Z[IDX,], obj, kernel = "linear.weighted")$p.value


}



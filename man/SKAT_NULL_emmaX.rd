 \name{SKAT_NULL_emmaX}
 \alias{SKAT_NULL_emmaX}
 \title{Get parameters and residuals from the null model with incorporating the kinship structure}
 \description{
     Compute model parameters and residuals for SKAT with incorporating the kinship structure. 
 }
 \usage{

SKAT_NULL_emmaX (formula, data=NULL, K=NULL, 
Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE) 
  
 }
\arguments{
      \item{formula}{an object of class ``formula'': a symbolic description of the NULL model to be fitted.}
      \item{data}{an optional data frame containing the variables in the model (default=NULL).  If it is NULL, the variables are taken from 'environment(formula)'}
      \item{K}{a kinship matrix. If K=NULL, the function reads Kin.File to get a kinship matrix.}
      \item{Kin.File}{an emmax-kin output file name. If K=NULL, the function reads this file. }
  	  \item{ngrids}{Number of grids to search for the optimal variance component}
      \item{llim}{Lower bound of log ratio of two variance components}
      \item{ulim}{Upper bound of log ratio of two variance components}
      \item{esp}{Tolerance of numerical precision error}
      \item{Is.GetEigenResult}{Return intermediate eigen-decomposition results}
}
\value{
	This function returns an object that has model parameters and residuals of the NULL model of no associations. 
	After obtaining it, use SKAT function to carry out the association test.

}
\details{

The Emma package code was used to implement this function.

Resampling is not implemented. 
                                                                      
}
\examples{


data(SKAT.fam.example)
attach(SKAT.fam.example)


obj<-SKAT_NULL_emmaX(y ~ X, K=K)
SKAT(Z, obj)$p.value

# SKAT-O
SKAT(Z, obj, method="optimal.adj")$p.value	

}



\author{Seunggeun Lee}

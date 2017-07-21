 \name{SKAT}
 \alias{SKAT}
 \alias{SKAT.SSD.OneSet}
 \alias{SKAT.SSD.OneSet_SetIndex}
 \title{SNP-set (Sequence) Kernel Association Test}
 \description{
     Test for association between a set of SNPS/genes and continuous or dichotomous phenotypes using kernel regression framework.      
 }
 \usage{

SKAT(Z, obj, kernel = "linear.weighted", 
  method="davies", weights.beta=c(1,25), weights=NULL, 
  impute.method="fixed", r.corr=0, is_check_genotype=TRUE,
  is_dosage = FALSE, missing_cutoff=0.15 , max_maf=1, estimate_MAF=1)

SKAT.SSD.OneSet(SSD.INFO, SetID, obj, \dots ,obj.SNPWeight=NULL)

SKAT.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, \dots ,obj.SNPWeight=NULL)

 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. 
      Missing genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based imputation. }
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{kernel}{a type of kernel (default= "linear.weighted"). See detail section. }
      \item{method}{a method to compute the p-value (default= "davies"). 
      "davies" represents an exact method that  computes the p-value by inverting the characteristic function of the mixture chisq, 
      "liu" represents an approximation method that matches the first 3 moments, 
      "liu.mod" represents modified "liu" method that matches kurtosis instead of skewness 
      to improve tail probability approximation,"SKATO" and "optimal.adj" represent a SKAT-O based on an unified approach, 
      and "optimal" is an old version of the implementation of SKAT-O. See details.}
      \item{weights.beta}{a numeric vector of parameters for the beta weights for the weighted kernels. 
      If you want to use your own  weights, please use the ``weights'' parameter. It will be ignored if ``weights'' parameter is not null.}
      \item{weights}{a numeric vector of weights for the weighted kernels. 
      It is \eqn{\sqrt{w}} in the SKAT paper. 
      So if you want to use the Madsen and Browning (2009) weight, you should set each element of weights as \eqn{1/ \sqrt{p(1-p)}}, 
      not \eqn{1/ p(1-p)}. When it is NULL, the beta weight with the ``weights.beta'' parameter is used. }
      \item{impute.method}{a method to impute missing genotypes (default= "fixed"). 
      "bestguess" imputes missing genotypes as most likely values (0,1,2), "random" imputes missing genotypes by generating binomial(2,p) random variables (p is the MAF), 
      and "fixed" imputes missing genotypes by assigning the mean genotype values (2p).}
      \item{r.corr}{the \eqn{\rho} parameter for the compound symmetric correlation structure kernels (default= 0). 
      If you give a vector value, SKAT will conduct the optimal test. It will be ignored if method=``optimal'' or method=``optimal.adj''. See details.}
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). 
      If Z has non-SNP data, please set it FALSE, otherwise you will get an error message. 
      If it is FALSE and you use weighted kernels, the weights should be given through the ``weights'' parameter.}
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{max_maf}{a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff). Any SNPs with MAF > cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. If estimate_MAF=1 (default), SKAT uses all samples to estimate MAFs. 
      If estimate_MAF=2, only samples with non-missing phenotypes and covariates are used to estimate MAFs. }
      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD. }
      \item{SetID}{a character value of Set ID. A set ID of each set can be found from SetInfo object in SSD.INFO.}
      \item{SetIndex}{a numeric value of Set index. A set index of each set can be found from SetInfo object in SSD.INFO.}
      \item{\dots}{further arguments to be passed to ``SKAT'' }
      \item{obj.SNPWeight}{ an output object of Read_SNP_WeightFile (default=NULL). 
      If NULL, the beta weight with the ``weights.beta'' parameter will be used.  }
      
}
\value{
	\item{p.value}{p-value of SKAT. }
	\item{p.value.resampling}{p-values from resampled outcomes. You can get it when you use obj from SKAT_Null_Model function with resampling. See the SKAT_Null_Model. }
	\item{p.value.noadj}{p-value of SKAT without the small sample adjustment. It only appears when small sample adjustment is applied.}
	\item{p.value.noadj.resampling}{p-values from resampled outcomes without the small sample adjustment.}
  	\item{pval.zero.msg}{(only when p.value=0) text message that shows how small the p.value is. ex. "Pvalue < 1.000000e-60" when the p.value is smaller than \eqn{10^{-60}} } 
  	\item{Q}{test statistic of SKAT. It has NA when method="optimal.adj" or "optimal".}
	\item{param}{estimated parameters of each method.}   
	\item{param$Is_Converged}{ (only with method="davies") an indicator of the convergence (1=convergence, 0=non-convergence). 
	When 0 (not converged), "liu" method will be used to compute p-values. }  
	\item{param$n.marker}{a number of SNPs in the genotype matrix.}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	
}
\details{


There are 6 types of pre-specified kernels:
"linear", "linear.weighted", "IBS", "IBS.weighted", "quadratic" and "2wayIX".
Among them, "2wayIX" is a product kernel consisting of main effects and SNP-SNP interaction terms. 

If users want to use dosage values instead of genotypes, set is_dosage=TRUE. 
Please keep in mind that plink formatted files (so SSD files) cannot used for dosages. Instead, 
you should make a genotype matrix Z to run SKAT.

The kernel matrix for the weighted linear kernel is
\eqn{K = G W W G }, where G is a genotype matrix and W is a diagonal weight matrix. 
Please note that it is different from the notation we used in the original SKAT paper, 
which was \eqn{K = G W G }. 
The Madsen and Browning (2009) weight is \eqn{w = 1/ \sqrt{p(1-p)}} in the current notation. 
By the previous notation, it is \eqn{w = 1/ p(1-p)}.

If you want to use the SSD file, you need to open it first using Open_SSD,
and then use either SKAT.SSD.OneSet  or SKAT.SSD.OneSet_SetIndex. 
Set index is a numeric value and automatically assigned to each set (from 1).

The r.corr represents a \eqn{\rho} parameter of the unified test, 
\eqn{Q_{\rho} = (1-\rho) Q_S + \rho Q_B}, where \eqn{Q_S} is a SKAT test statistic, 
and \eqn{Q_B} is a weighted burden test statistic. 
Therefore, \eqn{\rho=0} results in the original weighted linear kernel SKAT, 
and \eqn{\rho=1} results in the weighted burden test (default: \eqn{\rho=0}). 
If r.corr is a vector, SKAT-O will be conducted with adaptively seleting \eqn{\rho} from given r.corr values. 
\eqn{\rho} should be a value between 0 and 1. When method="optimal" or method="optimal.adj", the r.corr parameter will be ignored.

We slightly changed the implementation for SKAT-O to improve the estimation of p-values. You can run it by 
using method="optimal.adj" or "SKATO". It uses a grid of eight points \eqn{\rho=(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)} 
for the search of the optimal \eqn{\rho}.   
If you want to use the original implementation of SKAT-O, use method="optimal" that 
carries out SKAT-O with an equal sized grid of 11 points (from 0 to 1).

If the true p.value is very small, you can have p.value=0 due to numerical reasons. In this case, 
please see pval.zero.msg that shows how small it is. 
For example, if the p.value is smaller than \eqn{10^{-60}}, it has "Pvalue < 1.000000e-60".                                           

By default, SKAT uses impute.method="fixed" that imputes missing genotypes as the mean genotype values (2p). 
When variates are very rare and missing rates between cases and controls are highly unbalanced, 
impute.method="fixed" can yield inflated type I error rate. In this case, we recommend to use
impute.method="bestguess", which does not suffer the same problem.


}


\author{Seunggeun Lee, Micheal Wu}

\references{

Lee, S., Emond, M.J., Bamshad, M.J., Barnes, K.C., Rieder, M.J., Nickerson, D.A., 
NHLBI GO Exome Sequencing Project-ESP Lung Project Team, Christiani, D.C., Wurfel, M.M. and Lin, X. (2012) 
Optimal unified approach for rare variant association testing with application to small sample 
case-control whole-exome sequencing studies. \emph{American Journal of Human Genetics}, 91, 224-237.

Lee, S., Wu, M. C., and Lin, X. (2012)  Optimal tests for rare variant effects in sequencing association studies. \emph{Biostatistics}, 13, 762-775.

Wu, M. C.*, Lee, S.*, Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011)  Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). \emph{American Journal of Human Genetics}, 89, 82-93. \\
* contributed equally. 

Wu, M. C., Kraft, P., Epstein, M. P.,Taylor, D., M., Chanock, S. J., Hunter, D., J., and Lin, X. (2010)  Powerful SNP Set Analysis for Case-Control Genome-wide Association Studies. \emph{American Journal of Human Genetics}, 86, 929-942.


Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear
  Combination of chi-2 Random Variables, \emph{ Journal of the Royal
  Statistical Society. Series C }, 29, 323-333.

H. Liu, Y. Tang, H.H. Zhang (2009) A new chi-square approximation
  to the distribution of non-negative definite quadratic forms in
  non-central normal variables, \emph{Computational Statistics and Data Analysis}, 53, 853-856.

Duchesne, P. and Lafaye De Micheaux, P. (2010) Computing the distribution of quadratic forms: Further comparisons between the Liu-Tang-Zhang approximation and exact methods, \emph{Computational Statistics and Data Analysis}, 54, 858-862. 


}

\examples{


data(SKAT.example)
attach(SKAT.example)



#############################################################
#	SKAT with default Beta(1,25) Weights 
#		- without covariates

# continuous trait
obj<-SKAT_Null_Model(y.c ~ 1, out_type="C")
SKAT(Z, obj)$p.value

# dichotomous trait 
obj<-SKAT_Null_Model(y.b ~ 1, out_type="D")
SKAT(Z, obj)$p.value


##################################################
#	SKAT with default Beta(1,25) Weights
#		- with covariates

# continuous trait
obj<-SKAT_Null_Model(y.c ~ X, out_type="C")
SKAT(Z, obj)$p.value

obj.b<-SKAT_Null_Model(y.b ~ X, out_type="D")
SKAT(Z, obj.b)$p.value

##################################################
#	SKAT with default Beta(1,25) Weights
#		- Optimal Test

SKAT(Z, obj, method="optimal.adj")$p.value

# you can get the same p-value by using method="SKATO"
SKAT(Z, obj, method="SKATO")$p.value

#############################################################
#	SKAT with Beta(1,30) Weights

SKAT(Z, obj, weights.beta=c(1,30))$p.value


}



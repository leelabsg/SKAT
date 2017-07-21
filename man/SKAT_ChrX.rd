 \name{SKAT_ChrX}
 \alias{SKAT_ChrX}
 \title{SNP-set (Sequence) Kernel Association Test for  X chromosome variables}
 \description{
     Test for association between a set of SNPS/genes in the X chromosome and continuous or dichotomous outcomes using the kernel machine.      
 }
 \usage{

SKAT_ChrX(Z, obj, is_X.inact =TRUE
, kernel = "linear.weighted", method="davies", weights.beta=c(1,25)
, weights = NULL, impute.method = "fixed", r.corr=0, is_check_genotype=TRUE
, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1, estimate_MAF=1, SetID=NULL)


 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. 
      Missing genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based imputation. }
      \item{obj}{output object of the SKAT_Null_Model_ChrX function. }
      \item{is_X.inact}{an indicator variable for the X-inactivation coding (default=TRUE). 
      Male genotypes are coded as g=(0,2) when it is TRUE, and  g=(0,1) when it is false.}
      \item{kernel}{a type of kernel (default= "linear.weighted"). }
      \item{method}{a method to compute the p-value (default= "davies"). See SKAT page for details.}
      \item{weights.beta}{a numeric vector of parameters of beta weights. See SKAT page for details. }
      \item{weights}{a numeric vector of weights for the weighted kernels. See SKAT page for details. }
      \item{impute.method}{a method to impute missing genotypes (default= "fixed"). 
      "fixed" imputes missing genotypes by assigning the mean genotype value (2p), and "bestguess" uses best guess genotype values. }
      \item{r.corr}{the \eqn{\rho} parameter for the compound symmetric kernel. See  SKAT page for details. }
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). See SKAT page for details.}
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{max_maf}{a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff). Any SNPs with MAF > cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. See SKAT page for details. }
      \item{SetID}{Internal use only. }
     
      
}
\value{
	\item{p.value}{p-value of SKAT. }
	\item{p.value.resampling}{p-values from resampled outcome. You can get it when you use obj from SKAT_Null_Model function with resampling. See the SKAT_Null_Model. }
	\item{p.value.noadj}{p-value of SKAT without the small sample adjustment. It only appears when small sample adjustment is applied.}
	\item{p.value.noadj.resampling}{p-values from resampled outcome without the small sample adjustment. It only appears when small sample adjustment is applied. }
  	\item{pval.zero.msg}{(only when p.value=0) text message that shows how small the p.value is. ex. "Pvalue < 1.000000e-60" when p.value is smaller than \eqn{10^{-60}} } 
  	\item{Q}{the test statistic of SKAT. It has NA when method="optimal.adj" or "optimal".}
	\item{param}{estimated parameters of each method.}   
	\item{param$Is_Converged}{ (only with method="davies") an indicator of the convergence. 
	When 0 (not converged), "liu" method is used to compute p-value. }  
	\item{param$n.marker}{a number of SNPs in the genotype matrix}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	
}
\details{


For details of parameters, please see SKAT page.  
                                   

}


\author{Clement Ma and Seunggeun Lee}


\examples{



data(SKAT.example.ChrX)
attach(SKAT.example.ChrX)

#############################################################
#	Compute the P-value of SKAT 

# binary trait
obj.x<-SKAT_Null_Model_ChrX(y ~ x1 +x2 + Gender, SexVar="Gender", out_type="D")

# SKAT
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted")

# Burden
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted", r.corr=1)

# SKAT-O
SKAT_ChrX(Z, obj.x, kernel = "linear.weighted", method="SKATO")




}



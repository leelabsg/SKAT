 \name{SKATBinary_Robust}
 \alias{SKATBinary_Robust}
 \alias{SKATBinary_Robust.SSD.OneSet}
 \alias{SKATBinary_Robust.SSD.OneSet_SetIndex}
 \title{SNP set test for binary traits with robust region-based methods}
 \description{
     This function computes p-values of robust burden test, SKAT, and SKAT-O for binary traits using SPA and ER.   
 }
 \usage{
 
	SKATBinary_Robust(Z, obj, kernel = "linear.weighted", method="SKATO"
	, weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5), weights = NULL
	, r.corr.rare=0, r.corr.common=0, CommonRare_Cutoff=NULL, impute.method = "bestguess"
	, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=1
	, estimate_MAF=1)

	SKAT_CommonRare.SSD.OneSet(SSD.INFO, SetID, obj, \dots)

	SKAT_CommonRare.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, \dots )

 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, 
      where A is a major allele and a is a minor allele. }
      \item{obj}{output object from SKAT_Null_Model. }
      \item{kernel}{type of kernel (default= "linear.weighted"). The possible choices are "linear" and "linear.weighted".}
      \item{method}{type of gene based test (default= "SKATO"). The possible choices are
      "SKAT", "Burden" and "SKATO", which represents robust SKAT, Burden and SKAT-O tests, respectively. 
      \item{weights.beta.rare}{a numeric vector of parameters of beta weights for rare variants (MAF<=0.01). It is only used for weighted kernels. 
      If you want to use your own  weights, please specify the ``weights'' parameter.}
      \item{weights.beta.common}{a numeric vector of parameters of beta weights for common variants (MAF>0.01). It is only used for weighted kernels. 
      If you want to use your own  weights, please specify the ``weights'' parameter.}
      \item{weights}{a numeric vector of weights for the weighted kernels. See SKAT page for details.}
      \item{r.corr.rare}{the \eqn{\rho} parameter for rare variants (default= 0). \eqn{\rho} =0 and 1 indicate SKAT and Burden test, respectively}
      \item{r.corr.common}{the \eqn{\rho} parameter for common variants (default= 0). \eqn{\rho} =0 and 1 indicate SKAT and Burden test, respectively}
      \item{CommonRare_Cutoff}{MAF cutoff for common vs rare variants (default=NULL). It should be a numeric value between 
      0 and 0.5, or NULL. When it is NULL, \eqn{1/ \sqrt{2 SampleSize }} will be used. }	     	     
      \item{impute.method}{a method to impute missing genotypes (default= "bestguess"). "bestguess" imputes missing genotypes as most likely 
      values (0,1,2), "random" imputes missing genotypes by generating binomial(2,p) random variables (p is the MAF), 
      and "fixed" imputes missing genotypes by assigning the mean genotype value (2p).}
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). See SKAT page for details.}
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{max_maf}{a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff). Any SNPs with MAF > cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. See SKAT page for details. }

      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD. }
      \item{SetID}{a character value of Set ID. You can find a set ID of each set from SetInfo object of SSD.INFO. In SKATBinary_Robust function, this parameter is for the internal use only.}
      \item{SetIndex}{a numeric value of Set index. You can find a set index of each set from SetInfo object of SSD.INFO  }
      \item{\dots}{further arguments to be passed to ``SKATBinary_Robust'' }

	  
}
\value{
	\item{p.value}{p-value. It will be the mid p-value if ER or ER.A are used to compute the p-value.}
	\item{p.value.standard}{(ER and ER.A only) standard p-value.}
	\item{p.value.resampling}{p-values from resampled outcome. You can obtain it when n.Resampling (in SKAT_Null_Model) > 0. See SKAT_Null_Model page. }
	\item{p.value.standard.resampling}{(ER and ER.A only)standard p-values from resampled outcomes.}
	\item{m}{the number of individuals with minor alleles.}
	\item{MAP}{minimum possible p-values. It is available when the method.bin="ER" and m is sufficiently small.}
	\item{MAC}{total minor allele count (MAC).}
	\item{n.total}{(ER only) the number of resampling to be generated to get the p-value. 
	It can be smaller than N.Resampling when the total number of configurations of case-controls among individuals with minor alleles are smaller than
	N.Resampling.}
  	\item{is.accurate}{logical value for the accuracy of the p-value. If it is false, more resampling is needed to accurately estimate the p-value. }
	\item{param$n.marker}{a number of SNPs in the genotype matrix}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	\item{method.bin}{a type of method to be used to compute the p-value.}
}
\details{

This function implements six methods (method.bin) to compute p-values: 1) Efficient resampling (ER);
2) Quantile adjusted moment matching (QA); 3) Moment matching adjustment (MA);
4) No adjustment (UA);  5) Adaptive ER (ER.A); and 6) Hybrid. 
"Hybrid" selects a method based on the total minor allele count (MAC), the number of individuals with minor 
alleles (m), and the degree of case-control imbalance. When method.bin="ER" or "ER.A", SKATBinary compute mid-p-values and minimum achievable 
mid p-values. 

If seednum is not NULL, set.seed(seednum) function is used to specify seeds to get the same p-values 
of ER based methods for different runs. Therefore, please set seednum=NULL, if you do not want to set seeds. 

SKATBinary uses impute.method="bestguess" as a default method for the imputation, which is different from SKAT that uses 
impute.method="fixed" as a default method. We changed it because SKATBinary with impute.method="fixed" can yield false positives
when variates are very rare and missing rates between cases and controls are unbalanced. 
When missing rates between cases and controls are highly unbalanced, SKAT impute.method="fixed" can also yield false positives, 
but it happens less likely. So we did not change the default imputation method in SKAT.

}


\author{Seunggeun Lee}

\references{

Lee, S., Fuchsberger, C., Kim, S., Scott, L. (2015) 
An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies.
\emph{Biostatistics}, in press.

}

\examples{


data(SKATBinary.example)
attach(SKATBinary.example)


obj<-SKAT_Null_Model(y ~ x1 + x2, out_type="D")

# run SKAT (default method) with Hybrid
out = SKATBinary(Z, obj)

# p-value
out$p.value

# MAP
out$MAP

# method used to compute p-value (method.bin)
out$method.bin


#
#	Run burden and SKAT-O with Hybrid

SKATBinary(Z, obj, method="Burden")$p.value
SKATBinary(Z, obj, method="SKATO")$p.value

#
#	Run with SKAT-QA, -MA and -UA

SKATBinary(Z, obj, method.bin="QA")$p.value

SKATBinary(Z, obj, method.bin="MA")$p.value

SKATBinary(Z, obj, method.bin="UA")$p.value

# UA from SKAT function
SKAT(Z, obj)$p.value


#
#	Run with Adaptive ER

out =SKATBinary(Z, obj, method.bin="ER.A")

out$p.value

# the number of total resampling is smaller than 2*10^6 (default value)
out$n.total 

}

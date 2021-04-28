 \name{SKATBinary_Robust}
 \alias{SKATBinary_Robust}
 \alias{SKATBinary_Robust.SSD.OneSet}
 \alias{SKATBinary_Robust.SSD.OneSet_SetIndex}
 \title{SNP set test for binary traits with robust region-based methods}
 \description{
     This function computes p-values of robust burden test, SKAT, and SKAT-O for binary traits using SPA and ER.   
 }
 \usage{
 
	SKATBinary_Robust(Z, obj, kernel = "linear.weighted", method="SKAT"
	, r.corr=NULL, weights.beta=c(1,25), weights = NULL
	, impute.method = "bestguess",  is_check_genotype=TRUE
  , is_dosage = FALSE, missing_cutoff=0.15, max_maf=1
	, estimate_MAF=1)

	SKATBinary_Robust.SSD.OneSet(SSD.INFO
	, SetID, obj, \dots,obj.SNPWeight=NULL)

	SKATBinary_Robust.SSD.OneSet_SetIndex(SSD.INFO
	, SetIndex, obj, \dots ,obj.SNPWeight=NULL)
	

 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, 
      where A is a major allele and a is a minor allele. Now we support both matrix and sparse matrix. }
      \item{obj}{output object from SKAT_Null_Model. }
      \item{kernel}{type of kernel (default= "linear.weighted"). The possible choices are "linear" and "linear.weighted".}
      \item{method}{type of gene based test (default= "SKAT"). The possible choices are
      "SKAT", "Burden" and "SKATO", which represents robust SKAT, Burden and SKAT-O tests, respectively. }
      \item{r.corr}{the \eqn{\rho} parameter for all variants. \eqn{\rho} =0 and 1 indicate SKAT and Burden test, respectively.}
      \item{weights.beta}{a numeric vector of parameters of beta weights. 
      It is only used for weighted kernels. 
      If you want to use your own  weights, please specify the ``weights'' parameter.}
      \item{weights}{a numeric vector of weights for the weighted kernels. See SKAT page for details.}	     	     
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
      \item{obj.SNPWeight}{output object from Read_SNP_WeightFile (default=NULL). 
      If NULL, the beta weight with the ``weights.beta'' parameter will be used.  }
      \item{\dots}{further arguments to be passed to ``SKATBinary_Robust'' }
}
\value{
	\item{p.value}{p-value. It will be the p-value based on robust methods. }
  	\item{p.value_singlevariant}{p-value for each single variant in this region-based test.}
	\item{mac}{total minor allele count (MAC).}
  	\item{param$n.marker}{a number of SNPs in the genotype matrix.}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
  	\item{param$rho}{the \eqn{\rho} parameter for all variants. }
	\item{test.snp.mac}{a vector of minor allele count (MAC) of the snps tested. The name is SNP-ID. } 


}


\author{Zhangchen Zhao}

\references{

Zhao, Z., Bi, W., Zhou, W., VandeHaar, P., Fritsche, L. G., & Lee, S. (2019). UK-Biobank Whole Exome Sequence Binary Phenome Analysis with Robust Region-based Rare Variant Test. \emph{The American Journal of Human Genetics}, in press.

}

\examples{


data(SKATBinary.example)
Z<-SKATBinary.example$Z

obj<-SKAT_Null_Model(y ~ x1 + x2, out_type="D", data=SKATBinary.example)

# run SKAT (default method) with Hybrid
out = SKATBinary_Robust(Z, obj)

# p-value
out$p.value

#
#	Run burden and SKAT

SKATBinary_Robust(Z, obj, method="Burden")$p.value
SKATBinary_Robust(Z, obj, method="SKAT")$p.value



}

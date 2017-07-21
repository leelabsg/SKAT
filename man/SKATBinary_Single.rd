 \name{SKATBinary_Single}
 \alias{SKATBinary_Single}
 \title{Single variant tests for binary traits with Firth and efficient resampling methods}
 \description{
     This function computes p-values of single variant test using the firth and efficient resampling methods.   
 }
 \usage{
 
 
	SKATBinary_Single(Z, obj, method.bin="Hybrid"
	, impute.method = "bestguess", is_check_genotype=TRUE, is_dosage = FALSE
	, missing_cutoff=0.15, max_maf=1, estimate_MAF=1
	, N.Resampling=2*10^6, seednum=100, epsilon=10^-6)


 }
\arguments{
      \item{Z}{a numeric genotype vector. Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, 
      where A is a major allele and a is a minor allele. }
      \item{obj}{output object from SKAT_Null_Model. }
      \item{method.bin}{a type of method to compute a p-value (default="Hybrid"). See details.}
	  \item{impute.method}{a method to impute missing genotypes (default= "bestguess"). }
      \item{is_check_genotype}{a logical value indicating whether to check the validity of the genotype matrix Z (default= TRUE). See SKAT page for details. }
      \item{is_dosage}{a logical value indicating whether the matrix Z is a dosage matrix. If it is TRUE, SKAT will ignore ``is_check_genotype''. }
      \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{max_maf}{a cutoff of the maximum minor allele frequencies (MAF) (default=1, no cutoff). Any SNPs with MAF > cutoff will be excluded from the analysis.}
      \item{estimate_MAF}{a numeric value indicating how to estimate MAFs for the weight calculation and 
      the missing genotype imputation. See SKAT page for details }
      \item{N.Resampling}{a number of resampling to be conducted to get p-values (default=2 *10^6).}
      \item{seednum}{a seed number for random number generation (default=100). If NULL, no seed number will be assigned.}
      \item{epsilon}{a precision level  (default=10^-6).}
	  

}
\value{
	\item{p.value}{p-value. It will be the mid p-value if ER is used to compute the p-value.}
	\item{p.value.standard}{(ER only) standard p-value.}
	\item{p.value.resampling}{p-values from resampled outcome. You can obtain it when n.Resampling in SKAT_Null_Model was > 0. See the SKAT_Null_Model. }
	\item{p.value.standard.resampling}{(ER only) standard p-values from resampled outcome.}
	\item{m}{the number of individuals with minor alleles.}
	\item{MAP}{the minimum possible p-values. It is available when the method.bin="ER" and m is sufficiently small.}
	\item{MAC}{the total minor allele count (MAC).}
	\item{n.total}{(ER only) the number of resampling to be generated to get the p-value. 
	It can be smaller than N.Resampling when the total number of configurations of case-controls among individuals with minor alleles are smaller than
	N.Resampling.}
  	\item{is.accurate}{logical value for the accuracy of the p-value. If it is false, more resampling is needed to accurately estimate the p-value. }
	\item{method.bin}{a type of method to be used to compute the p-value.}	
}
\details{

This function implements three methods (method.bin) to compute p-values: 1) Efficient resampling (ER);
2) Firth biased adjusted likelihood ratio test (Firth);  and 3) Hybrid. 
"Hybrid" selects a method based on the total minor allele count (MAC), the number of individuals with minor 
alleles (m), and the degree of case-control imbalance. 

Adaptive ER (ER.A) is not implemented yet. 

If seednum is not NULL, set.seed(seednum) function is used to specify seeds to get the same p-values 
of ER based methods for different runs. Therefore, please set seednum=NULL, if you do not want to set seeds. 


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
out = SKATBinary_Single(Z[,1], obj)

# p-value
out$p.value

# MAP
out$MAP

# method used to compute p-value (method.bin)
out$method.bin


#
#	Use firth method to compute p-value
SKATBinary_Single(Z[,1], obj, method.bin="Firth")$p.value

}



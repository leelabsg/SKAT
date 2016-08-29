 \name{SKAT_CommonRare}
 \alias{SKAT_CommonRare}
 \alias{SKAT_CommonRare.SSD.OneSet}
 \alias{SKAT_CommonRare.SSD.OneSet_SetIndex}
 \alias{SKAT_CommonRare.SSD.All}
 \title{SKAT for the combined effect of common and rare variants}
 \description{
     Sequence Kernel association test for the combined effect of common and rare variants.      
 }
 \usage{

	SKAT_CommonRare(Z, obj, weights.beta.rare=c(1,25), weights.beta.common=c(0.5,0.5)
	, method="C", r.corr.rare=0, r.corr.common=0, CommonRare_Cutoff=NULL
	, test.type="Joint", is_dosage=FALSE, missing_cutoff=0.15
	, estimate_MAF=1, SetID1=NULL)


	SKAT_CommonRare.SSD.OneSet(SSD.INFO, SetID, obj, \dots)

	SKAT_CommonRare.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, \dots )

	SKAT_CommonRare.SSD.All(SSD.INFO, obj, \dots)
	
 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. 
      Missing genotypes will be imputed by the simple Hardy-Weinberg equilibrium (HWE) based imputation. }
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{weights.beta.rare}{a numeric vector of parameters of beta weights for rare variants (default=c(1,25)).}
      \item{weights.beta.common}{a numeric vector of parameters of beta weights for common variants (default=c(0.5,0.5)).}
            
      \item{method}{a method to combine common and rare variant effects (default="C"). "C" represents the combined sum test,
      and "A" represents the adaptive sum test. "AR" represents a different type of adaptive test in which 
      common variants are projected over rare variants. }
      \item{r.corr.rare}{the \eqn{\rho} parameter for rare variants (default= 0). \eqn{\rho} =0 and 1 indicate SKAT and Burden test, respectively}
      \item{r.corr.common}{the \eqn{\rho} parameter for common variants (default= 0). \eqn{\rho} =0 and 1 indicate SKAT and Burden test, respectively}
      \item{CommonRare_Cutoff}{MAF cutoff for common vs rare variants (default=NULL). It should be a numeric value between 
      0 and 0.5, or NULL. When it is NULL, \eqn{1/ \sqrt{2 SampleSize }} will be used. }
      \item{test.type}{a string to indicate test type (default="Joint"). "Joint" indicates the joint test of the 
      combined effects of common and rare variants. "Rare.Only" and "Common.Only" will conduct test only with rare and common variants, 
      respectively.  }
      \item{is_dosage}{see SKAT}
      \item{missing_cutoff}{see SKAT}
      \item{estimate_MAF}{see SKAT}
      \item{SetID1}{internal use only}
      
      
      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD. }
      \item{SetID}{a character value of Set ID. You can find a set ID of each set from SetInfo object of SSD.INFO}
      \item{SetIndex}{a numeric value of Set index. You can find a set index of each set from SetInfo object of SSD.INFO  }
      
      \item{\dots}{ furthuer arguments to be passed to ``SKAT_CommonRare'' }
      
}
\value{
	\item{p.value}{p-value. }
	\item{p.value.resampling}{p-values from resampled phenotypes. You can get it when you use obj from SKAT_Null_Model function with resampling. See the SKAT_Null_Model. }
	\item{n.rare}{the number of rare variants used for the test}
	\item{n.common}{the number of common variants used for the test}	
	\item{Cutoff}{the MAF cut-off to divide common and rare variants}	
		
  	\item{Q}{the test statistic of SKAT. It has NA when method="A" or "AR".}
	\item{param}{estimated parameters of each method.}   
	\item{param$Is_Converged}{an indicator of the convergence. 1 indicates the method is converged, and 0 indicates the method is not converged. 
	When 0 (not converged), "liu.mod" method is used to compute p-value. }  
	\item{param$n.marker}{a number of SNPs in the genotype matrix}  
	\item{param$n.marker.test}{a number of SNPs used for the test. It can be different from param$n.marker when 
	some markers are monomorphic or have higher missing rates than the missing_cutoff. } 
	
	\item{results}{(SKAT_CommonRare.SSD.All only) the dataframe that contains SetID, p-values (P.value), 
	the number of markers in the SNP sets (N.Marker.All), 
	the number of markers to test for an association after excluding non-polymorphic or high missing rates markers (N.Marker.Test), 
	and the number of rare (N.Marker.Rare) and common (N.Marker.Common) variants used for association tests.  }
	\item{P.value.Resampling}{(SKAT_CommonRare.SSD.All only) the matrix that contains p-values of resampled phenotypes. }
	
}
\details{

The small sample adjustment for binary traits is not implemented for "A" and "AR". 

}


\author{Seunggeun Lee}

\references{

Ionita-Laza, I.*, Lee, S.*, Makarov, V., Buxbaum, J. Lin, X. (2013). 
Sequence kernel association tests for the combined effect of rare and common variants.  
\emph{American Journal of Human Genetics}, 92, 841-853. 
* contributed equally. 

}

\examples{


data(SKAT.example)
attach(SKAT.example)



# continuous trait
obj<-SKAT_Null_Model(y.c ~ X, out_type="C")
SKAT_CommonRare(Z, obj)$p.value
SKAT_CommonRare(Z, obj, method="A")$p.value
SKAT_CommonRare(Z, obj, method="AR")$p.value


# dichotomous trait 
obj<-SKAT_Null_Model(y.b ~ X, out_type="D")

# Combined sum test in the manuscript (SKAT-C and Burden-C)
SKAT_CommonRare(Z, obj)$p.value
SKAT_CommonRare(Z, obj, r.corr.rare=1, r.corr.common=1 )$p.value


# Test only with common variant
SKAT_CommonRare(Z, obj, test.type="Common.Only")$p.value

# Test only with rare variant
SKAT_CommonRare(Z, obj, test.type="Rare.Only")$p.value


# Use CommonRare_Cutoff=0.01 instead of CommonRare_Cutoff = NULL
SKAT_CommonRare(Z, obj, CommonRare_Cutoff=0.01)$p.value



}



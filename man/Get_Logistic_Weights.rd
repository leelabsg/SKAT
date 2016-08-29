 \name{Get_Logistic_Weights}
 \alias{Get_Logistic_Weights}
 \alias{Get_Logistic_Weights_MAF}
 \title{Get the logistic weight}
 \description{

     Get logistic weights from either a genotype matrix (Z) or a vector of minor allele frequncies (MAF). 
     Users can apply this weights to SKAT by giving it as the ``weights'' parameter. 
     The logistic weight gives equal weights to rare variants and nearly zero weight to common variants. 

 }
 \usage{

	Get_Logistic_Weights(Z, par1=0.07, par2=150)


	Get_Logistic_Weights_MAF(MAF, par1=0.07, par2=150)
		
 }
\arguments{
       \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
       Each genotype should be coded as 0, 1, 2, and 9 for AA, Aa, aa, and missing, where A is a major allele and a is a minor allele. }
       \item{MAF}{a numeric vector of minor allele frequncies. }
       \item{par1}{a numeric value of the first parameter of the logistic weight (default= 0.07).}
       \item{par2}{a numeric value of the second parameter of the logistic weight(default= 150).}
}
\value{
	A vector of the logistic weight. 
}
\details{
	The formula for the weight is 
	\deqn{ weights = \frac{e^{(par1 - MAF) par2}}{1 + e^{(par1 - MAF) par2}}. }                                              
}


\author{Seunggeun Lee}


\examples{


data(SKAT.example)
attach(SKAT.example)


#############################################################
#	Compute the P-value of SKAT with the logistic Weight (par1=0.07, par2=150)

# Use logistic weight
obj<-SKAT_Null_Model(y.c ~ X, out_type="C")
weights<-Get_Logistic_Weights(Z, par1=0.07, par2=150)
SKAT(Z, obj, kernel = "linear.weighted", weights=weights)$p.value

# Weights function
MAF<-colMeans(Z)/2
plot(MAF,weights)


}



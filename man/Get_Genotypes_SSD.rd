 \name{Get_Genotypes_SSD}
 \alias{Get_Genotypes_SSD}
 \alias{Get_Genotypes_SSD_Sparse}
 \title{Get Genotype data from SSD file}
 \description{
	Read a SSD file and return a genotype matrix.
 }
 \usage{
	Get_Genotypes_SSD(SSD_INFO, Set_Index, is_ID = TRUE)
	
	Get_Genotypes_SSD_Sparse(SSD_INFO, Set_Index)
 }
\arguments{
      \item{SSD_INFO}{SSD_INFO object returned from Open_SSD.}
      \item{Set_Index}{a numeric value of Set index. The set index of each set can be found from SetInfo object of SSD.INFO. }
      \item{is_ID}{a logical value indicating whether to read SNP ID (default=TRUE). If TRUE, it reads SNP IDs and use them as column names.}
}
\value{
 	A genotype matrix with n rows and m columns, where n is the number of samples and m is the number of SNPs. Get_Genotypes_SSD_Sparse returns a sparse matrix. Get_Genotypes_SSD_Sparse always returns SNP IDs as the column names, so does not have is_ID parameter. 
}
\author{Seunggeun Lee, Larisa Miropolsky}


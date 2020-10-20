\name{Read_SNP_WeightFile}
\alias{Read_SNP_WeightFile}
\title{Read a file with custom weights}
\description{
     Read a file with custom weights
}
\usage{
	Read_SNP_WeightFile(FileName)

}
\arguments{
      \item{FileName}{input file name of a custom weight.}
     
}
\value{
	Output object has a hash table of SNP IDs and weights.
}
\details{
 The file should be a white-space (space or tab) delimitered file with 2 columns:
  SNP_ID and weight value.

 Please keep in mind that there should be no header!!
                                                        
}


\author{Seunggeun Lee}




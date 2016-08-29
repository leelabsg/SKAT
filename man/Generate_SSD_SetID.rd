 \name{Generate_SSD_SetID}
 \alias{Generate_SSD_SetID}
 \title{Generate SNP set data file (SSD) }
 \description{
	Generate a SNP set data file (SSD) from binary plink data files using user specified SNP sets. 

 }
 \usage{
	Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID,
	 File.SSD, File.Info, Is.FlipGenotype=TRUE)
 }
\arguments{
      \item{File.Bed}{name of the binary ped file (BED).}
      \item{File.Bim}{name of the binary map file (BIM).}
      \item{File.Fam}{name of the FAM file (FAM).}
      \item{File.SetID}{name of the SNP set ID file that defines SNP sets. 
      The first column must be Set ID, and the second column must be SNP ID. There should be no header!! }
      \item{File.SSD}{name of the SSD file generated. }
      \item{File.Info}{name of the SSD info file generated. }
      \item{Is.FlipGenotype}{internal use only, please do not change}
      
}

\details{
 The SetID file is a white-space (space or tab) delimitered file with 2 columns:
  SetID and SNP_ID.

 Please keep in mind that there should be no header!
 The SNP_IDs and SetIDs should be less than 50 characters, otherwise, it will return an error message.
      
 The SSD file is a binary formated file with genotypes. 
 The SSD info file is a text file with general information on data and SNP sets (first 6 rows), 
 and information on each set (after 8th row).
         
}


\author{Seunggeun Lee, Larisa Miropolsky}


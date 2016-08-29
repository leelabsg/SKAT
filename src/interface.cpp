 /*************************************************************
 *
 * Final Interface to R
 *
 **************************************************************/

#include <R.h>

void Generate_MWA_SetID_File(char* Bed, char* Bim, char* Fam, char* SetID, char* Mwa, char* Info, int MAFConvert,int* myerror) ;
void Kill_MWA_SetID_File();

void Open_MWA(char* MWA_File, char* Info,int* myerror);
void Close_MWA() ;
void Get_TotalNumberofSets(int* Total_Num_SNPSets);
void Get_TotalNumberofSNPs(int* Total_Num_SNP);
void Get_TotalNumberofInd( int* Total_Num_IND);
void Get_NumberofSnps(int SetID,int *Num_SNP,int* myerror);
void Get_Genotypes( int Set_number, int* Z, int size, int Is_MakeFile,int* myerror);
void Get_Genotypes_withID( int Set_number, int* Z, char * SNPID,int size, int Is_MakeFile, int* myerror);

extern "C" {

//===============================================================
//Generate_MWA_SetID_File(Bim, Bed, SetID_File)
	//Generate MWA file from SetID file. 
	//Bim : Bim file name
	//Bed : Bed file name
	//SetID_File : Set ID file. The first column is setid, and the second column is snp id

void R_Generate_MWA_SetID_File(
char** Bed, char** Bim, char** Fam, char** SetID, char** Mwa, char** Info, int*  MAFConvert, int * err) 
{

	//Rprintf("Bed[%s]\n",Bed[0]);
	//Rprintf("Bim[%s]\n",Bim[0]);
	Generate_MWA_SetID_File(Bed[0], Bim[0], Fam[0], SetID[0], Mwa[0], Info[0], *MAFConvert, err) ;
}

void R_Kill_MWA_SetID_File()
{
	Kill_MWA_SetID_File();
}

//===============================================================
//===============================================================
//===============================================================
//===============================================================
//	2. Open and Close MWA Files
//===============================================================
//Open_MWA(MWA_File, *MWA_File_ID)
	//Open an existing MWA file. It returns MWA_File_ID
	//MWA_File_Path : MWA File Name
	//MWA_File_ID : integer value of file ID. 

//USAGE FROM R:
//Open_MWA(MWA_File, MWA_FILE_ID)

void R_Open_MWA(char** MWA_File, char** Info, int * err)
{
	Open_MWA(MWA_File[0], Info[0],err);

}

void R_Close_MWA() 
{
	Close_MWA();
}
//===============================================================
//===============================================================
//===============================================================
//===============================================================
//	3. Get DATA from MWA files
//===============================================================
	
//Get_TotalNumberofSets(MWA_FILE_ID, *Total_Num_SNPSets)
	//Return Total number of SNP Sets
	//MWA_File_ID : integer value of file ID. 
	//Total_Num_SNPSets : Total Number of SNP Sets in the MWA file. 


void R_Get_TotalNumberofSets(int* Total_Num_SNPSets)
{
	Get_TotalNumberofSets(Total_Num_SNPSets);
}
//===============================================================
	
//Get_TotalNumberofSNPs(MWA_FILE_ID, *Total_Num_SNP)
	//Return Total number of SNP
	//MWA_File_ID : integer value of file ID. 
	//Total_Num_SNP : Total Number of SNPs in the MWA file. 


void R_Get_TotalNumberofSNPs(int* Total_Num_SNP)
{
	Get_TotalNumberofSNPs(Total_Num_SNP);
}

//===============================================================
//Get total number of individuals
//Get_TotalNumberofInd(MWA_FILE_ID, *Total_Num_IND) 
	//Return Total number of IND
	//MWA_File_ID : integer value of file ID. 
	//Total_Num_SNP : Total Number of individuals in the MWA file. 

void R_Get_TotalNumberofInd( int* Total_Num_IND)
{
	 Get_TotalNumberofInd( Total_Num_IND);
}



//===============================================================
	
//Get_NumberofSnps(MWA_FILE_ID, SetID,*Num_SNP)
	//Return a number of SNPs in the given SNP Set. 
	//MWA_File_ID : integer value of file ID. 
	//SetID : SetID
	//Num_SNP : number of SNPs in the given SNP Set. 


void R_Get_NumberofSnps(int* SetID,int *Num_SNP, int * err)
{
	Get_NumberofSnps(*SetID,Num_SNP,err);
} 

//===============================================================

//RETURN ONE LONG ARRAY - IT WILL BE MATRIX IN R


void R_Get_Genotypes( int *Set_number, int * Z , int * size, int *Is_MakeFile, int * err) // set_number base on INFO file. The result will be printed to file.
{
	Get_Genotypes( *Set_number, Z, *size,  * Is_MakeFile, err);
}	


void R_Get_Genotypes_withID( int *Set_number, int * Z , char * SNPID, int * size, int *Is_MakeFile, int * err) // set_number base on INFO file. The result will be printed to file.
{
	Get_Genotypes_withID( *Set_number, Z, SNPID,  *size,  * Is_MakeFile, err);
}	



} // extern "C"



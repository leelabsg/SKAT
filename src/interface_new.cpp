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

void Kernel_IBS_1(int * Z, int * pn, int * pp, double * Kernel);
void Kernel_IBS_Weight_1(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel);
void Kernel_2wayIX_1(int * Z, int * pn, int * pp, double * Kernel);

void SKAT_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);
void SKATO_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r, double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);
void GetProb(int k, int ngroup, int ncase, int * group, double * weight, double * prob);
void SKAT_Permu(double *Z, int *Y, int nSNP, int nSample, int nPermu, double * pval, double * pval_same, double epsilon);

void  qfc_1(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);

void SL_Binary_Boot_2(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *Z, int *err );
void ResampleSTAT_1(double * Z0, double *Z1, double * Z0_C, double * Z1_C, 
                  double * teststat_Z0, double *teststat_Z1, double *pteststat_Z0_C, double *pteststat_Z1_C,
                  double * r_corr, int *pn_r, int *pk, int *pm, int * pn,
                  int * total_k, int * ncase_k, double * p1,
                  int *buf1, int * buf2, int *buf3, double * teststat_one, /* buffers */
                  double * Q, int *err);

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

/* Kernel */

void Kernel_IBS(int * Z, int * pn, int * pp, double * Kernel){
	Kernel_IBS_1(Z, pn, pp, Kernel);
}
void Kernel_IBS_Weight(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel){
	Kernel_IBS_Weight_1(Z, pn, pp, UseGivenWeight ,  weight, Kernel);
}
void Kernel_2wayIX(int * Z, int * pn, int * pp, double * Kernel){
	Kernel_2wayIX_1(Z, pn, pp, Kernel);
}

/* R_Binary */
void RSKATExact(int * resarray, int * nres, int * nres_k, double * Z0, double *Z1, int * k, int * m, int * total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int * test_type, double * epsilon)
{
    SKAT_Exact(resarray, nres[0], nres_k, Z0, Z1, k[0], m[0], total[0], total_k, prob_k, odds, p1, IsExact, pval, pval_same, minP, test_type[0], epsilon[0]);
        
}

void RSKATOExact(int * resarray, int * nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int * n_r, double * param, int * k, int * m, int * total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int * test_type, double * epsilon)
{
    SKATO_Exact(resarray, nres[0], nres_k, Z0, Z1, r_corr, n_r[0], param, k[0], m[0], total[0], total_k, prob_k, odds, p1, IsExact, pval, pval_same, minP, test_type[0], epsilon[0]);
        
}
    
void   RGetProb(int* k, int* ngroup, int* ncase, int * group, double * weight, double * prob){
        
    GetProb(k[0], ngroup[0], ncase[0], group, weight, prob);
}

void   RSKATPermu(double *Z, int *Y, int *nSNP, int *nSample, int *nPermu,double * pval, double *pval_same,  double * epsilon){
        
    SKAT_Permu(Z, Y, nSNP[0], nSample[0], nPermu[0], pval, pval_same, epsilon[0]);
}

void  qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res){
	
	qfc_1(lb1, nc1, n1, r1, sigma, c1, lim1, acc,trace, ifault, res);

}

void SL_Binary_Boot(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *Z, int *err ){
	
	SL_Binary_Boot_2(pn, pm, pncase, pcase, buf1, buf2, Z, err );

}

void ResampleSTAT(double * Z0, double *Z1, double * Z0_C, double * Z1_C, 
                  double * teststat_Z0, double *teststat_Z1, double *pteststat_Z0_C, double *pteststat_Z1_C,
                  double * r_corr, int *pn_r, int *pk, int *pm, int * pn,
                  int * total_k, int * ncase_k, double * p1,
                  int *buf1, int * buf2, int *buf3, double * teststat_one, /* buffers */
                  double * Q, int *err){
	
	ResampleSTAT_1(Z0, Z1, Z0_C, Z1_C, teststat_Z0, teststat_Z1, pteststat_Z0_C, pteststat_Z1_C, r_corr, pn_r, pk, pm, pn,
                  total_k, ncase_k, p1, buf1, buf2, buf3, teststat_one, Q, err);                
                  
}

} // extern "C"


#include <Rinternals.h>
#include <R_ext/Rdynload.h>

       
static R_NativePrimitiveArgType type1[] = { STRSXP };
static R_NativePrimitiveArgType type2[] = { STRSXP, STRSXP, STRSXP, STRSXP,STRSXP,STRSXP,INTSXP,INTSXP };
static R_NativePrimitiveArgType type3[] = {STRSXP, STRSXP, INTSXP};
static R_NativePrimitiveArgType type4[] = { INTSXP };
static R_NativePrimitiveArgType type5[] = {INTSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType type6[] = {INTSXP, INTSXP, INTSXP, INTSXP,INTSXP};
static R_NativePrimitiveArgType type7[] = {INTSXP, INTSXP, RAWSXP, INTSXP, INTSXP,INTSXP};

static R_NativePrimitiveArgType type8[] = {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,  REALSXP, REALSXP, INTSXP, REALSXP,REALSXP,REALSXP,INTSXP, REALSXP};
static R_NativePrimitiveArgType type9[] = {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP,  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType type10[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType type11[] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,REALSXP, REALSXP,  REALSXP};
static R_NativePrimitiveArgType type12[] = {INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType type13[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType type14[] = {INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType type15[] = {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType type16[] = {INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType type17[] = {REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, 
			INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, 
                  REALSXP, INTSXP};


R_CMethodDef cMethods[] = {
 	{"R_Generate_MWA_SetID_File", (DL_FUNC) &R_Generate_MWA_SetID_File, 8, type2},
  	{"R_Kill_MWA_SetID_File", (DL_FUNC) &R_Kill_MWA_SetID_File, 0 },
   	{"R_Open_MWA", (DL_FUNC) &R_Open_MWA, 3, type3 },
   {"R_Close_MWA", (DL_FUNC) &R_Close_MWA, 0 }, 
   {"R_Get_TotalNumberofSets", (DL_FUNC) &R_Get_TotalNumberofSets, 1, type4 }, 
   {"R_Get_TotalNumberofSNPs", (DL_FUNC) &R_Get_TotalNumberofSNPs, 1, type4 }, 
   {"R_Get_TotalNumberofInd", (DL_FUNC) &R_Get_TotalNumberofInd, 1, type4 }, 
   {"R_Get_NumberofSnps", (DL_FUNC) &R_Get_NumberofSnps, 3, type5 }, 
//
   {"R_Get_Genotypes", (DL_FUNC) &R_Get_Genotypes, 5, type6},
   {"R_Get_Genotypes_withID", (DL_FUNC) &R_Get_Genotypes_withID, 6, type7},
//   
   {"RSKATExact", (DL_FUNC) &RSKATExact, 18, type8},
   {"RSKATOExact", (DL_FUNC) &RSKATOExact, 21, type9},
   {"RGetProb", (DL_FUNC) &RGetProb, 6, type10},
//
   {"RSKATPermu", (DL_FUNC) &RSKATPermu, 8, type11}, 
   {"Kernel_IBS", (DL_FUNC) &Kernel_IBS, 4, type12},
   {"Kernel_IBS_Weight", (DL_FUNC) &Kernel_IBS_Weight, 6, type13},   
//
   {"Kernel_2wayIX", (DL_FUNC) &Kernel_2wayIX, 4, type14},  
   {"qfc", (DL_FUNC) &qfc, 11, type15},  
//
	{"SL_Binary_Boot", (DL_FUNC) &SL_Binary_Boot, 8, type16},
	{"ResampleSTAT", (DL_FUNC) &ResampleSTAT, 22, type17},

  
   {NULL, NULL, 0}
};


extern "C" {

void R_init_SKAT(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

}


/*************************************************************
 *
 * Moving window project 
 * File: setid_bim_index.h	
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   create indexing - hash table between SetID.txt file and *.Bim file
 *	 [setid_row#] = bim_row#
 **************************************************************/

#ifndef _SETID_BIM_INDEX_H        
#define _SETID_BIM_INDEX_H 
#include <fstream>  
#include <iostream> 
#include <algorithm>
#include <vector>
#include <climits>
#include <cstring>
#include <memory>

#include "DArray.h"
#include "NPsort.h"
#include "error_messages.h"
//===============================================================
//===============================================================

using namespace std;

class SNP_info 
{
    
public:
    string snp_id;
    string A1;
    string A2;
    string Chr;
    long  Pos;
	int total_counter_per_letter[2]; //= {0,0};
	int line_counter_per_letter[2]; //= {0,0};
	int flag;
    
    SNP_info(){
        total_counter_per_letter[0] = 0;
        total_counter_per_letter[1] = 0;
        line_counter_per_letter[0] = 0;
        line_counter_per_letter[1] = 0;
        flag=0;
    }
};


//===============================================================
//===============================================================

class Hasht   // info for column i  the columns i = 0 ... get_noc()-1
{
private:
	 //int *index;
	
    std::string m_setidfile;
    std::string m_bimfile;

	char** m_bimf_snpsid;
	size_t* m_bimf_sorted;


	ofstream m_log;
	ifstream m_setid;
	ifstream m_bim;
	SNP_info* m_snp_sets;
	
	void upload_snpid_from_bim(int * myerror);
	void upload_snpid_from_setid_build_hash(int * myerror);

	int binsearch(const char* source); 

	void Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters = " ");

	int	Get_Num_of_SNPs_in_SetID(int * myerror);


public:	
	Hasht(char * setID, char * bim, char* mwa, int * myerror);
	
	~Hasht();
	SNP_info* get_snps_sets() {return this->m_snp_sets; }

    
	char** m_setidf_setid;
	size_t* m_hash_table;
	size_t m_num_of_snps;//number of snps
	size_t m_num_of_snps_insetid;

	size_t m_num_of_snps_insetid_org ; // added by SLEE 

};

#endif //_SETID_BIM_INDEX_H

/*************************************************************
 *
 * Moving window project 
 * File: bed_reader.h	
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Deal with *.bed file
 *
 * Revised by SLEE
 * Date: May 18, 2017
 **************************************************************/
#ifndef _BED_READER_H        
#define _BED_READER_H 

#include <fstream>  
#include <iostream> 
#include "setid_bim_index.h"

using namespace std;

//===============================================================
//===============================================================

#define MY_CHAR_BIT 8


class BedFileReader
{
public:
	BedFileReader(char* f, char* m, char* fm, char* o, Hasht* ht, int* myerror, char* info = NULL, int MAFConvert=1);

  // Added by SLEE
  BedFileReader(char* f, char* m, char* fm, int* myerror);
  void read_One_SNP(int SNPindex, int * genotype, int* myerror);
  void close_bed();
    
  
  
	~BedFileReader();
	size_t m_set_counter; 
	size_t m_num_of_snps_insetid;
	


private:

    // File names
	std::string m_filename;
	std::string m_filename_mwo;
	std::string m_filename_bim;
	std::string m_filename_fam;

    std::string m_file_temp_name;
    std::string m_info_file;
 	std::string m_info_rewritten;
    
    // streams for files
	std::ifstream m_file;
	std::ifstream m_bim;
	std::ifstream m_fam;
	std::ofstream m_file_temp;
	std::ofstream m_info;
	std::ifstream m_infoi;

	std::fstream m_file_mwo;


	std::ofstream m_info_rewr;
	size_t m_begin4rw; 
	
	// added by SLEE 

	//int m_encode_output;
  int m_MAFConvert; // =1 change coding for minor allele =0 no change

	size_t m_approx_line_lenght;//number of snps
	SNP_info* m_snp_sets;
	size_t m_line_counter; // 0 .... NumOfIndividuals
	//size_t m_setnumber;
	size_t m_size_of_esi;
	//size_t m_win_size;
	//size_t m_ovrlp_size;

  // Please add here

    void init(const char* bim_file, const char* fam_file, int* myerror);
	void read_data_and_create_mwo_used_hashtable(Hasht* ht, int* myerror);
	void upload_snpid_from_bim(int* myerror);
	void decode_byte(int* bits_val,size_t* individuals_counter, int* temp_snp_info0, int* temp_snp_info1,size_t snp_set_ind);
	void encode(int* temp_snp_info,char* encoded_snp_info);
};

#endif //_BED_READER_H


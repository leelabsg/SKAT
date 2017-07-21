/*************************************************************
 *
 * Moving window project
 * File: setid_bim_index.h
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   create indexing - hush table between SetID.txt file and *.bim file
 *	 [setid_row#] = bim_row#
 *
 * Ver 1.1
 **************************************************************/
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <bitset>
#include <math.h>
#include <fstream>
#include <iostream>

#include <algorithm>

#include "setid_bim_index.h"

//=======================================================================
// This function split "str" by "delimiters" and put the result to "tokens"
// Inputs:
// str - source string
// tokens - target vector of strings
// delimeters - string separator
// Notes:
// "tokens" - vector of strings, so to move inside use:
// tokens.at(0).c_str(); tokens.at(1).c_str().
// important clear "tokens" after each iteration: tokens.clear();
// otherwise it will append new results to existing.
//=======================================================================
void Hasht::Tokenize(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters )
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

//=======================================================================
// Constructor function
// Inputs:
// setID - name of file includes path(full or relative).
//		This file should contain two columns:
//		first setid, second snpid, delimited by "tab" or by "space"
// bim - name of file includes path(full or relative). "*.bim" file
// mwa - name of future file includes path(full or relative). "*.mwa" file.
//		here this "*.mwa" file not will be created and not in use
//		- it just to locate "*.LOG" file in the same place as future "*.mwa".
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=======================================================================
Hasht::Hasht(char * setID, char * bim, char* mwa, int * myerror)
{

	*myerror = NO_ERRORS;

    // Init pointers
    m_setidf_setid=NULL;
    m_hash_table=NULL;
    
    // Init files
	this->m_setidfile = setID;
	this->m_bimfile = bim;
	this->m_num_of_snps_insetid_org = 0;


    
    
    std::string log_filename;
    log_filename += mwa;
    log_filename += "_LOG.txt";
	this->m_log.open(log_filename.c_str());


	this->upload_snpid_from_bim(myerror);
	if (*myerror != 0)
		return;

	this->upload_snpid_from_setid_build_hash(myerror);
	if (*myerror != 0)
		return;

    // Job was done, remove id
	for (size_t i = 0; i < m_num_of_snps; ++ i)
		delete [] this->m_bimf_snpsid[i];
	delete [] this->m_bimf_snpsid;
	delete [] this->m_bimf_sorted;

	this->m_log.close();

}
//=======================================================================
// Destructor function
// Free memory that was dynamically allocated during the run
//=======================================================================
Hasht::~Hasht()
{

	delete [] this->m_hash_table;
	///delete [] this->m_snp_sets;  // Don't remove this line! This delete will be applyed from bed_reader module destructor

	for (int i = 0; i < m_num_of_snps_insetid; ++ i)
		delete [] this->m_setidf_setid[i];
	delete [] this->m_setidf_setid;
}

//=========================================================================
// This function creates lookup table between "setID" file and "*.bim" file
// Reads line by line "setID" file
// For every "snpID" finds ralated line in "*.bim" file
// Inputs:
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=========================================================================
void Hasht::upload_snpid_from_setid_build_hash(int * myerror)
{
	std::string line;
	std::vector<std::string> tokens;
	std::vector<std::string> tokens2, tokens3;

	int re = Get_Num_of_SNPs_in_SetID(myerror);
	if(re== 0){
		*myerror = CANT_OPEN_SETID_FILE4READ;
		return;
	}

	this->m_hash_table = new size_t[this->m_num_of_snps_insetid_org +1];
	this->m_setidf_setid = new char*[this->m_num_of_snps_insetid_org +1];

	this->m_setid.open(this->m_setidfile.c_str());
	if (!this->m_setid)
	{
		*myerror = CANT_OPEN_SETID_FILE4READ;
		return;
	}

	int index = 0;
	this->m_log << "#==== The following SNP IDs requested by SetId file not found in *.Bim file ====#" << std::endl;

	while (!this->m_setid.eof( ) )
    {
		tokens.clear();
		getline(this->m_setid, line);
		Tokenize(line, tokens, " \t\n\r");
		if (tokens.size() < 2)
		{
			tokens.clear();
			Tokenize(line, tokens, " ");
		}


		if (tokens.size() >= 2)
		{
			tokens2.clear();
			Tokenize(tokens.at(1).c_str(), tokens2, " ");
			tokens3.clear();
			Tokenize(tokens2.at(0).c_str(), tokens3, "\r");
			//if one or more spaces at the end of the line -
			//spaces at the end of the line don't be taken in account.
			//tokens2.at(0).c_str() instead tokens.at(1).c_str()
			int result_of_search = this->binsearch(tokens3.at(0).c_str());
			if (result_of_search == -1)
				this->m_log << line << std::endl;
			else
			{
                size_t nlength = tokens.at(0).length();
				this->m_hash_table[index] = result_of_search;
				m_setidf_setid[index] = new char[nlength +1];
                m_setidf_setid[index][nlength] = '\0';
				strncpy (m_setidf_setid[index] , tokens.at(0).c_str(), nlength); //copy the setId
				index ++;
			}
		}
	}
	m_num_of_snps_insetid = index;
	this->m_setid.close();

}

//=========================================================================
// Added by SLEE
//	This function get number of SNPs in the SNP set.
//=========================================================================
int Hasht::Get_Num_of_SNPs_in_SetID(int * myerror){

	std::ifstream temp_file;
	std::string line;
	temp_file.open(this->m_setidfile.c_str());
	if (!temp_file)
	{
		*myerror =  CANT_OPEN_SETID_FILE4READ;
		return 0;
	}

	this->m_num_of_snps_insetid_org = 0;//0;
	while (!temp_file.eof( ) ) {
		getline(temp_file, line);
		this->m_num_of_snps_insetid_org++;
	}

	temp_file.close();
    return 1;
}
//=======================================================================
// This function looking for "source" - snpID from "setID" file
// inside of array of "snpID"s from "*.bim" file
// Inputs:
// source - char* , snpID from "setID" file
// Outputs:
// Index of this "source" inside of array of "snpID"s from "*.bim" file
//=======================================================================
int Hasht::binsearch(const char* source)//int ar[],int size,
{
	int lb = 0,	ub = this->m_num_of_snps-1,	mid;             //lb=>lower bound,ub=>upper bound
	for( ;lb <= ub; )
	{
		mid = (lb + ub) / 2;
		if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) == 0)
			return m_bimf_sorted[mid];  // FOUND!!
		else if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) > 0)
			ub = mid - 1;
		else if(strcmp( this->m_bimf_snpsid[m_bimf_sorted[mid]], source) < 0)
			lb = mid + 1;
	}
	return -1 ; //cout<<"\n SEARCH UNSUCCESSFUL";
}
//=======================================================================
// This function reads "*.bim" file,
// creates "m_bimf_snpsid" array of char* - all snpIDs from "*.bim" file
// creates "m_bimf_sorted" int array of indexes to "m_bimf_snpsid" to sort it.
// creates "m_snp_sets" array of SNP_info objects that will hold info during future "*.mwa" preparation
// of how many genotypes of every type per snp - to calculate minor and major
// Inputs:
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=======================================================================
void Hasht::upload_snpid_from_bim(int * myerror)
{
	std::string line;
	std::vector<std::string> tokens;


	this->m_bim.open(this->m_bimfile.c_str());
	if (!this->m_bim)
	{
        //printf("Error: %s cannot be opend!\n", this->m_bimfile.c_str());
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}

	this->m_num_of_snps = 0;//0;
	while (!this->m_bim.eof( ) )
    {
		getline(this->m_bim, line);
        if(line.length() > 5){
            this->m_num_of_snps++; /* changed in 02/28/2015 */
        }
	}

	this->m_bim.close();
	//----------------------------------------------

 
	m_snp_sets = new SNP_info[this->m_num_of_snps];
/*	for (size_t j = 0; j < m_num_of_snps;++j)
	{

		this->m_snp_sets[j].total_counter_per_letter[0] = 0;
		this->m_snp_sets[j].total_counter_per_letter[1] = 0;
		this->m_snp_sets[j].line_counter_per_letter[0] = 0;
		this->m_snp_sets[j].line_counter_per_letter[1] = 0;

	}
*/
	//----------------------------------------------

	this->m_bim.open(this->m_bimfile.c_str());
	this->m_bim.seekg (0, std::ios::beg);
	if (!this->m_bim)
	{
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}


	m_bimf_snpsid = new char*[this->m_num_of_snps];
	m_bimf_sorted = new size_t[this->m_num_of_snps];

	for (size_t i = 0; i < m_num_of_snps;++i)
    {
		tokens.clear();
		getline(this->m_bim, line);
        //Tokenize(line, tokens, "	");
        Tokenize(line, tokens, " \t\n");
		if (tokens.size() < 6)
		{
			tokens.clear();
			Tokenize(line, tokens, " ");
		}
        if (tokens.size() < 6)
		{
			tokens.clear();
			Tokenize(line, tokens, "\t");
		}
        
        if(tokens.size() < 6){
            *myerror = BIM_FILE_TOKEN_WRONG;
            return;
            
        }
        
        this->m_snp_sets[i].Chr = tokens.at(0);
		this->m_snp_sets[i].snp_id =  tokens.at(1);
        this->m_snp_sets[i].Pos = atol(tokens.at(0).c_str());
        
        this->m_snp_sets[i].A1 = tokens.at(4);
        this->m_snp_sets[i].A1 = tokens.at(5);
        
        // Currently allele name should be one character.
        // For INDEL, (name > 1), it uses I
        /*
        if(tokens.at(4).length() == 1){
            this->m_snp_sets[i].letters[0] = tokens.at(4).c_str()[0];
        } else {
            this->m_snp_sets[i].letters[0] = 'I';
        }
        
        if(tokens.at(5).length() == 1){
            this->m_snp_sets[i].letters[1] = tokens.at(5).c_str()[0];
        } else {
            this->m_snp_sets[i].letters[1] = 'I';
        }
        */
        
        
		this->m_snp_sets[i].flag = 0;
        
        size_t nlength = tokens.at(1).length();
		m_bimf_snpsid[i] = new char[nlength+1];
        m_bimf_snpsid[i][nlength] = '\0';
		strncpy (m_bimf_snpsid[i] , tokens.at(1).c_str(), nlength);
        m_bimf_sorted[i] = i;

	}
	this->m_bim.close();
    
    sort_data::sort<char*, sort_data::char_ptr_less>(m_bimf_snpsid, m_bimf_sorted, m_num_of_snps);
    
	/*sort_data::sort((const void *)m_bimf_snpsid, m_bimf_sorted, m_num_of_snps, (DATA2SORT)D_CHARSTAR, (int)0, (int)0);*/

}
//=======================================================================


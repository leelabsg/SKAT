/*************************************************************
 *
 * Moving window project
 * File: mwo_reader.h
 * Date: Dec 9, 2010
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Deal with *.mwa file
 *
 * ver 1.1
 **************************************************************/


#include <bitset>
#include <fstream>
#include <iostream>
#include <cstring>
#include "mwo_reader.h"
//==========================================================
//Constructor
//Reads INFO file and takes all relevant info from there.
//Also based on INFO creates Offset table to move easily in ".mwa" file.
//Inputs:
// filename - path to ".mwa" file"
// myerror - allocated memory to put there the final information
//		if some error happen during the run.
// info - path to "INFO" file
//=======================================================================
MwoFileReader::MwoFileReader(char* filename, int* myerror, char* info)
{
	*myerror = NO_ERRORS;

	this->m_filename = filename;
	this->m_file.open(this->m_filename.c_str(), std::ios::binary);

	if (!this->m_file)
	{
		*myerror = CANT_OPEN_MWA_FILE4READ;
		return;
	}
    //===============================
	if (info == NULL)
	{
        std::string set_filename;
        set_filename += filename;
        set_filename += "bed.INFO.txt";
		this->m_infoin.open(set_filename.c_str());
	}
	else
		this->m_infoin.open(info);


	if (!this->m_infoin)
	{
		*myerror = CANT_OPEN_INFO_FILE4READ;
		return;
	}

	//===============================
	//upload this info from INFO file

	std::string line;
	std::vector<std::string> tokens;

	//WindowSize
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_win_size = atoi(tokens.at(0).c_str());  //  this->m_win_size = 20; //	WindowSize

	//OverlapSize
    //SLEE add: now it is used for MAF convert status (-999 or 1: convert, 0 no convert)
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_ovlp_size = atoi(tokens.at(0).c_str());  //  this->m_ovlp_size = 2; //	OverlapSize

	//NumberOfSNPs
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_num_of_different_snps = atoi(tokens.at(0).c_str());  //  this->m_num_of_different_snps = 1403896; //	NumberOfSNPs

	//NumberOfIndividuals
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_num_of_individuals = atoi(tokens.at(0).c_str());  //  this->m_num_of_individuals = 162; //	NumberOfIndividuals


	//DECODED NumberOfDECODEDbytesPerSNP(NotIncludesSnpID&SpaceAfter_Includes\n)
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_num_of_bytes_per_line = atoi(tokens.at(0).c_str());  // 	this->m_num_of_bytes_per_line = 42;//per line	DECODED NumberOfDECODEDbytesPerSNP(NotIncludedSnpID&SpaceAfter_Included\n)


	//TotalNumberOfSets
	tokens.clear(); getline(this->m_infoin, line);
	Tokenize(line, tokens, "	");
	this->m_total_num_of_sets = atoi(tokens.at(0).c_str());  //  this->m_total_num_of_sets = 1000; //

	if ( this->m_win_size != -999 && this->m_ovlp_size != -999 )
	{
		//this number includes the overlaping, it means that m_total_num_of_snps >= m_num_of_different_snps
		this->m_total_num_of_snps = this->m_win_size +
							((this->m_num_of_different_snps - this->m_win_size)/(this->m_win_size - this->m_ovlp_size)) * this->m_win_size +
							((this->m_num_of_different_snps - this->m_win_size)%(this->m_win_size - this->m_ovlp_size)) + this->m_ovlp_size;
	}
	else
	{
		this->m_total_num_of_snps = this->m_num_of_different_snps; //as in setid file
	}

	//=======================
	this->m_offsetarr = new size_t [this->m_total_num_of_sets];
	this->m_set_size = new size_t [this->m_total_num_of_sets];
	upload_offsets_table();

	this->m_infoin.close();
}
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
void MwoFileReader::Tokenize(const std::string& str,
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
//==========================================================
//Takes from INFO file:
//  1) offset information and put in into m_offsetarr array
//  2) set size for every set and put into m_set_size array
//==========================================================
void MwoFileReader::upload_offsets_table()
{
	std::string line;
	std::vector<std::string> tokens;

	//fill m_offsetarr
	getline(this->m_infoin, line); //skip the separator "===========" line
	for (size_t i = 0; i < this->m_total_num_of_sets; ++i)
	{
/*
[1] -> [0]
[3] -> [1]
[6] -> [2]
[8] -> [3]
*/		tokens.clear();
		getline(this->m_infoin, line);
		Tokenize(line, tokens, "	");
		//this->m_offsetarr[i] = atoi(tokens.at(3).c_str());
		this->m_offsetarr[i] = atoi(tokens.at(1).c_str());

		if( this->m_win_size != -999)
			this->m_set_size[i] = this->m_win_size;
		else
			this->m_set_size[i] = atol(tokens.at(3).c_str());
			//this->m_set_size[i] = atoi(tokens.at(8).c_str());
	}
}
//==========================================================
//Brings the genotype information from "*.mwa"  for specific set "set_num"
//Puts the information to Z - output array
//Prints the information to file if Is_MakeFile = 1
//Inputs:
//set_num - integer, number of set ( base on INFO file) that you are asking for information for it.
//Z - allocated in advance int array for output
//size - size of Z
//myerror - allocated memory to put there the final information
//		if some error happen during the run.
//Is_MakeFile - flag that shows if print info to file or not. Is_MakeFile = 1 print, Is_MakeFile = 0 don't
//		printed file will be created in same directory as ".mwa" file
//==========================================================
void MwoFileReader::get_set(size_t set_num,  int* Z,  size_t size, int* myerror, int Is_MakeFile, char * SNPID)
{
	*myerror = NO_ERRORS;
	if (set_num > 0 && set_num < this->m_total_num_of_sets + 1)
	{
		if (Z == NULL)
		{
			*myerror = WRONG_ALLOCATED_SIZE4OUTPUT_ARRAY;
			return; //ERROR: THE ARRAY NOT ALLOCATED YET
		}
	}
	else
	{
		*myerror = REQUESTED_SET_ID_DOESNOT_EXISTS;
		return;  //ERROR: - WRONG SET NUMBER
	}


	//base on this->m_offsetarr
    string temp_snp_id;
    temp_snp_id.reserve(1000);
	//char temp_snp_id[SNP_ID_SIZE];
	//memset(temp_snp_id,'\0',sizeof(temp_snp_id));
	bool end_of_file = false;
	char* buff = new char[1000];

	bool flag_snpid_done = false;
	bool flag_read_line_done = true;
	size_t snp_ind = 0;
	snpset* ss = new snpset;
	snp* msnp = new snp;
	//memset(msnp->m_name, '\0', sizeof(msnp->m_name));
	
	
	char* ch = new char;
	size_t char_counter = 0;
	size_t snp_id_ch_ind = 0;
	//==================================

	// Changed by Seunggeun
	if(this->m_file.eof()){
		this->m_file.clear();
	}

	this->m_file.seekg(this->m_offsetarr[set_num-1],std::ios::beg);
	//==================================

	while (!end_of_file)
	{
		memset(buff, '\0', 1000);
		this->m_file.read(buff,1000);

		for(int i=0; i<1000; i++)   //process the buff info
		{
			/* LSG changed */
			/*
			if(snp_ind == m_total_num_of_snps)
			{
 				end_of_file = 1;
				break;
			}
			else
			*/
			
			if(!flag_snpid_done)
			{
				if(buff[i] == '\n')  // if first char in line == '\n' - next_snp_set  will start at next char
				{
					prepare_out_array_print_snpset_to_file(ss,set_num, Z,size, Is_MakeFile, myerror, SNPID);
					if(*myerror != 0)
						return;
					//=============================
					end_of_file = true; // stop it after current set - read just one set in the middle of the file
					break;
					//=============================
					//ss = new snpset;
					//snp_id_ch_ind = 0;
					//continue;
				}//if first char in this line == '\n' - next_snp_set

 				/* changed by LSG */
				if(snp_ind == m_total_num_of_snps)
				{
	 				end_of_file = true;
					break;
				}

				if(buff[i] != ' ')
				{
					//temp_snp_id[snp_id_ch_ind] = buff[i];//TODO read + add to temp_snp_id until ' '
					//snp_id_ch_ind++;
                    temp_snp_id += buff[i];
				}
				else  // finished read snip_id - the string at the start of every line
				{
					//temp_snp_id[snp_id_ch_ind] = '\0';
					//strncpy (msnp->m_name,temp_snp_id, SNP_ID_SIZE-1);
					//memset(temp_snp_id,'\0',sizeof(temp_snp_id));
                    
                    msnp->m_name= temp_snp_id;
                    temp_snp_id.clear();
					snp_id_ch_ind = 0;
					flag_snpid_done = 1;
					flag_read_line_done = 0;

				}

			}
			else if(!flag_read_line_done)  //read line of specific set
			{
				if(char_counter == m_num_of_bytes_per_line-1 && buff[i] == '\n')  //last char in curr line '\n'
				{
					flag_read_line_done = true;
					flag_snpid_done = false;
					ss->m_snp.Add(msnp);
					msnp = new snp;
					snp_ind++;
					char_counter = 0;
				}
				else        // read items of specific line of specific set
				{
					*ch = buff[i];
					msnp->m_char.Add(ch);
                    ch = new char;
					char_counter++;
				}
			}

		}
	}

    
    // remove all 
    for (size_t i = 0; i < ss->m_snp.GetSize(); ++i){
        ss->m_snp.GetAt(i)->m_char.Free();        
    }
    ss->m_snp.Free();
    delete ss;
	return;



}

//=========================================================================
//This function prepares output array Z, and printing to file genotype info
//Inputs:
//ss - pointer to snpset
//set_num -
//Z - output array
//Zsise - size of Z
//Is_MakeFile - flag print to file (=1) or not (=0)
//myerror - allocated memory to put there the final information
//		if some error happen during the run.
//=========================================================================
void MwoFileReader::prepare_out_array_print_snpset_to_file(snpset* ss, int set_num, int* Z, size_t Zsize,
														   int Is_MakeFile, int* myerror, char * SNPID)
{
	if (Zsize != (this->m_num_of_individuals * this->m_set_size[set_num - 1]))
	{
		*myerror = WRONG_ALLOCATED_SIZE4OUTPUT_ARRAY;
		return;

	}
	size_t Zind = 0;
    std::string set_filename;
	std::ofstream myout;

	if (Is_MakeFile)
	{
        set_filename = this->m_filename;
        set_filename += ".SET";
		char buffer[1000];
		sprintf(buffer,"%d",set_num);
        set_filename += buffer;
		myout.open(set_filename.c_str() , std::ios::binary);
		if (!myout)
		{
			*myerror = WARNING_CANT_OPEN_FILE4WRITE_2PRINTSET;
			Is_MakeFile = 0;
		}

	}

	int bits_val[MY_CHAR_BIT];
	size_t ind_count = 0;
	size_t ind_count_prev=0;
	char buff[9];

    /* modified by LSG */

	for (size_t i = 0; i < ss->m_snp.GetSize(); ++i)
	{


		if (Is_MakeFile)
				myout << ss->m_snp.GetAt(i)->m_name << " ";

		if(SNPID != NULL){
			int start_id = SNP_ID_SIZE_MAX * i;
			strncpy(SNPID + start_id, ss->m_snp.GetAt(i)->m_name.c_str(), SNP_ID_SIZE_MAX-1);
			//printf("NAME: %s\n", ss->m_snp.GetAt(i)->m_name);
					
		}
		
		for (size_t j = 0; j < ss->m_snp.GetAt(i)->m_char.GetSize(); ++j)
		{
			//DECODE HERE - PRINT DECODED
			//===============================================================
			//=== This part converts Byte "buff[i]" to bits values "bits_val"
			//=== for example byte buff[0] = "w" ->  bits_val = 11101110
			memset((void *)bits_val, 0, sizeof(bits_val));
			int k = MY_CHAR_BIT;  //8
			while (k > 0)
			{
				-- k;
				bits_val[k] = (*(ss->m_snp.GetAt(i)->m_char.GetAt(j))&(1 << k) ? 1 : 0);
			}
			//here interpret Bit information "bits_val" to snps and count it - decode it
			ind_count_prev = ind_count;
			decode_byte(bits_val, buff, &ind_count );
			if (Is_MakeFile)
				myout << buff;

			for(size_t m = 0; m < ind_count - ind_count_prev; ++ m)
			{
				try
				{
					Z[Zind] = atoi(&buff[m * 2]);
					Zind ++;
				}
				catch(...)
				{
					*myerror = WRONG_ALLOCATED_SIZE4OUTPUT_ARRAY;
					return;

				}

			}


		}

        if (Is_MakeFile)
			myout << std::endl;
		ind_count = 0;

	}
	if (Is_MakeFile)
		myout.close();
}

//==========================================================
// This function interpret Bit information "bits_val" to snps
// and count it - decode it
// Inputs:
// bits_val	- bits representation of some character
//		for example bits_val = 01000100
//		bits_val shoul be read in reverse order in couples:
//		two last, two before last ... :   00(the last), 01(before the last), 00(next to first), 01(the first)
// buff - character buffer - decoded view of bits_val, this it an output of this function
//      for example if bits_val = 01000100 , then buff will be "0 9 0 9 "
//      buff allocated in advance for 4 individuals includes spaces and "\n" or "\0" at the end buff[9]
// ind_count - do until ind_count < this->m_num_of_individuals
//==========================================================
void MwoFileReader::decode_byte(int* bits_val,char* buff, size_t* ind_count)
{
	int g = 0;
//do until ind_count < this->m_num_of_individuals
//this->m_num_of_individuals

	for (int i = 3; i > -1; --i)
	{
		if (*ind_count == this->m_num_of_individuals)
		{
			buff[g] = '\0';
			return;
		}

		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
			buff[g] = '0' ;
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
			buff[g] =  '2' ;
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
			buff[g] =  '9';
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
			buff[g] =  '1';
		*ind_count  += 1;
		g ++;
		buff[g] = ' ';
		g ++;
	}
	buff[g] = '\0';
	return;

}

//==========================================================
//This function return number of snps in in specified set
//Inputs:
//SetID - the specified set
//Num_SNP - reserve place to put there the output -
//		number of snps in in specified set
//myerror - allocated memory to put there the final information
//		if some error happen during the run.
//==========================================================
size_t MwoFileReader::get_NumberofSnps(int SetID,int* myerror)
{
    size_t Num_SNP;
    *myerror = NO_ERRORS;
    if (SetID > 0 && SetID < this->m_total_num_of_sets + 1){
			Num_SNP = this->m_set_size[SetID - 1];
    } else {
        Num_SNP = -9999;
        *myerror = REQUESTED_SET_ID_DOESNOT_EXISTS;
    }
    return Num_SNP;


}
//==========================================================
//Destructor - free all dynamically allocated memory
//==========================================================
MwoFileReader::~MwoFileReader(){
		
		delete [] this->m_offsetarr;
		delete [] this->m_set_size;
		this->m_file.close();
}



/*************************************************************
*
* Moving window project 
* File: bed_reader.cpp	
* Date: Dec 9, 2010
* Author: Larissa Miropolsky
*
* Description:
*   Deal with *.bed file, create *.mwa file and INFO file
*
*************************************************************
* INFO file format:
* The previous format:
lines 1 - 6  - general information
line  7 - separator #=====#
Line  8 - last line (the data separated by tab): 
example:
SET#	1	OFFSET_FROM_BEG	0	BYTES	SET_ID	01	SET_SIZE	94


* New format:
lines 1 - 6  - general information
line  7 - separator #=====#
Line  8 - last line (the data separated by tab):
First column index[0] reffers to SET#, 
second column index [1] shows OFFSET_FROM_BEG in BYTES,
third column index [2] shows SET_ID as described in SetId file,
fourth column index [3] reffers to SET_SIZE: 
example:
1	0	01	94

so, in this code the indexes to read this file changed as following:
[1] -> [0]
[3] -> [1]
[6] -> [2]
[8] -> [3]

**/


#include <cstring>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <fstream>  
#include <iostream> 
#include "bed_reader.h"

//===========================================================================
//Constructor - to create BedFileReader object base on SetID - > lookup table
//Inputs:
//filename - path to "*.bed" file
//mapname - path to "*.bim" file
//famname - path to "*.fam" file
//outfile - path to "*.mwa" file
//ht - lookup table object that created in setid_bim_index.cpp  
//myerror - allocated memory to put there the final information 
//		if some error happen during the run.
//info - path to future ".INFO" file, can work without it - will be created 
//       in same directory as ".bed" file
//===========================================================================
BedFileReader::BedFileReader(char* filename, char* mapname, 
							 char* famname, char* outfile, Hasht* ht, int* myerror, char* info, int MAFConvert)
{
	*myerror = NO_ERRORS;

    // Added SLEE
    m_MAFConvert = MAFConvert;
    std::string str;
   
    
	this->m_filename = filename;
	this->m_filename_bim = mapname;
	this->m_filename_fam = famname;
	this->m_filename_mwo = outfile;
	//================================
	std::string line;
	this->m_line_counter = -1;
	this->m_fam.open(this->m_filename_fam.c_str());
	if (!this->m_fam)
	{
		*myerror = CANT_OPEN_FAM_FILE4READ;
		return;
	}
	while (!this->m_fam.eof( ) ) 
	{
		getline(this->m_fam, line);
		this->m_line_counter++;
	}
	this->m_fam.close();
	//================================
	if (info == NULL)
	{
        str = filename;
        str += ".INFO.txt";
		this->m_info_file = str;
	}
	else
		this->m_info_file = info;

	this->m_info.open(this->m_info_file.c_str());
	if (!this->m_info)
	{
		*myerror = CANT_OPEN_INFO_FILE4WRITE;
		return;
	}

    str = this->m_info_file;
    str += ".REWR";
    
	this->m_info_rewritten = str;
	this->m_info_rewr.open(this->m_info_rewritten.c_str());
	if (!this->m_info_rewr)
	{
		*myerror = CANT_OPEN_INFO_RWR_FILE4WRITE;
		return;
	}

    str = this->m_info_file;
    str = this->m_info_file +".TEMP.txt";
	this->m_file_temp_name = str;

	this->m_snp_sets = ht->get_snps_sets(); 
	this->m_approx_line_lenght = ht->m_num_of_snps; //CHECK IT

	this->m_info << "-999" << "\tWindowSize" << std::endl;
	this->m_info << m_MAFConvert << "\tMAFConvert" << std::endl;
	this->m_info << ht->m_num_of_snps_insetid << "\tNumberOfSNPs" << std::endl;
	this->m_info << this->m_line_counter << "\tNumberOfIndividuals" << std::endl;

	this->m_info_rewr << "-999" << "\tWindowSize" << std::endl;
	this->m_info_rewr << m_MAFConvert << "\tMAFConvert" << std::endl;
	this->m_info_rewr << ht->m_num_of_snps_insetid << "\tNumberOfSNPs" << std::endl;
	this->m_info_rewr << this->m_line_counter << "\tNumberOfIndividuals" << std::endl;

	//================================
	this->m_file.open(this->m_filename.c_str(), std::ios::binary);
	if (!this->m_file)
	{
		*myerror = CANT_OPEN_BED_FILE4READ;
		return;
	}
	this->read_data_and_create_mwo_used_hashtable(ht,myerror);
	if (*myerror != 0)
		return;
	this->m_file.close();
	this->m_info_rewr.close();
	//================================

	int result0 = remove( this->m_file_temp_name.c_str() );
	if (result0 < -1 ) // "-1" - file not exists; "0" - successfully removed
	{
		*myerror = CANT_REMOVE_PREV_INFOTEMP_FILE;
		return;
	}
	int result1 = rename( this->m_info_file.c_str() , this->m_file_temp_name.c_str() );
	if (result1 != 0)
	{
		*myerror = CANT_RENAME_INFO2INFOTEMP_FILE;
		return;
	}
	int result2 = rename( this->m_info_rewritten.c_str() , this->m_info_file.c_str() );
	if (result2 != 0)
	{
		*myerror = CANT_RENAME_INFOREWRITTEN2INFO_FILE;
		return;
	}
}

//==========================================================
//This function creates *.mwa file for SeiId application
//Inputs:
//ht - lookup table object
//myerror - allocated memory to put there the final information 
//		if some error happen during the run.
//Explanation:
//This function reads "*.bed" in chunks of m_size_of_esi - exactly info for ONE SPECIFIC SNP.
//m_size_of_esi = (this->m_line_counter+3)/4;  // Number of bytes that hold all genotype info per one snp = one line length!!!
//moves inside of *.bed to reach specified location of specified snp - based on lookup table: ht->m_hash_table[j]
//then read from there exactly one line: "m_size_of_esi" bytes.  
//==========================================================
void BedFileReader::read_data_and_create_mwo_used_hashtable(Hasht* ht,int* myerror)
{
	//relevant info:
	//ht->m_hash_table;
	//ht->m_num_of_snps_insetid;
	//ht->m_setidf_setid;
	//m_snp_sets;   m_snp_sets->snp_id; 

	this->m_file_mwo.open(this->m_filename_mwo.c_str(), std::ios::out|std::ios::binary);
	if (!this->m_file_mwo)
	{
		*myerror = CANT_OPEN_MWA_FILE4WRITE;
		return;
	}
	this->m_file_mwo.seekp(std::ios::beg); 
	size_t ii;
	//int count_win_size = 1;
	int bits_val[MY_CHAR_BIT];
	int* temp_snp_info0 = new int[this->m_line_counter];
	int* temp_snp_info1 = new int[this->m_line_counter];
	for (ii = 0; ii < this->m_line_counter; ++ii)
	{
		temp_snp_info0[ii] = 0; 
		temp_snp_info1[ii] = 0; 
	}

	size_t individuals_counter = 0;
	char* encoded_snp_info = new char[(this->m_line_counter+3)/4] ;
	this->m_size_of_esi = (this->m_line_counter+3)/4;  // !!!Number of bytes per one snp = one line length!!!
	char* buff = new char [this->m_size_of_esi]; 
	char tmpbuf[3]; memset(tmpbuf, '\0', sizeof(tmpbuf));
	this->m_file.read(tmpbuf,sizeof(tmpbuf)); // three first bytes - permanent in bed file

	//=========OFFSET===================
	size_t set_counter = 1; //   Will be changed based on "m_setidf_setid" if it's change - new set!
	size_t begin, current; 
	begin = this->m_file_mwo.tellp();
	current = this->m_file_mwo.tellp();
	//this->m_info << "#=================================================#" << std::endl;
	this->m_info <<"SET#\tOFFSET\tSET_ID\tSET_SIZE" << std::endl;
	this->m_begin4rw = this->m_info.tellp();


	//=========================================
	// 
	//Go in loop "i = 0 ... ht->m_num_of_snps_insetid-1" toward "ht->m_setidf_setid" &&
	//based on ht->m_hash_table go by OFFSET in this->m_file , read relevant line: "this->m_size_of_esi" bytes ! 
	//put the readed bytes into "encoded_snp_info"
	//

	size_t setSize = 0;
	this->m_num_of_snps_insetid = ht->m_num_of_snps_insetid;
	for (size_t j = 0; j < ht->m_num_of_snps_insetid; ++ j)
	{

        /* added by SLEE */
        /*
        if(j % 10000 ==0){
            
            std::cout << j <<"/" << ht->m_num_of_snps_insetid << " read" <<'\xd' << std::flush ;

        }
        */
		if (j == 0)
		{
			//start new set
			this->m_info << set_counter << "\t" << (current - begin) << "\t" << ht->m_setidf_setid[j]; 
			set_counter ++;
		}
		else if (strcmp( ht->m_setidf_setid[j], ht->m_setidf_setid[j - 1]) != 0)
		{	//start new set

			this->m_file_mwo << std::endl; //ENTER - empty line between every two snp sets.
			current = this->m_file_mwo.tellp();
			this->m_info << "\t" << setSize << std::endl;
			this->m_info <<  set_counter << "\t" << (current - begin) << "\t" << ht->m_setidf_setid[j]; // << std::endl;
			set_counter ++;
			setSize = 0;
                        
            
		}
		buff[0] = '\0';
		//moving inside of *.bed to reach specified location of specified snp - based on lookup table: ht->m_hash_table[j]
		//then read from there exactly one line 
		this->m_file.seekg(this->m_size_of_esi * ht->m_hash_table[j] + 3 ,std::ios::beg); // +3 because of first 3 bytes in the file
		this->m_file.read(buff,this->m_size_of_esi); 
		setSize ++;

        std::string temp_id = this->m_snp_sets[ht->m_hash_table[j]].snp_id;
		this->m_file_mwo << temp_id << " ";


		for(size_t i=0; i<this->m_size_of_esi; i++)   //process the buff info
		{	//===============================================================					
			//=== This part converts Byte "buff[i]" to bits values "bits_val"
			//=== for example byte buff[0] = "w" ->  bits_val = 11101110
			memset(bits_val, NULL, sizeof(bits_val));
			int k = MY_CHAR_BIT;  //8
			while (k > 0)
			{
				-- k;
				bits_val[k] = (buff[i]&(1 << k) ? 1 : 0);
			}
			//==========================================================
			//=== here interpret Bit information "bits_val" to snps and count it - decode it
			decode_byte(bits_val, &individuals_counter, temp_snp_info0, temp_snp_info1, ht->m_hash_table[j]);

		} //end of for(int i=0; i<this->m_size_of_esi; i++)   //process the buff info

        // debug
        /*
        
        for(size_t i=0; i<this->m_line_counter; i++){
            
            printf("%d-%d:[%d][%d]\n", j,i,temp_snp_info0[i], temp_snp_info1[i]);
        }
        */
		//if (individuals_counter >= this->m_line_counter)
		{
			//Check who is MAGORITY/MINORITY
			//write data to output file from temp_snp_info0, temp_snp_info1)
			//encode temp_snp_info0, temp_snp_info1;

			// ENCODE OUTPUT
			//if (this->m_encode_output == 1)
			{
				memset(encoded_snp_info,0,sizeof(encoded_snp_info));
				
                //printf("%d:%d-%d\n",ht->m_hash_table[j],  m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[0], this->m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[1]);
                /* If m_MAFConvert==1 && 0 allele is the major allele */
                /* Very important!*/
                /* Be careful, total_counter_per_letter[0] < total_counter_per_letter[1] indicates 0 is the minor allele, so
                 we need to encode for 0 (add minor allele vector for encode) */
                
                if (m_MAFConvert == 1 && this->m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[0] < this->m_snp_sets[ht->m_hash_table[j]].total_counter_per_letter[1])
				{
					//write snp information as encoded
                    // Add minor allele 
					this->encode(temp_snp_info0,encoded_snp_info );
					for (ii = 0; ii < this->m_size_of_esi; ++ii)
					{
						this->m_file_mwo << encoded_snp_info [ii] ; 
					}
				} else
				{
					this->encode(temp_snp_info1,encoded_snp_info );
					for (ii = 0; ii < this->m_size_of_esi; ++ii)
						this->m_file_mwo << encoded_snp_info [ii] ; 
				}
			}// END OF ENCODE OUTPUT

			this->m_file_mwo << std::endl; //ENTER at the end of every line
			individuals_counter = 0;
			for (ii = 0; ii < this->m_line_counter; ++ii)
			{
				temp_snp_info0[ii] = 0; 
				temp_snp_info1[ii] = 0; 
			}
			//snp_set_ind ++;
		}  //end of if (individuals_counter >= this->m_line_counter)
	} // end of for (int j = 0; j < ht->m_num_of_snps_insetid; ++ j)

	delete[] temp_snp_info0;
	delete[] temp_snp_info1;
	delete[] encoded_snp_info;

	this->m_file_mwo << std::endl; //ENTER at the end of every line
	this->m_file_mwo << '\0'; 
	this->m_file_mwo.close();

	this->m_info << "\t" << setSize  << std::endl;  // print last set size
	this->m_info << "#=================================================#" << std::endl;
	this->m_info << this->m_size_of_esi + 1 << 
		"\tDECODED NumberOfDECODEDbytesPerSNP(NotIncludesSnpID&SpaceAfter_Includes\\n) " << std::endl;
	this->m_info << set_counter - 1 << "\tTotalNumberOfSets" << std::endl;

	this->m_info_rewr << this->m_size_of_esi + 1 << 
		"\tDECODED NumberOfDECODEDbytesPerSNP(NotIncludesSnpID&SpaceAfter_Includes\\n) " << std::endl;
	this->m_info_rewr << set_counter - 1 << "\tTotalNumberOfSets" << std::endl;


	m_set_counter = set_counter - 1;

	this->m_info << '\0';
	this->m_info.close();

	//==========================================================
	//  REWRITE INFO FILE
	//==========================================================

	this->m_infoi.open(this->m_info_file.c_str());
	if (!this->m_infoi)
	{
		*myerror = CANT_OPEN_INFO_FILE4READ;
		return;
	}

	this->m_infoi.seekg (this->m_begin4rw, std::ios::beg);
	//this->m_info_rewr << "#=================================================#" << std::endl;
	this->m_info_rewr <<"SET#\tOFFSET\tSET_ID\tSET_SIZE" << std::endl;

	std::string line;
	size_t kk = 0;
	while (kk < m_set_counter ) 
	{
		getline(this->m_infoi, line);
		kk ++;
		this->m_info_rewr << line << std::endl;
	}
	this->m_infoi.close();
	
}

//==========================================================
// This function converting(encoding) every 4 bytes of temp_snp_info to "a" 
// "a"- array of 8 integers ("kind of" 8 bits) that will be converted to one encoded number
// possible values for temp_snp_info = 9,1,2,0
//Inputs:
//temp_snp_info - info of current snp
//encoded_snp_info - encoded info of current snp, this info will be written to "*.mwa" file
//==========================================================
void BedFileReader::encode(int* temp_snp_info,char* encoded_snp_info )
{

	size_t i = 0;
	unsigned int j = 0;
	int number = 0;
	size_t ind4enc = 0;
	int a[8]; 

	//=======================================================================================
	// converting every 4 bytes of temp_snp_info to "a" 
	// "a"- array of 8 integers ("kind of" 8 bits) that will be converted to one encoded number
	// possible values for temp_snp_info = 9,1,2,0
	while (i < this->m_line_counter)
	{
		memset(a, 0, sizeof(a)); 
		for (j = 0; j < 4; ++j)
		{	
			if (temp_snp_info[i] == 9) //missing value
			{
				a[j*2] = 1;
				a[j*2+1] = 0;
			}
			else if (temp_snp_info[i] == 1)//AG
			{
				a[j*2] = 0;
				a[j*2+1] = 1;
			}
			else if (temp_snp_info[i] == 0)//GG
			{
				a[j*2] = 0;
				a[j*2+1] = 0;
			}
			else if (temp_snp_info[i] == 2)//AA
			{
				a[j*2] = 1;
				a[j*2+1] = 1;
			}
			i ++;
		}

		ind4enc ++;
		if (ind4enc == this->m_size_of_esi+1)
			break; 
		else
		{	
			number = 0;
			//==============================================
			//converting 8 ints of "a" to one encoded number
			//for example a= 00100010 "number" will be: 34
			for (int  ii = 0; ii < 8; ++ii)
				number += a[ii] * (int)pow(2.0,(7-ii));
			//saving this encoded number to array that will be written into *.mwa file.
			encoded_snp_info[ind4enc-1] = (char)number;
			//=============================================
		}
	}
}

//==========================================================
// This function interpret Bit information "bits_val" to snps 
// and count it - decode it
// Inputs:
// bits_val	- bits representation of some character
//		for example bits_val = 01000100, read in couples
// individuals_counter - number of individuals that already counted in current snp line
// temp_snp_info0 - vector of integers that holds counting for first character 
//					of genotype(for example for this snp in use "G" and "A")
//					temp_snp_info0 will count number of "G") per individual 
// temp_snp_info1 - vector of integers that holds counting for second character 
//					of genotype(for example for this snp in use "G" and "A")
//					temp_snp_info0 will count number of "A") per individual
// snp_set_ind  - number of current set

//For example:
//					bits_val[0] = 0 and bits_val[1] = 0
//					temp_snp_info0[*individuals_counter - 1] = 2;
//					temp_snp_info1[*individuals_counter - 1] = 0;
//					//Homozegote 1 for example GG      write 20  (2 of "G", 0 of another character)

//					bits_val[0] = 1 and bits_val[1] = 1
//					temp_snp_info0[*individuals_counter - 1] = 0;
//					temp_snp_info1[*individuals_counter - 1] = 2;
//					//Homozegote 2 for example AA      write 02 (0 of "G", 2 of "A")

//					bits_val[0] = 1 and bits_val[1] = 0
//					temp_snp_info0[*individuals_counter - 1] = 9;
//					temp_snp_info1[*individuals_counter - 1] = 9;
//					//Missing value

//					bits_val[0] = 0 and bits_val[1] = 1
//					temp_snp_info0[*individuals_counter - 1] = 1;
//					temp_snp_info1[*individuals_counter - 1] = 1;
//					////Heterozegote for example AG  or GA     write 11 (1 of "G", 1 of "A")

						
//==========================================================
void BedFileReader::decode_byte(int* bits_val,size_t * individuals_counter, 
								int* temp_snp_info0, int* temp_snp_info1, size_t snp_set_ind)
{


	int flag = 0;
	for (int i = 0; i < 4; ++i)
	{
		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;
			if (*individuals_counter > this->m_line_counter)
				return;
			this->m_snp_sets[snp_set_ind].total_counter_per_letter[0] += 2;
			temp_snp_info0[*individuals_counter - 1] = 2;
			temp_snp_info1[*individuals_counter - 1] = 0;
			flag = 1 ;//Homozegote 1 for example GG      write 20 ; 00 - will point to [0] letter   + 2 to [0]

		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > this->m_line_counter)
				return;
			this->m_snp_sets[snp_set_ind].total_counter_per_letter[1] += 2;
			temp_snp_info0[*individuals_counter - 1] = 0;
			temp_snp_info1[*individuals_counter - 1] = 2;
			flag = 2 ;//Homozegote 2 for example AA      write 02 ; 11 - will point to [1] letter   + 2 to [1]
		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;		
			if (*individuals_counter > this->m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 9;
			temp_snp_info1[*individuals_counter - 1] = 9;
			flag = 3 ; //Missing value                   nothing to add - write 99 ;
		}
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > this->m_line_counter)
				return;
			this->m_snp_sets[snp_set_ind].total_counter_per_letter[0] ++;
			this->m_snp_sets[snp_set_ind].total_counter_per_letter[1] ++;
			temp_snp_info0[*individuals_counter - 1] = 1;
			temp_snp_info1[*individuals_counter - 1] = 1;

			flag = 4 ; //Heterozegote for example AG  or GA     write 11 ; 01 - will point to [0] and [1] letter   +1 +1 to [1]
		}
		else
			flag = 5 ; //Error
	}

    //printf("Dec[%d]%d-%d\n",snp_set_ind,  m_snp_sets[snp_set_ind].total_counter_per_letter[0], m_snp_sets[snp_set_ind].total_counter_per_letter[1]);
}

//==========================================================
//Destructor - free all dynamically allocated memory
//==========================================================
BedFileReader::~BedFileReader()
{
	delete[] m_snp_sets;
}

//==========================================================
//This function in use just in case of Moving Window
//1)counts lines in bim file - 
//  to know how many snps exists for future space allocation
//2)counts lines in fam file-
//  to know how many individuals exists for future space allocation

//Inputs:
//bim_file - "*.bim" file
//fam_file - "*.fam" file
//myerror - allocated memory to put there the final information 
//		if some error happen during the run.
//==========================================================
void BedFileReader::init(const char* bim_file, const char* fam_file,int* myerror)
{
	this->m_filename_bim = bim_file;
	//read mapfile, fill this->m_snp_sets[ii].snp_id; ii=0..this->m_approx_line_lenght
	std::string line;
	this->m_approx_line_lenght = -1;//0;
	this->m_bim.open(this->m_filename_bim.c_str());
	if (!this->m_bim)
	{
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}
	while (!this->m_bim.eof( ) ) 
	{
		getline(this->m_bim, line);
		this->m_approx_line_lenght++;
	}
	this->m_bim.close();


	this->m_filename_fam = fam_file;
	this->m_line_counter = -1;
	this->m_fam.open(this->m_filename_fam.c_str());
	if (!this->m_fam)
	{
		*myerror = CANT_OPEN_FAM_FILE4READ;
		return;
	}
	while (!this->m_fam.eof( ) ) 
	{
		getline(this->m_fam, line);
		this->m_line_counter++;
	}
	this->m_fam.close();
}

//=======================================================================
// This function reads "*.bim" file,
// creates "m_snp_sets" array of SNP_info objects that will hold info during future "*.mwa" preparation
// of how many genotypes of every type per snp - to calculate minor and major
// in use just in case of "Moving window", 
// otherwise this info comes from setid_bim_index.cpp 
// Inputs:
// myerror - allocated memory to put there the final information 
//		if some error happen during the run.
//=======================================================================
void BedFileReader::upload_snpid_from_bim(int* myerror)
{	
	std::string line;
	size_t ii = 0;
	this->m_bim.open(this->m_filename_bim.c_str());
	this->m_bim.seekg (0, std::ios::beg);
	if (!this->m_bim)
	{
		*myerror = CANT_OPEN_BIM_FILE4READ;
		return;
	}

	while (!this->m_bim.eof( ) ) 
	{
		getline(this->m_bim, line);
		for (size_t i = 0; i < line.size(); ++i)
		{
			if (line.at(i) == 9 || line.at(i) == ' ' || line.at(i) == ',' || line.at(i) == '\t')
			{	size_t j = i + 1;
                size_t k = 0;
			while(1)
			{
				this->m_snp_sets[ii].snp_id[k] = line.at(j);
				k++;
				j++;
				if (line.at(j)== 9 || line.at(j) == ' ' || line.at(j) == ',' || line.at(j) == '\t')
				{
					this->m_snp_sets[ii].snp_id[k] = '\0';
					break;
				}
			}
            
            // Added by SLEE
            // currently there is no need to save allele
			//this->m_snp_sets[ii].letters[1] = line.at(line.size()-1);
			//this->m_snp_sets[ii].letters[0] = line.at(line.size()-3);
			break;
			}
		}
		ii++;
	}


	this->m_bim.close();

}
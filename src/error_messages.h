/*************************************************************
 *
 * Moving window and Set-Bed merge. Sets creating project 
 * File: error_meassages.h	
 * Date: Jan 19, 2011
 * Author: Larissa Miropolsky
 *
 * Description:
 *   Definition of different types of known errors
 *
 **************************************************************/
#ifndef _ERROR_MESSAGES_H        
#define _ERROR_MESSAGES_H

// Added: no more size limitation
#define SNP_ID_SIZE_MAX 1024


#define NO_ERRORS					0			//all write
#define CANT_OPEN_BIM_FILE4READ		1			//check if it's exists
#define CANT_OPEN_FAM_FILE4READ		2			//check if it's exists
#define CANT_OPEN_BED_FILE4READ		3			//check if it's exists
#define CANT_OPEN_SETID_FILE4READ	4			//check if it's exists
#define CANT_OPEN_MWA_FILE4WRITE	5			//ckeck permissions
#define CANT_OPEN_MWA_FILE4READ		6			//check if it's exists
#define CANT_OPEN_INFO_FILE4WRITE	7			//ckeck permissions
#define CANT_OPEN_INFO_FILE4READ	8			//check if it's exists
#define CANT_OPEN_INFO_RWR_FILE4WRITE	9			//ckeck permissions
#define CANT_REMOVE_PREV_INFOTEMP_FILE		10  	//check if it's open, check permissions
#define CANT_RENAME_INFO2INFOTEMP_FILE		11  //check if it's open, check permissions
#define CANT_RENAME_INFOREWRITTEN2INFO_FILE 12 	//check if it's open, check permissions
#define WRONG_ALLOCATED_SIZE4OUTPUT_ARRAY	13  //check Number of individuals and Number of SNPs in for current set
#define REQUESTED_SET_ID_DOESNOT_EXISTS		14  //check setID in info file
#define WARNING_CANT_OPEN_FILE4WRITE_2PRINTSET		15 // warning - if can't print - skip printing
#define BIM_FILE_TOKEN_WRONG 16 
#define BED_FILE_READ_ERROR 17 

#endif //_ERROR_MESSAGES_H





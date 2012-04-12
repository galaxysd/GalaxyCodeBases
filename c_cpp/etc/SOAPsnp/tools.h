/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-14
 *CLASS NAME: 
 *FUNCTION : some useful function to the SOAPsnp
 *FILE NAME : tool.h
 *UPDATE DATE : 2010-8-3
 *UPDATE BY : Bill Tang
 *UPDATE DATE : 2010-9-6
 *UPDATE BY :Zhou Guyue
 *UPDATE: 2010-10-15 change the struct SFS_PARA, add a member qs.
 *******************************************************************************
 */
#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <strstream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
using namespace std;
typedef ifstream my_ifstream;

#define NOUSE_ALIGNMENT "alignment is no use"
#define HIT_COUNT_X0		"X0:i:"
#define HIT_COUNT_H0		"H0:i:"
#define SFS_PARA_ERROR   -1              //if sfs para is wrong return SFS_PARA_ERROR
#define SFS_PARA_SUCC    1				 //if sfs para is right return SFS_PARA_SUCC


/* integer to string*/
std::string myitoa(int num);
// A function to spilt string s into vector vec according to char splitchar
void StringSplit(std::string s, char splitchar, std::vector<std::string>& vec);
/* count the indel length from the cigar string. return the indel length and the position*/
int count_indel_len(const std::string cigar, int &pos);
/* format the sam text to the soap text*/
std::string alignment_format(const std::string &sam_ali);

/*a struct for sfs parament*/
typedef struct
{
	int doBay;			// 0:not Bay£¨1:do Bay£¨default is 0
	int doJoint;			//0:not Joint£¨1:do Joint£¨default is 0
	int writeFr; 		// 0:not write Fr, 1: write Fr, default is 0
	std::string outfiles;		//out files path
	int start_pos;	// start position
	int stop_pos;	// stop position
	int under_FP;		// 0: not under Flow Protect £¨1:under Flow Protect£¨default is 0
	int allow_PathZ;		// 0: notallow Phat Zero£¨1:allow Phat Zero «£¨default is 0
	double bias;	// default bias is 0.03
	double pvar;	// The probability for beeing variable in the population£¨default is 0.015
	double eps;	//error rate  ,default is 0.005
	double pThres ;		//select the out put threshold of snp £¨default is 0.01
	int sigle_MJ;	 // 0:not sigle minor joint£¨1:sigle minor joint£¨default is 0
	int sigle_MB;	 //0: not sigle minor bay£¨1: sigle minor bay£¨default is  0
	int alternative; //
	int jointFold;	//
	int qs;	// the quality score threshold.
	std::string faiPath; // fai file's path.
	//update 11-29
	int sfs_first_last; // 0: output h[0]~h[2k] , 1:output h[0] and h[2k] 
}SFS_PARA ;

/*initialize a sfs_para struct*/
#define sfs_init()(new SFS_PARA)
/*destroy a sfs_para struct*/
#define sfs_destroy(b) do {if (b) { delete b; }} while (0)	
/*sfs help information*/
void sfs_info(void);
/* initialize sfs parament*/
int init_sfs_para(const string sfs_path, SFS_PARA *sfs );
/*change string to bool*/
int str_bool( const char* str);
/* count soft clip number at the beginning*/
int count_soft_clip(const std::string cigar, int &last_S);
 
#endif

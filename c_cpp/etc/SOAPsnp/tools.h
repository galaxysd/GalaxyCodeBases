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
 *******************************************************************************
 */
#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <strstream>
#include <vector>
#include <string>
#include <sstream>


#define NOUSE_ALIGNMENT "alignment is no use"
#define HIT_COUNT_X0		"X0:i:"
#define HIT_COUNT_H0		"H0:i:"

/* integer to string*/
std::string myitoa(int num);
// A function to spilt string s into vector vec according to char splitchar
void StringSplit(std::string s, char splitchar, std::vector<std::string>& vec);
/* count the indel length from the cigar string. return the indel length and the position*/
int count_indel_len(const std::string cigar, int &pos);
/* format the sam text to the soap text*/
std::string alignment_format(const std::string &sam_ali);
#endif

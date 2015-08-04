/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-14
 *CLASS NAME:
 *FUNCTION : definition of the useful functions
 *FILE NAME : tool.cpp
 *UPDATE DATE : 2010-8-3
 *UPDATE BY : Bill Tang
 UPDATE DATE : 2010-9-9
 *UPDATE BY : Bill Tang
  UPDATE DATE : 2010-9-15
 *UPDATE BY : Zhou Guyue
 *******************************************************************************
 */

#include <string.h>
#include <stdlib.h>
#include "tools.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

/* integer to string*/
std::string myitoa(int num) {
	char num_str[256] = {0};
	sprintf(num_str, "%d", num);
	return num_str;
}

// A function to spilt string s into vector vec according to char splitchar
void StringSplit(std::string s, char splitchar, std::vector<std::string>& vec) {
	// promise the vec is empty
	if(vec.size()>0)
		vec.clear();
	// chomp the string
	while(s.size() > 0 && s[0] == splitchar)
		s = s.substr(1, s.length() - 1);
	while(s.size() > 0 && s[s.length() - 1] == splitchar)
		s = s.substr(0, s.length() - 1);
	int length = s.length();
	int start=0;
	for(int i=0;i<length;i++)
	{
		if(s[i] == splitchar)
		{
			vec.push_back(s.substr(start,i - start));
			while(s[i + 1] == splitchar)
				i ++;
			start = i+1;
		}
		else if(i == length-1)// attach last
		{
			vec.push_back(s.substr(start,i+1 - start));
		}
	}
}

/* count the indel length from the cigar string. return the indel length and the position*/
int count_indel_len(const std::string cigar, int &pos) {
	const char *tmp = cigar.c_str();
	int begin = 0;
	int end;
	int len = 0;
	pos = 0;
	int indel_flag = 0;
	for (end = 0; end < cigar.size(); end ++) {
		switch (tmp[end]) {
			case 'M':
			{
				pos += atoi(cigar.substr(begin, end - begin).c_str());
				begin = end + 1;
				break;
			}
			case 'I':
			{
				if (indel_flag != 0)
				{
					return 300;
				}
				//len = 200 + atoi(cigar.substr(begin, end - begin).c_str());  //for new soap format
				len = 100 + atoi(cigar.substr(begin, end - begin).c_str()); //for old soap format
				begin = end + 1;
				indel_flag = 1;
				break;
			}
			case 'D':
			{
				if (indel_flag != 0)
				{
					return 300;
				}
				//len = 100 + atoi(cigar.substr(begin, end - begin).c_str());  //for new soap format
				len = 200 + atoi(cigar.substr(begin, end - begin).c_str()); //for old soap format
				begin = end + 1;
				indel_flag = 1;
				break;
			}
			case 'N':
			{
				pos += atoi(cigar.substr(begin, end - begin).c_str());
				begin = end + 1;
				break;
			}
			case 'S':
			{
				pos += atoi(cigar.substr(begin, end - begin).c_str());
				begin = end + 1;
				break;
			}
			case 'H':
			{
				pos += atoi(cigar.substr(begin, end - begin).c_str());
				begin = end + 1;
				break;
			}
			case 'P':
			{
				pos += atoi(cigar.substr(begin, end - begin).c_str());
				begin = end + 1;
				break;
			}
			default : break;
		}
		if (len != 0) {
			break;
		}
	}
	return len;
}

/* format the sam text to the soap text*/
std::string alignment_format(const std::string &sam_ali) {
	if(sam_ali.empty()) {
		return NOUSE_ALIGNMENT;
	}
	std::vector<std::string> vec;
	StringSplit(sam_ali, '\t', vec);
	if ((vec.size() < 11) || (vec[5] == "*"))
		return NOUSE_ALIGNMENT;

	std::string format;
	int pos = 0;
	int best_hit = 1;
	int flag = atoi(vec[1].c_str());
	int last_clip_num;
	int clip_num = count_soft_clip(vec[5], last_clip_num);

	if (clip_num + last_clip_num > vec[9].size())
		return NOUSE_ALIGNMENT;

	if (flag & (0x1 << 6)) {
		format = vec[0] + "/1" + "\t"; // add query name
	} else if (flag & (0x1 << 7)) {
		format = vec[0] + "/2" + "\t"; // add query name
	} else {
		format = vec[0] + "\t";  // add query name
	}

	format += vec[9].substr(clip_num, vec[9].size() - last_clip_num - clip_num) + "\t"; // add sequence
	format += vec[10].substr(clip_num, vec[10].size() - last_clip_num - clip_num) + "\t"; // add quality

	for (int i = 11; i < vec.size(); i++) {
		if (vec[i].find(HIT_COUNT_H0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		} else if (vec[i].find(HIT_COUNT_X0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		}
	}
	format += myitoa(best_hit)  + "\t"; // add number of best hit.

	format += ((flag & (0x1 << 6)) ? "a" : "b"); // add a/b
	format += "\t";
	format += myitoa(vec[9].size() - clip_num - last_clip_num); // add length
	format += "\t";
	format += ((flag & (0x1 << 4)) ? "-" : "+");  // add strand
	format += "\t";
	format += vec[2] + "\t";  // add chr
	format += vec[3] + "\t";  // add location
	int mismatch = count_indel_len(vec[5], pos);
	if (mismatch == 300)
	{
		return NOUSE_ALIGNMENT;
	}
	format += myitoa(mismatch);  // add number of mismatch
	format += "\t";						//for old soap format
	format += myitoa(pos - clip_num);  // add indel position
	//format += myitoa(pos) + "\t";  // add indel position. changed at 2010-11-22, to compate to the new version soap.
	//format += vec[5];               //for new soap format
	return format;
}

/**
 * DATE: 2010-9-9
 * FUNCTION: count the soft clip number at the read's beginning.
 * PARAMETER: cigar: the cigar string. last_S: last soft clip's length.
 * RETURN:	soft clip number.
 */
int count_soft_clip(const std::string cigar, int &last_S)
{
	const char *tmp = cigar.c_str();
	last_S = 0;
	// count the last soft clip number.
	if (tmp[cigar.size() - 1] == 'S')
	{
		int j = cigar.size() - 2;
		while (tmp[j] >= '0' && tmp[j] <= '9')
		{
			j--;
		}
		last_S = atoi(cigar.substr(j + 1, cigar.size() - j - 1).c_str());
	}

	int i = 0;
	for (; i < cigar.size(); ++i)
	{
		if (tmp[i] >= '0' && tmp[i] <= '9')
		{
			// jump to the next position.
			continue;
		}
		else if (tmp[i] == 'S')
		{
			// return the soft clip number.
			int num = atoi(cigar.substr(0, i).c_str());
			return num;
		}
		else
		{
			// no soft clip infront of the reads.
			break;
		}
	}
	return 0;
}

/**
 * DATE: 2010-9-7
 * FUNCTION: change string to bool
 * PARAMETER: str: a string which will change to bool. len: the str length.
 * RETURN:	bool type of str
 */
int str_bool(const char* str)
{
	if ( *str == '0')
	{
		 return 0;
	}
	else
	{
			return 1;
	}
}

/**
 * DATE: 2010-9-20
 * FUNCTION: intialize sfs parament
 * PARAMETER: sfs_path: sfs configuration file path. sfs: a sfs_para struct.
 * RETURN:	 SFS_PARA_ERROR:if sfs parament is wrong ,SFS_PARA_SUCC : if sfs parament is right.
 * UPDATE: 2010-10-15
 */
int init_sfs_para(const string sfs_path, SFS_PARA  *sfs )
{
	string in;

	/*set default vaule for sfs_par struct*/
	sfs->allow_PathZ = 0;

	sfs->bias = 0.03;
	sfs->doBay = 0;
	sfs->doJoint = 0;
	sfs->eps = 0.005;
	sfs->outfiles = "";
	sfs->pThres = 0.01;
	sfs->pvar = 0.015;
	sfs->sigle_MB = 0;
	sfs->sigle_MJ = 0;
	sfs->start_pos = 0;
	sfs->stop_pos = 0;
	sfs->under_FP = 0;
	sfs->writeFr = 0;
	sfs->alternative = 0;
	sfs->jointFold = 1;
	sfs->qs = 52;
	//update 11-29 for control the column
	sfs->sfs_first_last = 0 ;

	my_ifstream sfsfile(sfs_path.c_str()); //open sfs configuration file

	if (!sfsfile.is_open())         //can not open sfsfile
	{
		cerr << "Cannot open file:" << sfs_path << endl;
		return SFS_PARA_ERROR;
	}

	stringstream temp_ss;
	while (sfsfile >> in)
	{
		if (in.empty())
		{
			continue;
		}
		if (in == "doBay")
		{
			sfsfile >> in ;
			sfs->doBay = str_bool(in.c_str());
			continue;
		}
		else if (in == "doJoint")
		{
			sfsfile >> in ;
			sfs->doJoint = str_bool(in.c_str());
			continue;
		}
		else if (in == "writeFr")
		{
			sfsfile >> in ;
			sfs->writeFr = str_bool(in.c_str());
			continue;
		}
		else if (strcmp(in.c_str() , "outfiles") == 0)
		{
			sfsfile >> sfs->outfiles;
			continue;
		}
		else if (in == "start")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->start_pos;
			temp_ss.clear();
			continue;
		}
		else if (in == "stop")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->stop_pos;
			temp_ss.clear();
			continue;
		}
		else if (in == "underFlowProtect")
		{
			sfsfile >> in ;
			sfs->under_FP = str_bool(in.c_str());
			continue;
		}
		else if (in == "allowPhatZero")
		{
			sfsfile >> in ;
			sfs->allow_PathZ = str_bool(in.c_str());
			continue;
		}
		else if ( in == "bias")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->bias;
			temp_ss.clear();
			continue;
		}
		else if (in == "pvar")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->pvar;
			temp_ss.clear();
			continue;
		}
		else if (in == "eps")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->eps;
			temp_ss.clear();
			continue;
		}
		else if (in == "pThres")
		{
			sfsfile >> in ;
			temp_ss << in;
			temp_ss >> sfs->pThres;
			temp_ss.clear();
			continue;
		}
		else  if (in == "sigleMinorJoint")
		{
			sfsfile >> in ;
			sfs->sigle_MJ = str_bool(in.c_str());
			continue;
		}
		else if (in == "sigleMinorBay" )
		{
			sfsfile >> in ;
			sfs->sigle_MB = str_bool(in.c_str());
			continue;
		}
		else if (in == "alt" )
		{
			sfsfile >> in ;
			sfs->alternative = atoi(in.c_str());
			continue;
		}
		else if (in == "qs" )
		{
			sfsfile >> sfs->qs;
			continue;
		}
		//else if (in == "sfsfirstlast")
		else if ( in == "lightoutput")
		{
			sfsfile >> in ;
			sfs->sfs_first_last = str_bool(in.c_str());
			continue;
		}
		else
		{
			sfsfile.close();
			return SFS_PARA_ERROR;
		}
	}

	sfsfile.close();
	return SFS_PARA_SUCC;
}

/*sfs help information*/
/**
 * DATE: 2010-9-7
 * FUNCTION: return help information
 * PARAMETER: void
 * RETURN:	 void
 */
void sfs_info()
{
	cerr << "\t-> (outfiles) -outfiles \n" << endl;
	cerr << "\t-> (which analysis) -doBay -doJoint\n" << endl;
	cerr << "\t-> (optional variables)  -eps  -bias  -pThres -underFlowProtect\n" << endl;
	cerr << "\t-> (optional variables) -start -stop \t(format is chromo:position)\n" << endl;
	cerr << "\t-> (optional variables) -singleMinor[Joint,Bay] [0 or 1] should we only use a singleminor?\n" << endl;
	cerr << "\t-> (optional variables) -sfsfirstlast [0 or 1] should we only output sfs the first column and the last column \n" <<endl;
}

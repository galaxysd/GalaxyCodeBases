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
 *******************************************************************************
 */

#include "tools.h"


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
				len = 100 + atoi(cigar.substr(begin, end - begin).c_str()); 
				begin = end + 1;
				break;
			}
			case 'D':
			{
				len = 200 + atoi(cigar.substr(begin, end - begin).c_str()); 
				begin = end +1;
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
	char num[10];
	int pos = 0;
	int best_hit = 1;
	int flag = atoi(vec[1].c_str());
	if (flag & (0x1 << 6)) {
		format = vec[0] + "/1" + "\t"; // add query name
	} else if (flag & (0x1 << 7)) {
		format = vec[0] + "/2" + "\t"; // add query name
	} else {
		format = vec[0] + "\t";  // add query name
	}
	format += vec[9] + "\t"; // add sequence
	format += vec[10] + "\t"; // add quality
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
	format += myitoa(vec[9].size()); // add length
	format += "\t";
	format += ((flag & (0x1 << 4)) ? "-" : "+");  // add strand
	format += "\t";
	format += vec[2] + "\t";  // add chr
	format += vec[3] + "\t";  // add location
	format += myitoa(count_indel_len(vec[5], pos));  // add number of mismatch
	format += "\t";
	format += myitoa(pos);  // add indel position
	return format;
}

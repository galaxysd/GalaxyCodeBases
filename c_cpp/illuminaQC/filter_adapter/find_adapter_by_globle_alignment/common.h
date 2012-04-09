//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22

#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

using namespace std;


void split(string &strLine, vector<string>& tokens, const char* delim);
void numberStatic(vector<float> &vecNum, float &mean, float &median, float &SD, float &max, float &min);
string int2str( int InNum);
string float2str( float InNum);

#endif


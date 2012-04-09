//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22

#ifndef _LOCATE_ADAPTER_H_
#define _LOCATE_ADAPTER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

using namespace std;


int locate_adapter ( string &seq, string adapter, int insert_size, int &align_len, int &mis_match );

#endif


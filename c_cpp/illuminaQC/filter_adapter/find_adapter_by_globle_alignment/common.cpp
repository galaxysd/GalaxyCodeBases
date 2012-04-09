//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22

#include "common.h"

using namespace std;


//simulate the perl "split" function, split a string into an array
void split(string &strLine, vector<string>& tokens, const char* delim) 
{
	int count = 0;
	for(;;) {
		//erase delimiter
	   int i = strLine.find_first_not_of(delim);
	   if(i == -1)
	    	break;
	   strLine.erase(0, i);
	   
	   i = strLine.find_first_of(delim);
	   if(i == -1) {
		    tokens.push_back(strLine);
		    break;
	   } else {
		    string token = strLine.substr(0, i);
		    strLine.erase(0, i);
		    tokens.push_back(token);
	   }
	}
}

//caculate mean median SD, min, max for a given set of numbers 
void numberStatic(vector<float> &vecInput, float &mean, float &median, float &SD, float &max, float &min)
{	
	vector<float> vecNum = vecInput;
	int N = vecNum.size();
	float sum = 0.0, sum_2 = 0.0;
	sort( vecNum.begin(), vecNum.end() );
	min = vecNum[0];
	max = vecNum[N-1];
	median = vecNum[int(N/2)];

	for (int i=0; i<N; i++)
	{	sum += vecNum[i];
		sum_2 += vecNum[i] * vecNum[i];
	}
	
	mean = sum / N;
	SD = sqrt( (sum_2 - sum * sum / N) / (N - 1) );
}


//convert int number to string 
string int2str( int InNum)
{	string OutStr;
	char MidStr[1000]; //最大为1000位的数
	sprintf(MidStr, "%d",InNum);
	OutStr = MidStr;
	return OutStr;
}

//convert float number into string, 保留六位小数
string float2str( float InNum)
{	string OutStr;
	char MidStr[1000]; //最大为1000位的数
	sprintf(MidStr, "%f",InNum);
	OutStr = MidStr;
	return OutStr;
}


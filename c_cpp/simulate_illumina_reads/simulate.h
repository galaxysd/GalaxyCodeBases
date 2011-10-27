#ifndef __SIMULATE_H_
#define __SIMULATE_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>

using namespace std;

extern char alphabet[128];
extern char Bases[5];
extern char c_bases[5];
extern double GC_bias_abundant[101];

//for split string
void split(string strLine, vector<string>& tokens, const char* delim);

//simulate normal distribution insertsize by Box-muller method
int simulate_insertsize(int mean, int sd);

//simulate quality value
string simulate_quality(vector<int> error_pos, int read_len, double** correct_base_quality, double** error_base_quality);

//simulate GC bias
int simulate_GC_bias(string insert_str);

//Simulate illumina error distribution on different cycles,
//with the model function f(x)=0.00001*x**4
vector <double> error_distribution(float error_rate, int read_length);

//get the reverse and complement sequence
string reversecomplementary (string read);

//Realization of error sequencing
char get_error_match(char base);

//Realization of snp 
char get_match(char base);

//Produce heterozygous SNPs in multiploid
string Get_snp(string &seq,ofstream &snp,string id, float hetersnp_rate);

//Getting insert sequence when heterozygous indel rate is bigger than 0
string get_insertion(int num);

//Produce heterozygous indels in multiploid
string Get_indel(string &seq,ofstream &indel,string id1,float heterindel_rate);

#endif

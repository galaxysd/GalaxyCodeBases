//filter illumina raw reads from two directions:
//(1) filter the low quality data, by using the Ns,Bs,and quality information.
//    rank the reads by average error rate.
//(2) filter the major parts of adapters, by alignment and position judgement.
//    use dynamic programming alignment (non-gapped) to find adapters, slow but sensitive. 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"

using namespace std;

int Is_filter_quality = 0;
double Error_rate_cutoff = 0.03;
int Min_read_len = 50;
int Qual_shift = 64;
double Qual2Err[128];
uint64_t *StatRank;
uint64_t AdptInsert[250]; //insert <= read_length
uint64_t Nrate[250]; //N number <= read length

int Is_filter_adapter = 0;
string adpt1="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";//pDNA-3+
int Given_read_length = 100;
int *DPscore;
int Adapter_align_cutoff = 50;  //minimum alignment score 

int Stat_only_mode = 0;

char alphabet[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//
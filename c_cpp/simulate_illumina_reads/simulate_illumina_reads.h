#ifndef __SIMULATE_ILLUMINA_READS_H_
#define __SIMULATE_ILLUMINA_READS_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include "simulate.h"
#include "gzstream.h"

using namespace std;

//read in quality distribution file and initialize the quality distribution table.
void get_quality_distr(string exe_path);

//get genome sequence then simulate reads
void Get_genome(igzstream &inf,igzstream &inf2,ofstream &log1);

//get fragment from genome sequence and simulate reads 
long long get_reads(string id_line,string id,string &sequ,string &sequ2,ofstream &log2,long long read_genome);

//output simulated reads
long long output_reads(vector <double> err_dist,string &seq,int seqlen, int rd_pair, string id_seq,ofstream &log3,long long reads_all);

#endif

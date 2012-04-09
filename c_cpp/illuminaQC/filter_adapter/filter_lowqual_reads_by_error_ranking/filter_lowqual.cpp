// filter the low quality data, by using the Ns,Bs,and quality information.
// rank the reads by average error rate.

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
#include "gzstream.h"

using namespace std;

double Error_rate_cutoff = 0.03;
int Min_read_len = 50;
int Quality_shift = 64;
int Given_read_length = 100;
int Is_trimming = 0;
int Stat_only_mode = 0;
double Qual2Err[128];
uint64_t *StatRank;

void usage() 
{	cout << "filter_lowqual_by_rank <reads_x.fq>  <reads_x.fq.clean>  <result_stat_file>" << endl;
	cout << "   -e <float>  set the error rate cutoff for trimming, default=" << Error_rate_cutoff << endl;
	cout << "   -m <int>    set the minimum trimmed read length, default=" << Min_read_len << endl;
	cout << "   -q <int>    set the quality shift value, default=" << Quality_shift << endl;
	cout << "   -r <int>    set the sequencing read length, default=" << Given_read_length << endl;
	cout << "   -t <int>    whether trim, 0:no; 1:yes; default=" << Is_trimming << endl;
	cout << "   -s <int>    do statistic only, 1:yes; 0:no; default=" << Stat_only_mode << endl;
	cout << "   -h        get help information" << endl;
	exit(0);
}


void get_Ncount_errorRate(string &seq, string &qual, int &N_count, double &error_rate)
{	
	int seq_len = seq.size();
	N_count = 0;
	error_rate = 0.0;
	
	for (int i=0; i<seq_len; i++)
	{	if (seq[i] == 'N')
		{	N_count ++;
			qual[i] = Quality_shift;   // Quality_shift is 64
		}
		error_rate += Qual2Err[(unsigned char)qual[i]];
		
	}
}


//trim的时候,每次trim掉5个bases, 比较reads左边和右边的连续5个bases,哪边错误率高就trim哪边
void trim_low_quality_reads(string &seq, string &qual, double &error_rate)
{	
	
	//trim the qualities of reads from both heads and tails, take 5-bases as one trimming unit
	int seq_len = seq.size();
	int left_pos = 0;
	int right_pos = seq_len - 1;
	while (1)
	{	double left_bases_error = Qual2Err[(unsigned char)qual[left_pos]] + Qual2Err[(unsigned char)qual[left_pos+1]] + Qual2Err[(unsigned char)qual[left_pos+2]]  + Qual2Err[(unsigned char)qual[left_pos+3]] + Qual2Err[(unsigned char)qual[left_pos+4]];
		double right_bases_error = Qual2Err[qual[right_pos]] + Qual2Err[qual[right_pos-1]] + Qual2Err[qual[right_pos-2]] + Qual2Err[qual[right_pos-3]] + Qual2Err[qual[right_pos-4]];
		if (left_bases_error > right_bases_error)
		{	error_rate -= left_bases_error;
			left_pos += 5;
			seq_len -= 5;
		}else
		{	error_rate -= right_bases_error;
			right_pos -= 5;
			seq_len -= 5;
		}
		if (error_rate <= Error_rate_cutoff * seq_len || seq_len <= Min_read_len)
		{	break;
		}
	}

	if (seq_len >= Min_read_len && error_rate <= Error_rate_cutoff*seq_len)
	{	seq = seq.substr(left_pos,seq_len);
		qual = qual.substr(left_pos,seq_len);
		error_rate /= (double)seq_len;
	}else{
		seq = "";
		qual = "";
	}
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "e:m:q:r:t:s:h")) !=-1) {
		switch(c) {
			case 'e': Error_rate_cutoff=atof(optarg); break;
			case 'm': Min_read_len=atoi(optarg); break;
			case 'q': Quality_shift=atoi(optarg); break;
			case 'r': Given_read_length=atoi(optarg); break;
			case 't': Is_trimming=atoi(optarg); break;
			case 's': Stat_only_mode=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 4) usage();
	
	string in_reads1_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string out_reads1_file = argv[optind++];
	string out_stat_file = argv[optind++];

	clock_t time_start, time_end;
	time_start = clock();
	
	time_end = clock();
	cerr << "\nProgram starting\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	uint64_t *Nrate = new uint64_t[Given_read_length+1]; //N number <= read length
	
	//assign initial values to Qual2Err
	for(int i=0; i<41; i++)
	{	Qual2Err[i+Quality_shift] = pow(10.0,-i/10.0);
	}
	int statRankNum = 199; //100 ranks, to make sure the array is big enough to store the data
	StatRank = new uint64_t[statRankNum];
	for (int i=0; i<statRankNum; i++)
	{	StatRank[i] = 0;
	}
	
	igzstream infile1 (in_reads1_file.c_str());
	ogzstream cleanfile1 (out_reads1_file.c_str());
	ofstream statfile (out_stat_file.c_str());
	
	uint64_t total_raw_reads = 0;
	uint64_t total_raw_bases = 0;
	uint64_t total_clean_reads = 0;
	uint64_t total_clean_bases = 0;

	string id1, seq1, qid1, q1;
	while(getline(infile1,id1,'\n'))
	{
		if(id1[0] == '@')
		{
			getline(infile1,seq1,'\n');
			getline(infile1,qid1,'\n');
			getline(infile1,q1,'\n');
			
			int seq_len = seq1.size();
			total_raw_reads += 1;
			total_raw_bases += seq_len;
			
			double err_rate1;
			int N_count1;
			get_Ncount_errorRate(seq1,q1, N_count1, err_rate1);
			
			Nrate[N_count1] ++;
			int rank1 = int(err_rate1/seq_len*100)+1;
			StatRank[rank1] ++;
			
			if (err_rate1 > Error_rate_cutoff * seq_len )
			{	if (Is_trimming)
				{	trim_low_quality_reads(seq1,q1, err_rate1);
					if (seq1.size())
					{	total_clean_reads ++;
						total_clean_bases += seq1.size();
						if (!Stat_only_mode)
						{	cleanfile1 << id1   << "\n" << seq1 << "\n" << qid1 << "\n" <<  q1 << endl;
						}
					}
				}
			}else
			{	
				total_clean_reads ++;
				total_clean_bases += seq1.size();
				if (!Stat_only_mode)
				{	cleanfile1 << id1   << "\n" << seq1 << "\n" << qid1 << "\n" <<  q1 << endl;
				}
			}

		}
	}



	//output the statistic result
	statfile << "#total_raw_reads:   " << total_raw_reads << endl;
	statfile << "#total_raw_bases:   " << total_raw_bases << endl;
	statfile << "#total_clean_reads: " << total_clean_reads << "\t" << total_clean_reads/(double)total_raw_reads*100 << "%" << endl;
	statfile << "#total_clean_bases: " << total_clean_bases << "\t" << total_clean_bases/(double)total_raw_bases*100 << "%" << endl;
	
	uint64_t filtered_lowqual_reads = total_raw_reads - total_clean_reads;
	uint64_t filtered_lowqual_bases = total_raw_bases - total_clean_bases;
	statfile << "#filtered_lowqual_reads: " << filtered_lowqual_reads << "\t" << (double)filtered_lowqual_reads/total_raw_reads*100 << "%"  << endl;
	statfile << "#filtered_lowqual_bases: " << filtered_lowqual_bases << "\t" << (double)filtered_lowqual_bases/total_raw_bases*100 << "%"  << endl;
	
	statfile << "\n#ErrorRank\tNumber_reads\tPercent_of_reads_in_raw_data\n";
	for (int i=1; i<statRankNum; i++)
	{	if (StatRank[i] > 0)
		{	statfile << (double)i << "%\t" << StatRank[i] << "\t" << StatRank[i]/(double)total_raw_reads*100 << "%"  << endl;
		}
	}
	
	statfile << "\n#NsInReads\tNumber_of_reads\tPercent_of_reads_in_raw_data\n";
	for (int i=0; i<Given_read_length+1; i++)
	{	if (Nrate[i] > 0)
		{	statfile << i << "\t" << Nrate[i] << "\t" << Nrate[i]/(double)total_raw_reads*100 << "%" << endl;
		}
	}
	

	time_end = clock();
	cerr << "\nAll jobs finishd done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}


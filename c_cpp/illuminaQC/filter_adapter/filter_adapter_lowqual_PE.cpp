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
string adpt2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";//pDNA-5-
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

//用2-bit表示一个base：A(0), C(1), G(2), T(3), N(4)
int scoreMatrix[5][5] = 
{	{5, -10, -10, -10, -10},
	{-10, 5, -10, -10, -10},
	{-10, -10, 5, -10, -10},
	{-10, -10, -10, 5, -10},
	{-10, -10, -10, -10, 5}
}; //matrix for DNA alignment with 95% identity


void usage() 
{	cout << "filter_raw_reads <reads_1.fq> <reads_2.fq> <reads_1.fq.clean> <reads_2.fq.clean> <result_stat_file>" << endl;
	cout << "   -q <int>    wheter to filter low quality, default=" << Is_filter_quality << endl;
	cout << "   -e <float>  set the error rate cutoff for trimming, default=" << Error_rate_cutoff << endl;
	cout << "   -m <int>    set the minimum trimmed read length, default=" << Min_read_len << endl;
	cout << "   -a <int>    wheter to filter adapter, default=" << Is_filter_adapter << endl;
	cout << "   -1 <str>    give adaper1 sequence for read1, default=" << adpt1 << endl;
	cout << "   -2 <str>    give adaper2 sequence for read12, default=" << adpt2 << endl;
	cout << "   -r <int>    set the sequencing read length, default=" << Given_read_length << endl;
	cout << "   -s <int>    do statistic only, 1:yes; 0:no; default=" << Stat_only_mode << endl;
	cout << "   -h        get help information" << endl;
	exit(0);
}


//外部调用函数
void local_ungapped_aligning (string &seq_i, string &seq_j, int &max_score, int &align_i_start, int &align_j_start, int &align_i_end, int &align_j_end)
{	
//	if (!seq_i.size())
//	{	return;
//	}
	
	int i_size = seq_i.size() + 1;
	int j_size = seq_j.size() + 1;
	max_score = 0;
	
	//i: y-axis; j: x-axis;
	int pre_i;
	int pre_j;
	for (int i=1; i<i_size; i++)
	{	for (int j=1; j<j_size; j++)
		{	pre_i = i - 1;
			pre_j = j - 1;
			int this_score = DPscore[pre_i * j_size + pre_j] + scoreMatrix[ alphabet[ seq_i[pre_i] ] ][ alphabet[ seq_j[pre_j] ] ];
			if (this_score < 0)
			{	this_score = 0;
			}
			DPscore[i * j_size + j] = this_score; //(i+1), (j+1)代表1-based order
			if (this_score > max_score)
			{	max_score = this_score;
				align_i_end = i;
				align_j_end = j;
			}
		}
	}
	
	//trace back by looping instead of recursion
	if (max_score >= Adapter_align_cutoff)
	{
		int pos_i = align_i_end;
		int pos_j = align_j_end;
		int score = max_score;
		while (score > 0)
		{	
			//align_i.push_back( seq_i[pos_i - 1] );
			//align_j.push_back( seq_j[pos_j - 1] );
			pos_i --;
			pos_j --;
			score = DPscore[pos_i * j_size + pos_j];
		}
		align_i_start = pos_i + 1;
		align_j_start = pos_j + 1;
		//reverse( align_i.begin(), align_i.end() );
		//reverse( align_j.begin(), align_j.end() );
		//cerr << "after trace back\n";
	}

}

//对含有N的reads加强处理：超过read长度一半是N的直接扔掉
//trim的时候,每次trim掉5个bases
//本处理逻辑仍有漏洞，可能多删reads，需要进一步完善
void filter_trim_low_quality(string &seq, string &qual, double &error_rate)
{	
	int seq_len = seq.size();
	if (seq_len == 0)
	{	return;
	}
	int N_count = 0;
	error_rate = 0.0;
	
	for (int i=0; i<seq_len; i++)
	{	if (seq[i] == 'N')
		{	N_count ++;
			qual[i] = '@';  //质量值设为0
		}
		error_rate += Qual2Err[qual[i]];
	}
	Nrate[N_count] ++;

	//deal with reads with more than one N bases
	if (seq_len - N_count < Min_read_len)
	{	seq_len = 0;
		seq = "";
		qual = "";
		return;
	}
	
	//trim the qualities of reads from both heads and tails, take 5-bases as one trimming unit
	if (error_rate > Error_rate_cutoff * seq_len )
	{	
		int left_pos = 0;
		int right_pos = seq_len - 1;
		while (1)
		{	double left_bases_error = Qual2Err[qual[left_pos]] + Qual2Err[qual[left_pos+1]] + Qual2Err[qual[left_pos+2]]  + Qual2Err[qual[left_pos+3]] + Qual2Err[qual[left_pos+4]];
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
		}else{
			seq = "";
			qual = "";
		}
	}

	if (seq_len)
	{	error_rate /= (double)seq_len;
	}
	
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "q:e:m:a:1:2:r:s:h")) !=-1) {
		switch(c) {
			case 'q': Is_filter_quality=atoi(optarg); break;
			case 'e': Error_rate_cutoff=atof(optarg); break;
			case 'm': Min_read_len=atoi(optarg); break;
			case 'a': Is_filter_adapter=atoi(optarg); break;
			case '1': adpt1=optarg; break;
			case '2': adpt2=optarg; break;
			case 'r': Given_read_length=atoi(optarg); break;
			case 's': Stat_only_mode=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 6) usage();
	
	string in_reads1_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string in_reads2_file = argv[optind++];
	string out_reads1_file = argv[optind++];
	string out_reads2_file = argv[optind++];
	string out_stat_file = argv[optind++];

	clock_t time_start, time_end;
	time_start = clock();
	
	time_end = clock();
	cerr << "\nProgram starting\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	cerr << "Is_filter_quality: " << Is_filter_quality << endl;
	cerr << "Error_rate_cutoff: " << Error_rate_cutoff << endl;
	cerr << "Min_trimmed_read_len:      " << Min_read_len << endl;
	cerr << "Is_filter_adapter: " << Is_filter_adapter << endl;
	cerr << "Given_read_length: " << Given_read_length << endl;
	cerr << "adapter1_sequence: " << adpt1 << endl;
	cerr << "adapter2_sequence: " << adpt2 << endl;


	//assign initial values to Qual2Err
	for(int i=0; i<41; i++)
	{	Qual2Err[i+Qual_shift] = pow(10.0,-i/10.0);
	}
	int statRankNum = 1001;
	StatRank = new uint64_t[statRankNum];
	for (int i=0; i<statRankNum; i++)
	{	StatRank[i] = 0;
	}
	
	//assign initial values for adapter alignment
	int adapter_size = (adpt1.size() >= adpt2.size()) ? adpt1.size() : adpt2.size();
	int array_size = (Given_read_length + 1) * (adapter_size + 1);
	DPscore = new int[array_size];
	//construct DP matrix
	DPscore[0] = 0; //stands for score
	for (int j=1; j<=adapter_size; j++)
	{	DPscore[j] = 0;  
	}
	
	for (int i=1; i<=Given_read_length; i++)
	{	DPscore[i * (adapter_size+1)] = 0;
	}
	
	igzstream infile1 (in_reads1_file.c_str());
	igzstream infile2 (in_reads2_file.c_str());
	ogzstream cleanfile1 (out_reads1_file.c_str());
	ogzstream cleanfile2 (out_reads2_file.c_str());
	ofstream statfile (out_stat_file.c_str());
	
	uint64_t total_raw_reads = 0;
	uint64_t total_raw_bases = 0;
	uint64_t total_clean_reads = 0;
	uint64_t total_clean_bases = 0;
	uint64_t total_clean_pairs = 0;
	uint64_t filtered_reads_lowqual = 0;
	uint64_t stat_reads1_lowqual = 0;
	uint64_t stat_reads2_lowqual = 0;
	uint64_t filtered_adapter_reads = 0;
	uint64_t stat_read1_adapter = 0;
	uint64_t stat_read2_adapter = 0;

	string id1, seq1, qid1, q1, id2, seq2, qid2, q2;
	while(getline(infile1,id1,'\n'))
	{
		if(id1[0] == '@')
		{
			getline(infile1,seq1,'\n');
			getline(infile1,qid1,'\n');
			getline(infile1,q1,'\n');
			
			getline(infile2,id2,'\n');
			getline(infile2,seq2,'\n');
			getline(infile2,qid2,'\n');
			getline(infile2,q2,'\n');
			
			total_raw_reads += 2;
			total_raw_bases += seq1.size() + seq2.size();
			int is_contain_adapter = 0;

			if (Is_filter_adapter)
			{	
				int align_score1, seq1_start, seq1_end, adpt1_start, adpt1_end;
				int align_score2, seq2_start, seq2_end, adpt2_start, adpt2_end;
				local_ungapped_aligning (seq1, adpt1, align_score1, seq1_start, adpt1_start, seq1_end, adpt1_end);
				local_ungapped_aligning (seq2, adpt2, align_score2, seq2_start, adpt2_start, seq2_end, adpt2_end);
				
				if (align_score1 >= Adapter_align_cutoff || align_score2 >= Adapter_align_cutoff)
				{	if (align_score1 >= Adapter_align_cutoff)
					{	stat_read1_adapter ++;
					}
					if (align_score2 >= Adapter_align_cutoff)
					{	stat_read2_adapter ++;
					}
					
					if (align_score1 >= Adapter_align_cutoff && align_score2 >= Adapter_align_cutoff)
					{	is_contain_adapter = 1;
						filtered_adapter_reads += 2;
						seq1 = "";
						seq2 = "";
						int insert1_size = seq1_start - adpt1_start;
						int insert2_size = seq2_start - adpt2_start;
						if ( insert1_size == insert2_size && insert1_size >= 0)
						{	AdptInsert[insert1_size] ++;
						}
						
					}
				}
			}
			
			if (Is_filter_quality && !is_contain_adapter) 
			{
				double err_rate1, err_rate2;
				filter_trim_low_quality(seq1,q1,err_rate1);
				filter_trim_low_quality(seq2,q2,err_rate2);

				if (seq1.size())
				{	int rank1 = int(err_rate1*1000)+1;
					StatRank[rank1] ++;
				}else
				{	stat_reads1_lowqual ++;
				}
				if (seq2.size())
				{	int rank2 = int(err_rate2*1000)+1;
					StatRank[rank2] ++;
				}else
				{	stat_reads2_lowqual ++;
				}
			}
			
			//output the filtered clean result
			if (seq1.size() && seq2.size())
			{	
				if (seq1.size())
				{	total_clean_reads ++;
					total_clean_bases += seq1.size();
				}
				if (seq2.size())
				{	total_clean_reads ++;
					total_clean_bases += seq2.size();
				}
				if (seq1.size() && seq2.size())
				{	total_clean_pairs ++;
				}
				
				if (!Stat_only_mode)
				{	cleanfile1 << id1   << "\n" << seq1 << "\n" << qid1 << "\n" <<  q1 << endl;
					cleanfile2 << id2   << "\n" << seq2 << "\n" << qid2 << "\n" <<  q2 << endl;
				}
			}

		}
	}
	delete []DPscore;
	
	//output the statistic result
	statfile << "#total_raw_reads:   " << total_raw_reads << endl;
	statfile << "#total_raw_bases:   " << total_raw_bases << endl;
	statfile << "#total_clean_reads: " << total_clean_reads << "\t" << total_clean_reads/(double)total_raw_reads << endl;
	statfile << "#total_clean_bases: " << total_clean_bases << "\t" << total_clean_bases/(double)total_raw_bases << endl;
	statfile << "#total_clean_pairs: " << total_clean_pairs << "\t" << total_clean_pairs/((double)total_raw_reads/2) << endl;
	
	if (Is_filter_adapter)
	{	
		statfile << "\n#stat_read1_adapter: " << stat_read1_adapter << "\t" << (double)stat_read1_adapter/(total_raw_reads/2) << endl;
		statfile << "#stat_read2_adapter: " << stat_read2_adapter << "\t" << (double)stat_read2_adapter/(total_raw_reads/2) << endl;
		statfile << "#filtered_adapter_reads: " << filtered_adapter_reads << "\t" << (double)filtered_adapter_reads/total_raw_reads << endl;
		uint64_t pairs_with_insert = 0;
		for (int i=0; i<250; i++)
		{	if (AdptInsert[i] > 0)
			{	pairs_with_insert += AdptInsert[i];
			}
		}
		statfile << "#reads_with_insert_size: " << pairs_with_insert*2 << "\t" << (double)pairs_with_insert*2/filtered_adapter_reads << endl;
		statfile << "#InsertSize\tNumber_read_pairs\tRatio_of_reads_with_insert_size\n";
		for (int i=0; i<250; i++)
		{	if (AdptInsert[i] > 0)
			{	statfile << i << "\t" << AdptInsert[i] << "\t" << AdptInsert[i]/(double)pairs_with_insert << endl;
			}
		}
	}

	if (Is_filter_quality)
	{	uint64_t filtered_reads_lowqual = stat_reads1_lowqual + stat_reads2_lowqual;
		statfile << "\n#filtered_reads1_lowqual: " << stat_reads1_lowqual << "\t" << (double)stat_reads1_lowqual/(total_raw_reads/2) << endl;
		statfile << "#filtered_reads2_lowqual: " << stat_reads2_lowqual << "\t" << (double)stat_reads2_lowqual/(total_raw_reads/2)  << endl;
		statfile << "#filtered_reads_lowqual: " << filtered_reads_lowqual << "\t" << (double)filtered_reads_lowqual/total_raw_reads << endl;
		
		statfile << "#QualRank\tNumber_reads\tRatio_of_reads_in_clean_data\n";
		for (int i=1; i<statRankNum; i++)
		{	if (StatRank[i] > 0)
			{	statfile << (double)i/10 << "%\t" << StatRank[i] << "\t" << StatRank[i]/(double)total_clean_reads << endl;
			}
		}
		
		statfile << "\n#NsInReads\tNumber_of_reads\tRatio_of_reads_in_raw_data\n";
		for (int i=0; i<250; i++)
		{	if (Nrate[i] > 0)
			{	statfile << i << "\t" << Nrate[i] << "\t" << Nrate[i]/(double)total_raw_reads << endl;
			}
		}
	}

	time_end = clock();
	cerr << "\nAll jobs finishd done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}


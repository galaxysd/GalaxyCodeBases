//filter illumina reads that contained adapters, and make statistics
//the algorithm used is local ungapped alignment, it is a little slow but most sensitive.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <zlib.h>
#include "gzstream.h"

using namespace std;

string adapter_type = "gDNA-3"; // default adapter type
string adpt_3end="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";//pDNA-3+ for read_1
string adpt_5end="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";//pDNA-5- for read_2
int Given_read_length = 100;
int Adapter_align_cutoff = 50;  //minimum alignment score
int Stat_only_mode = 0;
int Is_trimming = 0;
int Min_read_len = 50;
int *DPscore;

unsigned char alphabet[128] =
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
{	cout << "filter_adapter <reads_x.fq>  <reads_x.fq.clean>  <result_stat_file>" << endl;
	cout << "   -a <str>    set adaper type gDNA-3 or gDNA-5, default=" << adapter_type << endl;
	cout << "   -r <int>    set the sequencing read length, default=" << Given_read_length << endl;
	cout << "   -c <int>    set the minimum alignment score, default=" << Adapter_align_cutoff << endl;
	cout << "   -m <int>    set the minimum trimmed read length, default=" << Min_read_len << endl;
	cout << "   -t <int>    whether trim, 0:no; 1:yes; default=" << Is_trimming << endl;
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


int main(int argc, char *argv[])
{
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "a:r:c:m:t:s:h")) !=-1) {
		switch(c) {
			case 'a': adapter_type=optarg; break;
			case 'r': Given_read_length=atoi(optarg); break;
			case 'c': Adapter_align_cutoff=atoi(optarg); break;
			case 'm': Min_read_len=atoi(optarg); break;
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

	string adapter;
	if (adapter_type == "gDNA-3")
	{	adapter = adpt_3end;
	}else if (adapter_type == "gDNA-5")
	{	adapter = adpt_5end;
	}

	clock_t time_start, time_end;
	time_start = clock();

	time_end = clock();
	cerr << "\nProgram starting\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	//assign initial values for adapter alignment
	int adapter_size = adapter.size();
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
	ogzstream cleanfile1 (out_reads1_file.c_str());
	ofstream statfile (out_stat_file.c_str());

	uint64_t total_raw_reads = 0;
	uint64_t total_raw_bases = 0;
	uint64_t total_clean_reads = 0;
	uint64_t total_clean_bases = 0;
	uint64_t *AdptInsert = new uint64_t[Given_read_length]; //insert <= read_length
	uint64_t total_reads_with_adapter = 0;

	string id1, seq1, qid1, q1;
	while(getline(infile1,id1,'\n'))
	{
		if(id1[0] == '@')
		{
			getline(infile1,seq1,'\n');
			getline(infile1,qid1,'\n');
			getline(infile1,q1,'\n');

			total_raw_reads += 1;
			total_raw_bases += seq1.size();

			int align_score1, seq1_start, seq1_end, adapter_start, adapter_end;
			local_ungapped_aligning (seq1, adapter, align_score1, seq1_start, adapter_start, seq1_end, adapter_end);

			if (align_score1 >= Adapter_align_cutoff)
			{
				int insert1_size = seq1_start - adapter_start;
				if ( insert1_size < 0) { insert1_size = 0; }
				AdptInsert[insert1_size] ++;
				total_reads_with_adapter ++;
				if (Is_trimming && insert1_size >= Min_read_len)
				{	total_clean_reads ++;
					total_clean_bases += insert1_size;
					if (!Stat_only_mode)
					{	cleanfile1 << id1   << "\n" << seq1.substr(0,insert1_size) << "\n" << qid1 << "\n" <<  q1 << endl;
					}
				}

			}else
			{	total_clean_reads ++;
				total_clean_bases += seq1.size();
				if (!Stat_only_mode)
				{	cleanfile1 << id1   << "\n" << seq1 << "\n" << qid1 << "\n" <<  q1 << endl;
				}
			}
		}
	}
	delete []DPscore;

	//output the statistic result
	statfile << "#total_raw_reads:   " << total_raw_reads << endl;
	statfile << "#total_raw_bases:   " << total_raw_bases << endl;
	statfile << "#total_clean_reads: " << total_clean_reads << "\t" << total_clean_reads/(double)total_raw_reads*100 << "%" << endl;
	statfile << "#total_clean_bases: " << total_clean_bases << "\t" << total_clean_bases/(double)total_raw_bases*100 << "%" << endl;
	statfile << "\n#reads_with_adapter: " << total_reads_with_adapter << "\t" << (double)total_reads_with_adapter/total_raw_reads*100 << "%" << endl;

	statfile << "\n#InsertSize\tNumber_reads\tPercent_in_reads_with_adapters\n";
	for (int i=0; i<Given_read_length; i++)
	{	if (AdptInsert[i] > 0)
		{	statfile << i << "\t" << AdptInsert[i] << "\t" << AdptInsert[i]/(double)total_reads_with_adapter*100 << "%" << endl;
		}
	}

	time_end = clock();
	cerr << "\nAll jobs finishd done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}


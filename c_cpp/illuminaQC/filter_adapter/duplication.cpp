#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <time.h>
#include <cmath>
#include "gzstream.h"
//author: shiyujian

using namespace std;
void usage ()
{
        cout<<"\nUsage: stat_dup [options] <read1.fq> <read2.fq>\n"
        		<<"  -n <int>   loading reads number\n"
        		<<"  -e <float> the max error rate cutoff\n"
        		<<"  -s <str>   the output stat file\n"
        		<<"  -a <int>   delete from the initiation of read1\n"
        		<<"  -b <int>   delete from the end of read1\n"
        		<<"  -c <int>   delete from the initiation of read2\n"
        		<<"  -d <int>   delete from the end of read2\n"
        		<<"  -h         outpput help information\n" 
        		<<endl;
        exit(1);
}

typedef unsigned long long bit64_t; // 64-bit
typedef unsigned int bit8_t; //8-bit
typedef struct read1
{
	bit64_t AA[7];
}Read1;

struct Compare
{
	bool operator()(const read1 &a, const read1 &b)
	{
		for (int i=0;i<7 ;i++ )
		{
			if (a.AA[i] < b.AA[i])
			{
				return true;
			}
			else if (a.AA[i] > b.AA[i])
			{
				return false;
			}
		}
		return false;
	}
};


Read1 in_read1(string &seq,int n)//n<7
{
	Read1 base;
	base.AA[0]=0;
	base.AA[1]=0;
	base.AA[2]=0;
	base.AA[3]=0;
	base.AA[4]=0;
	base.AA[5]=0;
	base.AA[6]=0;
	int x=0;
	for (int m=0;m<n ;m++ )
	{
		x=m+1;
		for (int i=32*m;i<32*x ;i++ )
		{
			int k=((seq[i]&0x06)>>1);
			base.AA[m]=base.AA[m]<<2;
			base.AA[m]+=k;
		}
	}
	for (int i=32*n;i<seq.size() ;i++ )
	{
		int k=((seq[i]&0x06)>>1);
		base.AA[n]=base.AA[n]<<2;
		base.AA[n]+=k;
	}
	return base;
}

int main (int argc, char *argv[])
{
	//deal with command line options
	int c; //must be int type
	int total_reads_num,seq1_head,seq1_tail,seq2_head,seq2_tail,Q_shift=64;
	double max_err_rate;
	string stat_file;
	while((c=getopt(argc,argv,"n:a:b:c:d:e:s:h")) != -1) { // :表示有参数值
		switch ( c )
		{	// return optarg, type is char*, use atoi or atof to convert
			case 'n' : total_reads_num = atoi(optarg); break;
			case 'a' : seq1_head = atoi(optarg); break;
			case 'b' : seq1_tail = atoi(optarg); break;
			case 'c' : seq2_head = atoi(optarg); break;
			case 'd' : seq2_tail = atoi(optarg); break;
			case 'e' : max_err_rate = atof(optarg); break;
			case 'q' : Q_shift = atoi(optarg); break;
			case 's' : stat_file = optarg; break;
			case 'h' : usage(); break;
			default  : usage(); 
		}
	}
	if (argc<2) usage();
	string in_seq1=argv[optind++];
	string in_seq2=argv[optind++];

	string id1="", seq1="", qid1="", q1="", id2="", seq2="", qid2="", q2="";
	int read1_len=0, read2_len=0, total_len=0;

	igzstream infile1( in_seq1.c_str());
	if ( ! infile1)
	{
		cerr<<"fail to open input file"<<in_seq1<<endl;
	}
	getline ( infile1, id1, '\n');
	getline ( infile1, seq1, '\n');
	read1_len=seq1.length();

	igzstream infile2( in_seq2.c_str());
	if ( ! infile2)
	{
		cerr<<"fail to open input file"<<in_seq2<<endl;
	}
	getline ( infile2, id2, '\n');
	getline ( infile2, seq2, '\n');
	read2_len=seq2.length();
	
	int read1_use_len = read1_len - seq1_head - seq1_tail;
	int read2_use_len = read2_len - seq2_head - seq2_tail;
	total_len=read1_use_len + read2_use_len;
	if (total_len>224)
	{
		cerr<<"max pair reads length is 224"<<endl;
		exit(1);
	}
	if (total_len<1)
	{
		cerr<<"format error"<<endl;
		exit(1);
	}
	double q2e[2][100];
	for(int i=0;i<100;i++)
	{
		q2e[0][i] = pow(10.0,double(i)/(-10.0));
		q2e[1][i] = 1;
	}
	double max_err_num = max_err_rate * (read1_use_len + read2_use_len);
	float size_len=total_len/32.0;
	cout<<"size_len"<<size_len<<endl;
	int size_num=(int)size_len;
	if (size_num==size_len)
	{
		size_num--;
	}
	cout<<"size_num:"<<size_num<<endl;
	infile1.close();
	infile2.close();
	id1="", id2="", seq1="", seq2="";

//	infile1.open ( in_seq1.c_str());
//	infile2.open ( in_seq2.c_str());

	ofstream stat (stat_file.c_str());
	if (! stat)
	{
		cerr << "fail to create stat file" <<stat_file<< endl;
	}
	
	igzstream infile11( in_seq1.c_str());
	if ( ! infile1)
	{
		cerr<<"fail to open input file"<<in_seq1<<endl;
	}

	igzstream infile22( in_seq2.c_str());
	if ( ! infile2)
	{
		cerr<<"fail to open input file"<<in_seq2<<endl;
	}
	map<Read1, int, Compare> read_count;
	int clean=0,total_num=0,dup=0;
	time_t start_time,end_time;
	time(&start_time);
	while ( getline ( infile11, id1, '\n') )
	{
		if (id1[0] == '@')
		{
			getline ( infile11, seq1, '\n');
			getline ( infile11, qid1, '\n');
			getline ( infile11, q1, '\n');
			getline ( infile22, id2, '\n');
			getline ( infile22, seq2, '\n');
			getline ( infile22, qid2, '\n');
			getline ( infile22, q2, '\n');
			string read=seq1.substr(seq1_head,read1_use_len) + seq2.substr(seq2_head,read2_use_len);
			string qual=q1.substr(seq1_head,read1_use_len) + q2.substr(seq2_head,read2_use_len);
			double err_num = 0;
			for(int i=0;i<read.length();i++)
				err_num += q2e[!(read[i]-'N')][qual[i]-Q_shift];
			if(err_num > max_err_num)
				continue;
				
			total_num++;
			Read1 read_seq=in_read1(read,size_num);
			//cout<<read_count[read_seq]<<endl;
			if (read_count[read_seq]>0)
			{
				read_count[read_seq]++;
				//cout<<id1<<endl<<seq1<<endl<<qid1<<endl<<q1<<endl;
				//cout<<id2<<endl<<seq2<<endl<<qid2<<endl<<q2<<endl;
			}
			else
			{
				read_count[read_seq]=1;
//				clean++;
//				outfile1<<id1<<endl<<seq1<<endl<<qid1<<endl<<q1<<endl;
//				outfile2<<id2<<endl<<seq2<<endl<<qid2<<endl<<q2<<endl;
			}
			if(total_num >= total_reads_num)
				break;
		}
	}
	
	int freq[256]={0},max_index=0;
	map<Read1, int>::const_iterator map_it = read_count.begin();
	int max_freq_num=0;
	while (map_it != read_count.end())
	{
		if(map_it->second > max_index)
			max_index = map_it->second;
		if(map_it->second>255)
		{
			max_freq_num += map_it->second;
			freq[255]++;
			max_index = 255;
		}
		else
			freq[map_it->second]++;
		if (map_it->second >1)
		{
			dup+=map_it->second;
		}
		//cout<<map_it->second<<endl;
		map_it++;
	}
//	for(int i=2;i<=max_index;i++)
//		dup += freq[i];
	time(&end_time);
	stat<<"#total_reads_num: "<<total_num<<endl
		<<"#Duplicate_reads: "<<dup<<endl
		<<"#Duplicate_rate: "<<double(dup)/double(total_num)<<endl
		<<"#Used_time: "<<end_time-start_time<<endl
		<< "#freq\tSpecies_Number\tIndividual_Number\tratio" << endl;
	for(int i=1;i<max_index;i++)
		stat<<i<<'\t'<<freq[i]<<'\t'<<i*freq[i]<<'\t'<<double(i*freq[i])/double(total_num)<<endl;
	if(max_index<255)
		stat<<max_index<<'\t'<<freq[max_index]<<'\t'<<max_index*freq[max_index]<<'\t'<<double(max_index*freq[max_index])/double(total_num)<<endl;
	else
		stat<<">=255\t"<<freq[255]<<'\t'<<max_freq_num<<'\t'<<double(max_freq_num)/double(total_num)<<endl;
//	stat<<"Total_reads:"<<total_num<<endl<<"Duplicate_reads:"<<dup<<endl<<"Clean_reads:"<<clean<<endl<<"Used_time:"<<end_time-start_time<<endl;
	infile11.close();
	infile22.close();
	stat.close();
}


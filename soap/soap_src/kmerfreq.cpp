#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<unistd.h>


using namespace std;

typedef unsigned char bit8_t;
typedef unsigned long bit64_t;
typedef unsigned int bit32_t;

int Ele=sizeof(bit64_t)*4;
const bit8_t nv=0;  //convert unknown char as 'A'
bit8_t alphabet[256] =
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */ 
 0,nv, 1,nv,nv,nv, 2,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'P' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'a' */ 
 0,nv, 1,nv,nv,nv, 2,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'p' */
nv,nv,nv,nv, 3,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};
char sbin[4] =
{
	'A', 'C', 'G', 'T'
};
struct Inf
{
	bit32_t f;
	bit64_t mer;
};

bool more_abund(Inf a, Inf b)
{
	return a.f>=b.f;
};

int kmer=15;
int initial=2;
float min_abund=0;
int min_num=0;
char *infile=NULL;
bit64_t Mer;

bool ConvertKmer(string &s)
{
	Mer=0;
	for(bit32_t i=initial-1; i<initial-1+kmer; i++) {
		if(alphabet[s[i]]==4)
			return 0;
		Mer=(Mer<<2)|alphabet[s[i]];
	}
	return 1;
}

string ConvertSeq(bit64_t e)
{
	string s;
	for(int i=0; i<kmer; i++) {
		s.push_back(sbin[(e>>(kmer-1-i)*2)&0x3]);
	}
	return s;
}

void usage(void)
{
	cout<<"Usage:	kmerfreq [options] <reads.fa|fq> <stdout>\n"
		<<"	-k	<int>	kmer size, default="<<kmer<<"\n"
		<<"	-i	<int>	initial bp, default="<<initial<<"\n"
		<<"	-b	<float>	minimal frequency threshold, default="<<min_abund<<"\n"
		<<"	-m	<int>	minimal appearance threshold, default="<<min_num<<"\n"
		<<"	-h	help\n\n"
		<<"Output format: kmer,number,% of total.\n\n"
		<<"Note: it will need 4^15*4=4Gb RAM to build hash for 15-mer, then 4GB RAM to store refseq.\n\n";
}
int main(int argc, char *argv[])
{
	if(argc<2)
		usage();
	int c;
	while((c=getopt(argc, argv, "k:i:b:m:h")) !=-1) {
		switch(c) {
			case 'k': kmer=atoi(optarg); break;
			case 'i': initial=atoi(optarg); break;
			case 'b': min_abund=atof(optarg); break;
			case 'm': min_num=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if(optind<argc)
		infile=argv[optind];
	ifstream fin(infile);
	if(!fin) {
		cerr<<"fatal error: failed to open file "<<infile<<endl;
		exit(1);
	}
	bit32_t *freq;
	bit64_t total=0;
	bit64_t num_kmers=0;
	for(int i=0; i<kmer; i++) {
		total=(total<<2)|0x3;
	}
	freq=new bit32_t[total];
	for(bit64_t i=0; i<total; i++) {
		freq[i]=0;
	}
	string id, seq, qual;
	char ch[1000];
	while(!fin.eof()) {
		fin>>id;
		if(fin.eof())
			break;
		fin.getline(ch, 1000);
		fin>>seq;
		if(id[0] == '@') {
			fin.getline(ch, 1000);
			fin.getline(ch, 1000);
			fin>>qual;
			fin.getline(ch, 1000);
		}
		if(ConvertKmer(seq)) {
			freq[Mer]++;
			num_kmers++;
		}
	}
	fin.close();
	//
	if(min_num==0)
		min_num=int(num_kmers*min_abund);
	vector<Inf> a;
	Inf e;
	for(bit64_t i=0; i<total; i++) {
		if(freq[i] >= min_num) {
			e.f=freq[i];
			e.mer=i;
			a.push_back(e);
		}
	}
	sort(a.begin(), a.end(), more_abund);
	for(vector<Inf>::iterator p=a.begin(); p!=a.end(); p++) {
		cout<<ConvertSeq(p->mer)<<"\t"<<p->f<<"\t"<<float(100*p->f)/num_kmers<<"%"<<"\n";
	}
	return 0;
}

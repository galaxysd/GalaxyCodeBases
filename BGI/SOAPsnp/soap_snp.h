#ifndef SOAP_SNP_HH_
#define SOAP_SNP_HH_
#include <iostream>
#include "gzstream.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <stdlib.h> 
#include <semaphore.h>
#include "tools.h"
#include "SamCtrl.h"
#include "SfsMethod.h"
#include "ThreadManager.h"

#define OPENSFS_ERROR -1;
#define OPENSFS_SUCC  1;
#define POINTER_NULL  -1;

typedef unsigned long long ubit64_t;
typedef unsigned int ubit32_t;
typedef double rate_t;
typedef unsigned char small_int;
using namespace std;
typedef ifstream my_ifstream;
typedef ofstream my_ofstream;
typedef myfstream gzoutstream;  //update by zhukai on 2010-12-09

const size_t capacity = sizeof(ubit64_t)*8/4;
const char abbv[17]={'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};
const ubit64_t glf_base_code[8]={1,2,8,4,15,15,15,15}; // A C T G
const ubit64_t glf_type_code[10]={0,5,15,10,1,3,2,7,6,11};// AA,CC,GG,TT,AC,AG,AT,CG,CT,GT
const int global_win_size = 1000;

// Some global variables
class Files 
{
public:
	igzstream soap_result;
	my_ifstream ref_seq, dbsnp, region;
	my_ofstream /*consensus,*/ baseinfo, o_region;
	//ogzstream consensus;//update by zhukai 2010-11-30
	gzoutstream * consensus;
	fstream matrix_file;
	SamCtrl sam_result;

	FILE *sfsfile; // sfs file handle 

	FILE *jointSfsfile; // joint file handle
	FILE *freqfile; // freq file handle
	
/*	ofstream sfsfile;
	ofstream jointSfsfile;
	ofstream freqfile;
*/
	string fSFSall;
	string fFreq;
	string fJoint;

	Files()
	{
		soap_result.close();
		ref_seq.close();
		dbsnp.close();
		//consensus.clear();
		baseinfo.close();
		matrix_file.close();
		region.close();
		o_region.close();
		sfsfile = NULL;
		jointSfsfile = NULL;
		freqfile = NULL;
		consensus = NULL;
	};
	
public:
	virtual ~Files()
	{
		if (freqfile != NULL)
			fclose(freqfile);
		if (sfsfile != NULL)
			fclose(sfsfile);
		if (jointSfsfile != NULL)
			fclose(jointSfsfile);
		/*if (consensus != NULL)
			delete consensus;*/
	};
	/*open sfsfile for output*/
	int OpenSfsfile(const string outfiles, const int writeFr, const int doBay, const int doJoint);
	/*get file in mode */
	FILE* getFILE(const char*fname,const char* mode);
};

class Parameter 
{
public:
	//guyue
	std::string alignment_name, consensus_name,soapresult_name;
	std::string matrix_file;
	std::string sfs_path;
	bool is_matrix_in ; // Generate the matrix or just read it?
	bool file_list ;
	int c;
	int ret;
	char in_mode[5];
	int CPU;

	sem_t sem_read;
	sem_t sem_call_cns;
	sem_t sem_readwin_return;
	
	char q_min; // The char stands for 0 in fastq
	char q_max; // max quality score
	int read_length; // max read length
	bool is_monoploid; // Is it an monoploid? chrX,Y,M in man.
	bool is_snp_only;  // Only output possible SNP sites?
	bool refine_mode; // Refine prior probability using dbSNP
	bool rank_sum_mode; // Use rank sum test to refine HET quality
	bool binom_mode; // Use binomial test to refine HET quality
	bool transition_dominant; // Consider transition/transversion ratio?
	int glf_format; // Generate Output in GLF format: New File Definitions! Since May 1, 2009
	bool region_only; // Only report consensus in specified region
	std::string glf_header; // Header of GLF format
	rate_t althom_novel_r, het_novel_r; // Expected novel prior
	rate_t althom_val_r, het_val_r; // Expected Validated dbSNP prior
	rate_t althom_unval_r, het_unval_r; // Expected Unvalidated dbSNP prior
	rate_t global_dependency, pcr_dependency; // Error dependencies, 1 is NO dependency

	SFS_PARA *sfs; // the parameters for sfs.
	bool is_cns_gz;
	
// Default onstruction
	Parameter()
	{
		//guyue
		alignment_name = "" ;
		consensus_name = "" ;
		soapresult_name = "" ;
		matrix_file = "" ;
		sfs_path = "" ;
		is_matrix_in = false; // Generate the matrix or just read it?
		file_list = false;
		ret = 0;
		c = 0;
		in_mode[0] = 0;
		in_mode[1] = 0;
		in_mode[2] = 0;
		in_mode[3] = 0;
		in_mode[4] = 0;
		CPU = 4;
		sem_init(&sem_read, 0, 2);
		sem_init(&sem_call_cns, 0, 0);
		sem_init(&sem_readwin_return, 0, 1);

		q_min = 64;
		q_max = 64+40;
		read_length = 45;
		is_monoploid = is_snp_only = refine_mode = rank_sum_mode = binom_mode = transition_dominant = region_only =false;
		glf_format = 0;
		glf_header = "";
		althom_novel_r=0.0005, het_novel_r=0.0010;
		althom_val_r=0.05, het_val_r=0.10;
		althom_unval_r=0.01, het_unval_r=0.02;
		global_dependency= log10(0.9), pcr_dependency= log10(0.5); // In Log10 Scale
		sfs = NULL;
		is_cns_gz = false;
	};
	~Parameter()
	{
		/*if (sfs != NULL)
			sfs_destroy(sfs);*/
		sem_destroy(&sem_read);
		sem_destroy(&sem_call_cns);
	};
};

class Soap_format 
{
	// Soap alignment result
	std::string read_id, read, qual, chr_name;
	int hit, read_len, mismatch;
	int position;
	char ab, strand;
public:
	Soap_format()
	{
		hit = 0;
		position = 0;
		read_len = 0;
		mismatch = 0;
		chr_name = "";
	}

	// changed by Bill 2010-10-11
	friend istringstream & operator>>(std::istringstream & alignment, Soap_format & soap)
	{
		alignment>>soap.read_id>>soap.read>>soap.qual>>soap.hit>>soap.ab>>soap.read_len>>soap.strand>>soap.chr_name>>soap.position>>soap.mismatch;
		//cerr<<soap<<endl;
		//exit(1);
		//if (soap.mismatch>200) 
		//{
		//	int indel_pos,indel_len;
		//	alignment>>indel_pos;
		//	// add by Bill 2010-11-22
		//	int clip_num = 0, last_clip_num;
		//	string cigar;
		//	alignment >> cigar;
		//	clip_num = count_soft_clip(cigar, last_clip_num);
		//	indel_pos -= clip_num;

		//	indel_len = soap.mismatch-200;
		//	soap.read = soap.read.substr(0,indel_pos) + soap.read.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
		//	soap.qual = soap.qual.substr(0,indel_pos) + soap.qual.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
		//}
		//else if(soap.mismatch>100) 
		//{
		//	int indel_pos,indel_len;
		//	string temp("");
		//	alignment>>indel_pos;
		//	// add by Bill 2010-11-22
		//	int clip_num = 0, last_clip_num;
		//	string cigar;
		//	alignment >> cigar;
		//	clip_num = count_soft_clip(cigar, last_clip_num);
		//	indel_pos -= clip_num;

		//	indel_len = indel_len - 100;
		//	for(int i = 0; i != indel_len; i++) 
		//	{
		//		temp = temp+'N';
		//	}
		//	soap.read = soap.read.substr(0,indel_pos)+temp+soap.read.substr(indel_pos,soap.read_len-indel_pos);
		//	soap.qual = soap.qual.substr(0,indel_pos)+temp+soap.qual.substr(indel_pos,soap.read_len-indel_pos);
		//}
		if(soap.mismatch>200) 
		{ //deletion
			int indel_pos,indel_len;
			string temp("");
			alignment>>indel_pos;
			indel_len = soap.mismatch-200;
			for(int i=0; i!=indel_len; i++) 
			{
				temp = temp+'N';
			}
			soap.read = soap.read.substr(0,indel_pos) + temp + soap.read.substr(indel_pos,soap.read_len-indel_pos);
			soap.qual = soap.qual.substr(0,indel_pos) + temp + soap.qual.substr(indel_pos,soap.read_len-indel_pos);
			//cerr<<soap<<endl;
		}
		else if (soap.mismatch>100) 
		{ //insertion
			int indel_pos,indel_len;
			alignment>>indel_pos;
			indel_len = soap.mismatch-100;
			soap.read = soap.read.substr(0,indel_pos) + soap.read.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			soap.qual = soap.qual.substr(0,indel_pos) + soap.qual.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			//cerr<<soap<<endl;
		}

		// add by Bill 2010-10-11
		if (soap.read.size() != soap.read_len)
		{
			soap.read_len = soap.read.size();
		}

		soap.position -= 1;
		return alignment;
	}
	friend std::ostream & operator<<(std::ostream & o, Soap_format & soap) 
	{
		o<<soap.read_id<<'\t'<<soap.read<<'\t'<<soap.qual<<'\t'<<soap.hit<<'\t'<<soap.ab<<'\t'<<soap.read_len<<'\t'<<soap.strand<<'\t'<<soap.chr_name<<'\t'<<soap.position<<'\t'<<soap.mismatch;
		return o;
	}
	char get_base(std::string::size_type coord) 
	{
		return read[coord];
	}
	char get_qual(std::string::size_type coord) 
	{
		return qual[coord];
	}
	bool is_fwd(){
		return (strand=='+');
	}
	int get_read_len(){
		return read_len;
	}
	inline int get_pos()
	{
		return position;
	}
	std::string get_chr_name()
	{
		return chr_name;
	}
	std::string get_read_id()
	{
		return read_id;
	}
	int get_hit()
	{
		return hit;
	}
	bool is_unique()
	{
		return (hit==1);
	}
	bool is_N(int coord) 
	{
		return (read[coord] == 'N');
	}
};

// dbSNP information
class Snp_info 
{
	bool validated;
	bool hapmap_site;
	bool indel_site;
	rate_t * freq; // Arrary of 4 elements recording frequency of ACTG
public:
	Snp_info()
	{
		validated=hapmap_site=indel_site=false;
		freq = new rate_t [4];
		memset(freq,0,sizeof(rate_t)*4);
	}
	Snp_info(const Snp_info & other) 
	{
		validated = other.validated;
		hapmap_site = other.hapmap_site;
		indel_site = other.indel_site;
		freq = new rate_t [4];
		memcpy(freq, other.freq, sizeof(rate_t)*4);
	}
	~Snp_info(){
		delete [] freq;
	}
	friend std::istringstream& operator>>(std::istringstream & s, Snp_info & snp_form)
	{
		s>>snp_form.hapmap_site>>snp_form.validated>>snp_form.indel_site>>snp_form.freq[0]>>snp_form.freq[1]>>snp_form.freq[2]>>snp_form.freq[3];
		return s;
	}
	Snp_info & operator=(Snp_info& other) 
	{
		this->validated = other.validated;
		this->hapmap_site = other.hapmap_site;
		this->indel_site = other.indel_site;
		this->freq = new rate_t [4];
		memcpy(this->freq, other.freq, sizeof(rate_t)*4);
		return *this;

	}
	bool is_validated()
	{
		return validated;
	}
	bool is_hapmap()
	{
		return hapmap_site;
	}
	bool is_indel()
	{
		return indel_site;
	}
	rate_t get_freq(char bin_base_2bit)
	{
		return freq[bin_base_2bit];
	}
};

// Chromosome(Reference) information
class Chr_info 
{
	ubit32_t len;
	ubit64_t* bin_seq; // Sequence in binary format
	ubit64_t* region_mask;
	ubit64_t* region_win_mask;
	// 4bits for one base: 1 bit dbSNPstatus, 1bit for N, followed two bit of base A: 00, C: 01, T: 10, G:11,
	// Every ubit64_t could store 16 bases
	map<ubit64_t, Snp_info*> dbsnp;
	//update by guyue 11-26 the first postion not N 
	int m_start_position ; 
public:
	// add by Bill.
	unsigned m_region_len;
	Chr_info()
	{
		len = 0;
		m_region_len = 0;
		bin_seq = NULL;
		region_mask = NULL;
		region_win_mask = NULL;
		m_start_position = 0; 
	};
	Chr_info(const Chr_info & other);
	~Chr_info(){
		delete [] bin_seq;
		delete [] region_mask;
		delete [] region_win_mask;
	}
	ubit32_t length() {
		return len;
	}
	ubit64_t get_bin_base(std::string::size_type pos)
	{
		return (bin_seq[pos/capacity]>>(pos%capacity*4))&0xF; // All 4 bits
	}
	int binarize(std::string & seq);
	int insert_snp(std::string::size_type pos, Snp_info & new_snp);
	int region_mask_ini();

	/* 
		be changed by Bill
		date : 2010.7.23
		added the judgement
	*/
	bool is_in_region(std::string::size_type pos)
	{
		if (pos > len)
			return 0;
		return ((region_mask[pos/64]>>(63-pos%64))&1);
	}
	bool is_in_region_win(std::string::size_type pos)
	{
		if (region_win_mask == NULL) {
			return 0;
		}
		pos /= global_win_size; // Calculate in which windows the site is
		//cerr<<pos<<endl;
		//exit(1);
		return ((region_win_mask[pos/64]>>(63-pos%64))&1);
	}
	int set_region(int start, int end);
	Snp_info * find_snp(ubit64_t pos) 
	{
		return dbsnp.find(pos)->second;
	}
	ubit64_t * get_region()
	{
		return region_mask;
	}

	// return the non N position
	int getStartPos(void)
	{
		return m_start_position;
	}
};

typedef std::string Chr_name;

class Genome 
{
public:
	map<Chr_name, Chr_info*> chromosomes;

	Genome(my_ifstream & fasta, my_ifstream & known_snp);
	~Genome();

	bool add_chr(Chr_name &);
	int read_region(my_ifstream & region, Parameter * para);
};

class Prob_matrix 
{
public:
	rate_t *p_matrix, *p_prior; // Calibration matrix and prior probabilities
	rate_t *base_freq, *type_likely, *type_prob; // Estimate base frequency, conditional probability, and posterior probablity
	rate_t *p_rank, *p_binom; // Ranksum test and binomial test on HETs
	//update11-11
	bool ran_sum_mode;
	Prob_matrix(bool flag);
	~Prob_matrix();
	int matrix_gen(igzstream & alignment, Parameter * para, Genome * genome);
	int matrix_read(std::fstream & mat_in, Parameter * para);
	int matrix_write(std::fstream & mat_out, Parameter * para);
	int prior_gen(Parameter * para);
	int rank_table_gen();

	// new function. developed by Bill Tang
	int matrix_gen(SamCtrl &alignment, Parameter * para, Genome * genome);
	void count_qual(ubit64_t *count_matrix, Parameter *para);
	void deal_reads(ubit64_t *count_matrix, Genome *genome, Soap_format &soap, map<Chr_name, Chr_info*>::iterator &current_chr);

};
//update 2010-12-20 by guyue
class Pos_info
{
public:
	unsigned char ori;
	//small_int *base_info;
	//update by guyue
	//int *count_uni, *q_sum, depth, dep_uni, repeat_time, *count_all
	small_int base_info[4*2*64*256];
	int depth, dep_uni, repeat_time;
	int count_uni[4], q_sum[4], count_all[4];
	// add by bill,update by guyue
	//int *count_sfs;
	int count_sfs[4];
	int pos;

	Pos_info()
	{
		ori = 0xFF;
		//base_info = new small_int [4*2*64*256]; // base info : 4x2x64x64 matrix, base x strand x qual x read_pos
		memset(base_info,0,sizeof(small_int)*4*2*64*256);
		pos = -1;
		//count_uni = new int [4]; // Count of unique bases
		memset(count_uni,0,sizeof(int)*4);
		//q_sum = new int [4]; // Sum of quality of unique bases
		memset(q_sum,0,sizeof(int)*4);
		depth = 0;
		dep_uni = 0;
		repeat_time = 0;
		//count_all = new int [4]; // Count of all bases
		memset(count_all,0,sizeof(int)*4);
		// add by bill
		//count_sfs = new int [4]; // Count of bases that used to sfs.
		memset(count_sfs,0,sizeof(int)*4);
	}
	~Pos_info()
	{
		//delete [] base_info;
		//delete [] count_uni;
		//delete [] q_sum;
		//delete [] count_all;
		//// add by bill
		//delete [] count_sfs;
	}
};
 
//update 12-27 bu guyue
class Call_win 
{
public:
	ubit64_t win_size;
	ubit64_t read_len;
	Pos_info * sites;
	//update 12-20
	int pos_size;
	int last_start;
	bool recycled;
	//update 11-23
	bool done_pro_win;
	// update 11-26 judge if the laststart had been set.
	bool m_is_set_ls;
	//update 12-27 
	std::string::size_type coord;
	small_int k;
	small_int base_i;
	ubit64_t o_base, strand;
	char allele1, allele2, genotype, type, type1/*best genotype*/, type2/*suboptimal genotype*/, base1, base2, base3;
	int i, q_score, q_adjusted, qual1, qual2, qual3, q_cns, all_count1, all_count2, all_count3;
	int global_dep_count;
	int *pcr_dep_count;
	//pcr_dep_count = new int [para->read_length*2];
	double  rank_sum_test_value, binomial_test_value;
	bool is_out;
	double real_p_prior[16];
	double likelihoods[10];

	map<Chr_name, Chr_info*>::iterator current_chr;
	Call_win(ubit64_t read_length, ubit64_t window_size=global_win_size) 
	{
/*		if (win_size < read_length)
		{
			win_size = (read_length / win_size + 1) * win_size;
			global_win_size = win_size;
		}
		{
			win_size = window_size;
		}*/
		//update 12-20
		pos_size = sizeof(Pos_info);
		//update 11-25
		done_pro_win = false;
		m_is_set_ls = false;
		//update 12-27
		pcr_dep_count = new int [read_length*2];
		for(int i =0; i<16; i++)
		{
			real_p_prior[i] = 0;
		}
		for(int i =0; i<10; i++)
		{
			likelihoods[i] = 0;
		}

		win_size = window_size;
		sites = new Pos_info [win_size+read_length];
		read_len = read_length;
		last_start = 0;
		recycled = false;

		//count_pro_win = 0;
	}
	~Call_win()
	{
		delete [] sites;
		delete [] pcr_dep_count;
	}

	int initialize(ubit64_t start);
	int deep_init(ubit64_t start);
	int recycle();
	int call_cns(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, gzoutstream * consensus, my_ofstream & baseinfo, SfsMethod &sfsMethod, const int id);
	int soap2cns(igzstream & alignment, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id);
	int snp_p_prior_gen(double * real_p_prior, Snp_info* snp, Parameter * para, char ref);
	double rank_test(Pos_info & info, char best_type, double * p_rank, Parameter * para);
	double normal_value(double z);
	double normal_test(int n1, int n2, double T1, double T2);
	double table_test(double *p_rank, int n1, int n2, double T1, double T2);

	// new function. developed by Bill Tang
	int soap2cns(SamCtrl &alignment, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id);
	void pro_win(gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter *para, SfsMethod &sfsMethod, const int id);
	void deal_read(Soap_format &soap, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id);

	// processed tail of chromosome
	virtual void deal_tail(gzoutstream * consensus, my_ofstream& baseinfo, Genome* genome, Prob_matrix* mat, Parameter* para, SfsMethod &sfsMethod, const int id);
};

#endif /*SOAP_SNP_HH_*/

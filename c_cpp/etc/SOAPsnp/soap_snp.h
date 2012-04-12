#ifndef SOAP_SNP_HH_
#define SOAP_SNP_HH_
#include <iostream>
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
typedef unsigned long long ubit64_t;
typedef unsigned int ubit32_t;
typedef double rate_t;
typedef unsigned char small_int;
using namespace std;
const size_t capacity = sizeof(ubit64_t)*8/4;
const char abbv[17]={'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};
const ubit64_t glf_base_code[8]={1,2,8,4,15,15,15,15}; // A C T G
const ubit64_t glf_type_code[10]={0,5,15,10,1,3,2,7,6,11};// AA,CC,GG,TT,AC,AG,AT,CG,CT,GT
const int global_win_size = 1000;

// Some global variables
class Files {
public:
	ifstream soap_result, ref_seq, dbsnp, region;
	ofstream consensus, baseinfo, o_region;
	fstream matrix_file;
	Files(){
		soap_result.close();
		ref_seq.close();
		dbsnp.close();
		consensus.close();
		baseinfo.close();
		matrix_file.close();
		region.close();
		o_region.close();
	};
};

class Parameter {
public:
	char q_min; // The char stands for 0 in fastq
	char q_max; // max quality score
	small_int read_length; // max read length
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
// Default onstruction
	Parameter(){
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
	};
};

class Soap_format {
	// Soap alignment result
	std::string read_id, read, qual, chr_name;
	int hit, read_len, position, mismatch;
	char ab, strand;
public:
	Soap_format(){;};
	friend std::istringstream & operator>>(std::istringstream & alignment, Soap_format & soap) {
		alignment>>soap.read_id>>soap.read>>soap.qual>>soap.hit>>soap.ab>>soap.read_len>>soap.strand>>soap.chr_name>>soap.position>>soap.mismatch;
		//cerr<<soap<<endl;
		//exit(1);
		if(soap.mismatch>200) {
			int indel_pos,indel_len;
			string temp("");
			alignment>>indel_pos;
			indel_len = soap.mismatch-200;
			for(int i=0; i!=indel_len; i++) {
				temp = temp+'N';
			}
			soap.read = soap.read.substr(0,indel_pos)+temp+soap.read.substr(indel_pos,soap.read_len-indel_pos);
			soap.qual = soap.qual.substr(0,indel_pos)+temp+soap.qual.substr(indel_pos,soap.read_len-indel_pos);
			//cerr<<soap<<endl;
		}
		else if (soap.mismatch>100) {
			int indel_pos,indel_len;
			alignment>>indel_pos;
			indel_len = soap.mismatch-100;
			soap.read = soap.read.substr(0,indel_pos) + soap.read.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			soap.qual = soap.qual.substr(0,indel_pos) + soap.qual.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			//cerr<<soap<<endl;
		}
		//cerr<<soap.position<<'\t';
		soap.position -= 1;
		//cerr<<soap.position<<endl;
		return alignment;
	}
	friend std::ostream & operator<<(std::ostream & o, Soap_format & soap) {
		o<<soap.read_id<<'\t'<<soap.read<<'\t'<<soap.qual<<'\t'<<soap.hit<<'\t'<<soap.ab<<'\t'<<soap.read_len<<'\t'<<soap.strand<<'\t'<<soap.chr_name<<'\t'<<soap.position<<'\t'<<soap.mismatch;
		return o;
	}
	char get_base(std::string::size_type coord) {
		return read[coord];
	}
	char get_qual(std::string::size_type coord) {
		return qual[coord];
	}
	bool is_fwd(){
		return (strand=='+');
	}
	int get_read_len(){
		return read_len;
	}
	inline int get_pos(){
		return position;
	}
	std::string get_chr_name(){
		return chr_name;
	}
	int get_hit(){
		return hit;
	}
	bool is_unique(){
		return (hit==1);
	}
	bool is_N(int coord) {
		return (read[coord] == 'N');
	}
};

// dbSNP information
class Snp_info {
	bool validated;
	bool hapmap_site;
	bool indel_site;
	rate_t * freq; // Arrary of 4 elements recording frequency of ACTG
public:
	Snp_info(){
		validated=hapmap_site=indel_site=false;
		freq = new rate_t [4];
		memset(freq,0,sizeof(rate_t)*4);
	}
	Snp_info(const Snp_info & other) {
		validated = other.validated;
		hapmap_site = other.hapmap_site;
		indel_site = other.indel_site;
		freq = new rate_t [4];
		memcpy(freq, other.freq, sizeof(rate_t)*4);
	}
	~Snp_info(){
		delete [] freq;
	}
	friend std::istringstream& operator>>(std::istringstream & s, Snp_info & snp_form) {
		s>>snp_form.hapmap_site>>snp_form.validated>>snp_form.indel_site>>snp_form.freq[0]>>snp_form.freq[1]>>snp_form.freq[2]>>snp_form.freq[3];
		return s;
	}
	Snp_info & operator=(Snp_info& other) {
		this->validated = other.validated;
		this->hapmap_site = other.hapmap_site;
		this->indel_site = other.indel_site;
		this->freq = new rate_t [4];
		memcpy(this->freq, other.freq, sizeof(rate_t)*4);
		return *this;

	}
	bool is_validated(){
		return validated;
	}
	bool is_hapmap(){
		return hapmap_site;
	}
	bool is_indel(){
		return indel_site;
	}
	rate_t get_freq(char bin_base_2bit) {
		return freq[bin_base_2bit];
	}
};

// Chromosome(Reference) information
class Chr_info {
	ubit32_t len;
	ubit64_t* bin_seq; // Sequence in binary format
	ubit64_t* region_mask;
	ubit64_t* region_win_mask;
	// 4bits for one base: 1 bit dbSNPstatus, 1bit for N, followed two bit of base A: 00, C: 01, T: 10, G:11,
	// Every ubit64_t could store 16 bases
	map<ubit64_t, Snp_info*> dbsnp;
public:
	Chr_info(){
		len = 0;
		bin_seq = NULL;
		region_mask = NULL;
		region_win_mask = NULL;
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
	ubit64_t get_bin_base(std::string::size_type pos) {
		return (bin_seq[pos/capacity]>>(pos%capacity*4))&0xF; // All 4 bits
	}
	int binarize(std::string & seq);
	int insert_snp(std::string::size_type pos, Snp_info & new_snp);
	int region_mask_ini();
	bool is_in_region(std::string::size_type pos) {
		return ((region_mask[pos/64]>>(63-pos%64))&1);
	}
	bool is_in_region_win(std::string::size_type pos) {
		pos /= global_win_size; // Calculate in which windows the site is
		//cerr<<pos<<endl;
		//exit(1);
		return ((region_win_mask[pos/64]>>(63-pos%64))&1);
	}
	int set_region(int start, int end);
	Snp_info * find_snp(ubit64_t pos) {
		return dbsnp.find(pos)->second;
	}
	ubit64_t * get_region() {
		return region_mask;
	}
};

typedef std::string Chr_name;
class Genome {
public:
	map<Chr_name, Chr_info*> chromosomes;

	Genome(ifstream & fasta, ifstream & known_snp);
	~Genome();

	bool add_chr(Chr_name &);
	int read_region(std::ifstream & region, Parameter * para);
};

class Prob_matrix {
public:
	rate_t *p_matrix, *p_prior; // Calibration matrix and prior probabilities
	rate_t *base_freq, *type_likely, *type_prob; // Estimate base frequency, conditional probability, and posterior probablity
	rate_t *p_rank, *p_binom; // Ranksum test and binomial test on HETs
	Prob_matrix();
	~Prob_matrix();
	int matrix_gen(std::ifstream & alignment, Parameter * para, Genome * genome);
	int matrix_read(std::fstream & mat_in, Parameter * para);
	int matrix_write(std::fstream & mat_out, Parameter * para);
	int prior_gen(Parameter * para);
	int rank_table_gen();

};

class Pos_info {
public:
	unsigned char ori;
	small_int *base_info;
	int pos, *count_uni, *q_sum, depth, dep_uni, repeat_time, *count_all;

	Pos_info(){
		ori = 0xFF;
		base_info = new small_int [4*2*64*256]; // base info : 4x2x64x64 matrix, base x strand x qual x read_pos
		memset(base_info,0,sizeof(small_int)*4*2*64*256);
		pos = -1;
		count_uni = new int [4]; // Count of unique bases
		memset(count_uni,0,sizeof(int)*4);
		q_sum = new int [4]; // Sum of quality of unique bases
		memset(q_sum,0,sizeof(int)*4);
		depth = 0;
		dep_uni = 0;
		repeat_time = 0;
		count_all = new int [4]; // Count of all bases
		memset(count_all,0,sizeof(int)*4);
	}
	~Pos_info(){
		delete [] base_info;
		delete [] count_uni;
		delete [] q_sum;
		delete [] count_all;
	}
};

class Call_win {
public:
	ubit64_t win_size;
	ubit64_t read_len;
	Pos_info * sites;
	Call_win(ubit64_t read_length, ubit64_t window_size=global_win_size) {
		sites = new Pos_info [window_size+read_length];
		win_size = window_size;
		read_len = read_length;
	}
	~Call_win(){
		delete [] sites;
	}

	int initialize(ubit64_t start);
	int deep_init(ubit64_t start);
	int recycle();
	int call_cns(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, std::ofstream & consensus, std::ofstream & baseinfo);
	int soap2cns(std::ifstream & alignment, std::ofstream & consensus, std::ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para);
	int snp_p_prior_gen(double * real_p_prior, Snp_info* snp, Parameter * para, char ref);
	double rank_test(Pos_info & info, char best_type, double * p_rank, Parameter * para);
	double normal_value(double z);
	double normal_test(int n1, int n2, double T1, double T2);
	double table_test(double *p_rank, int n1, int n2, double T1, double T2);
};

#endif /*SOAP_SNP_HH_*/

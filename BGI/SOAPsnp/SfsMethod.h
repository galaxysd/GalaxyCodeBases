/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-9-13
 *CLASS NAME: SfsMethod
 *FUNCTION : a class that contain method about sfs.
 *FILE NAME : SfsMethod.h
 *UPDATE DATE : 2010-9-13
 *UPDATE BY : Bill Tang
 *UPDATE DATE : 2010-10-13
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */

#pragma once

#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include "assert.h"
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <semaphore.h>
//#include <string.h>

using namespace std;

class Files;
class Parameter;
class Pos_info;
class Prob_matrix;

/*soaplikelihood to sfslikelihood*/
const int sfs_type[10]={0,1,3,2,5,7,6,15,11,10}; // AA,AC,AG,AT,CC,CG,CT,GG,GT,TT

typedef struct _loci{
	string chromo;
	int position;
/*	bool operator <(const _loci lc2) const
	{
		if (chromo == lc2.chromo)
			return position < lc2.position;
		else
			return chromo < lc2.chromo;
	}
	*/
}loci;
/*
struct loci_hash
{
	size_t operator()(const loci& lc) const
	{
		return size_t(lc.chromo + lc.position);
	}
};*/

/*
This struct is all the data we need to have for each site for all individuals

*/
typedef struct {
	int* locus; //4 times number of individuals {A,C,T,G} +4 the last 4 will be the sum
	//unsigned char  *lk; //3 times the number of individuals {AA,Aa,aa}
	int  *lk; //10 times the number of individuals {AA,Aa,aa}
	short int major; //either ACGT {0...3} a
	short int minor; //either ACGT {0...3} A
	int ref;
	//it might be overzealous to keep this in the datum...
	//but it might prove handy for debug purposes
	double phat;
	//  int numHits; //counter for keeping track of which files has data for a loci
//	int position;	// allele's position.
//	std::string chromo;	// the chromone ID.
//	bool is_be_record; // record that if the pos have bean in.
}datum;

//comparison operator used when comparing loci
struct cmp_loci {
	bool operator()(const loci& first,const  loci& second) const {
		if(first.chromo==second.chromo)
		{
			return first.position<second.position;
		}
		else
		{
			return first.chromo<second.chromo;
		}
	}
};
/*
struct cmp_loci {
	bool operator()(const loci& first,const  loci& second) const {
		return ((first.chromo == second.chromo) && (first.position == second.position));
	}
};
*/
struct cmp_char {
	bool operator()(const char* first,const  char* second) const {
		int tmp = std::strcmp(first, second);
		return tmp<0;
	}
};

//typedef boost::unordered_map<loci, datum, loci_hash, cmp_loci> aMap; // mapping of the reads per locus
typedef std::map<loci, datum, cmp_loci> aMap; // mapping of the reads per locus
typedef std::map<char *,int,cmp_char> cMap; // mapping of the chromoname -> int id
typedef datum	aVector;	// the vector that store the datum.

class SfsMethod
{
public:
	SfsMethod();
	SfsMethod(int numInds);
	virtual ~SfsMethod(void);
	// the main function for the sfs, update 11-29 add int sfsfirstlast
	virtual void algo(aMap & asso, int numInds, double eps, double pvar, double B, FILE * sfsfile, int normalize, int alternative, int singleMinor, double pThres, int allowPhatZero, int firstlast);
	virtual void algoJoint(aMap & asso, int numInds, double eps, double pvar, double B, FILE * sfsfile, int underFlowProtect, int alternative, int singleMinor, int fold);
	// write the freqfile
	void writeFreq(FILE * pFile, int numInds, aMap & asso);
	// allocation functions, this shouldn't look magical for anyone
	virtual double* allocDoubleArray(int len);
	//update 11-16
	sem_t sem_call_sfs;
	//sem_t sem_map_change;
	sem_t sem_call_sfs_return ;
	//sem_t sem_map_number;
	int file_end_flag;
private:
	// genome likelihood Format look up table
	int glfLookup[4][4];
	//used for offsetting each frame in the array of loglikeratios for a site
	enum {AA = 0, Aa = 1, aa = 2};
	//used for offsetting each frame in the array of basereads for a site
	enum {A = 0, C = 1, G = 2, T = 3};
	//to speed up calculations, we will generate a likelookup table for all values 0-255.
	double likeLookup[256];
	//The index will then be the value from one of the columns in the glf file
	double *binomLookup;
	//this is the result array
	double *bayes;
	double *likes;
	// aMap pointer.
	aMap asso[3];
	//the map index
	int m_map_idx;
	int m_map_idx_process; // the index of the process map
	// aVector *asso;
	// individual number
	int m_numInds;
	// the chr index map.
	cMap m_faiIndex;
	// the thread synchronization lock.
	static pthread_mutex_t m_pthreadMutex;
public:
	// test for infinity
	int isinf_local(double x);
	// test for a NaN
	int isnan_local(double x);
	// add up all the value form the ary
	double sum(const double* ary, int len);
	// init the like look up table.
	void generate_likeLookup();
	// count ratio
	double calcLikeRatio(int a)
	{
		return (pow(10.0,-a/10.0));
	}
	// get Max id in the array d.
	int getMaxId(double *d,int len);
	// input should be number of individuals
	void generate_binom(int k);
	// FUNCTION: generate factln
	double factln(int n);
	// To calculate the nCk ( n combination k)
	double bico(int n, int k);
	// generate gammln
	double gammln(double xx);
	// clean map
	void cleanUpMap(aMap& asso);
	// get char from int
	char getChar(int i);
	// add protect
	double addProtect(double *ary,int len);
	// add three protect
	double addProtect3(double a,double b, double c);
	// add two protect
	double addProtect2(double a,double b);
	// print array
	void print_array(double *ary,int len);
	// alloc datum memery.
	datum allocDatum(int numInds);
	// alloc array memery.
	int* allocIntArray(int len);
	// will calculate the sums and input the values in the last 4 slots, and set the minor, major accordingly
	void calcSumBias(aMap & asso,int numInds);
	// FUNCTION:  soaplikelihood to sfslikelihood
	void soaplk_sfslk(int *likelihood, double *type_likely, int id);
	//copy source's coverage to the dest
	void cov_Cpy(int *dest, int const *source, int id);
	// alloc map memery
	virtual void allocVec(void);
	// free the map memery
	virtual void delMap(aMap &amap);
	// do SFS
	int call_SFS(Parameter* para, Files* files);
	// get the data what sfs want
	int getMapData(const Pos_info& site, const Prob_matrix* mat, const std::string & chr, const int id);
	// initialization
	int init(const int numInds);
	// This builds a map from: char* -> int
	virtual void buildMap(const char* fname);
	// free the vector's data memery
	virtual void delVector(aVector * avec);
	// get data what sfs want
	virtual int getSFSData(const Pos_info & site, const Prob_matrix * mat, const std::string & chr, const int id);
	// clean up the vector's member.
	void cleanVec(void);
	//update 11-16
	// change map
	void mapChange(void);
	//get the map index can be process
	int getidxProcess(void);
	//wait the map to be process
	void setidxProcess(void);
};

//the call sfs structor
typedef struct _big_call_sfs_args
{
	SfsMethod * sfsMethod;
	Files  *files;
	Parameter * para;

	inline _big_call_sfs_args(SfsMethod * a, Files * b, Parameter * c)
	{
		sfsMethod = a;
		files = b;
		para = c;
	};

	inline _big_call_sfs_args()
	{
		sfsMethod = NULL;
		files = NULL;
		para = NULL;
	};

}BIG_CALL_SFS_ARGS;

void *_sfsMethod_callsfs(void * args);

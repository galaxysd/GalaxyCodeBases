/*****************************************************************
 *****************************************************************/


#ifndef  __Tpye_C_
#define  __Tpye_C_

/*typedef unsigned int ubit32_t;*/

#define MAX_FILENAME_LEN 1024
#define READ_NAME_LEN 256
#define SMALL_QUERY_ALLOCATION_UNIT	65536
#define LARGE_QUERY_ALLOCATION_UNIT	1048576
#define BATCH_NUM 10000
#define MAX_SUFFIX_LEN 32
#define MAX_HITS_NUM 65535
#define INI_BUF_SIZE 10000
#define MAX_MISMATCH 2
#define MAX_GAP 3

#define FORWARD 1
#define REVERSE 0
#define SUBSTITUTION 0  // 0 1 2 3 A C G T

enum fileformat{FASTA = 0,
	FASTQ 
};

typedef struct _ReadInfo_{
	unsigned id;
	unsigned nameLength	:16;
	unsigned seqLength	:16;
	unsigned numOfExact	:16;
	unsigned numOf1Mis	:16;
	unsigned numOf2Mis	:16;
	unsigned numOfHits	:16;
	unsigned mm_n;
	char name[READ_NAME_LEN];
	char file;
	char *seq;
	char *qual;
	char *revSeq;
	char *revQual;
}ReadInfo;

typedef struct {
	unsigned int : 12;
}uint12;

typedef struct _HitInfo_{
	int chrID;
	unsigned int pos;
	uint12 mm[MAX_MISMATCH];
	unsigned int l		:12;		// 1st mismatch or start of gap		----------  --
	unsigned int r 		:12;		// 2nd mismatch or end of gap		|    pos  |base|
	unsigned int type 	:4;		// 0: exact, 1: mismatch, 2: ins; 3: del;
	unsigned int chain 	:1;		// 1: f; 0: r;
	unsigned int file 	:1;		// 0: a; 1: b;
	unsigned int n_mm;
}HitInfo;

typedef void (*FormattedOutput)(ReadInfo *, HitInfo *, const unsigned int, const char *, FILE *);
typedef int (*LoadSeq)(FILE *, ReadInfo *);

typedef struct _Para_{
	char QueryFileName_A[MAX_FILENAME_LEN+1];
	char QueryFileName_B[MAX_FILENAME_LEN+1];
	char DataBaseName[MAX_FILENAME_LEN+1];
	char OutputFileName[MAX_FILENAME_LEN+1];
	char UnmappedFileName[MAX_FILENAME_LEN+1];
	char UnpairedFileName[MAX_FILENAME_LEN+1];
	int OutputFormat;
	int Unmapped;
	int MismatchNumber;
	int GapSize;
	int MaxHitsNum;
	int IgnoreMismatch;
	int GapForbidden;
	int ZeroQual;
	int TrimLowQual;
	int FilterNs;
	int ReportRepeat;
	int Chains;
	int OutputID;
	int ProcsNum;
	int MappingPair;
	int MappingQual;
	int MinInsertSize;
	int MaxInsertSize;
	int Unpaired;
	int MatchMode;
	int Cutoff;
	int Trim5;
	int Trim3;
	FormattedOutput seqout;
}Parameters;

typedef struct _PairInfo_{
	int flag;		//0: readsA	1: readsB
	HitInfo hit_L;
	HitInfo hit_R;
	int insertSize;
	int mismatch;
}PairInfo;

/*
typedef struct _Allele_{
	int pos;
	char allele;
}AlleleInfo;

typedef struct _PosTable_{
	char chrName1[MAX_SEQ_NAME_LENGTH+1];
	unsigned int startPos1;
	char chrName2[MAX_SEQ_NAME_LENGTH+1];
	unsigned int startPos2;
}PosTable;

typedef struct _HitInfo_{
	int chain;
	int gapPos;
	int gapSize;
	int numOfMismatch;
	unsigned int pos;
	unsigned int mappingQual;
	int allelPos[MAX_MISMATCH];
	char alleles[MAX_MISMATCH];
	char gap[MAX_GAP];
}HitInfo;
//*/


#endif /* __Tpye_C_*/

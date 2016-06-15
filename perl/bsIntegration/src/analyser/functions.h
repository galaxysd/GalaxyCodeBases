#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <err.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h> //offsetof
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "klib/khash.h"
#include "klib/kvec.h"
#include "klib/kstring.h"

#include <htslib/sam.h>

typedef struct {
	const char * ProjectID;
	const char * WorkDir;
	const char * RefileName;
	uint16_t minGrepSlen;
	uint16_t minHumMapQ;
	uint8_t samples;
} __attribute__ ((__packed__)) Config_t;

Config_t myConfig;

typedef struct {
	int8_t isHum;
	uint32_t ChrLen;
	//int32_t tid;	// bam1_core_t.tid; https://github.com/samtools/htslib/blob/f58922d3ad64185890cd01cdf99cc649919afcfa/htslib/sam.h#L150
} __attribute__ ((__packed__)) ChrInfo_t;
typedef struct {
	uint16_t insertSize;
	uint16_t SD;
	const char * fileName;
} __attribute__ ((__packed__)) BamInfo_t;

KHASH_INIT(chrNFO, kh_cstr_t, ChrInfo_t, 1, kh_str_hash_func, kh_str_hash_equal)
KHASH_INIT(bamNFO, kh_cstr_t, BamInfo_t, 1, kh_str_hash_func, kh_str_hash_equal)

khash_t(chrNFO) *chrNFOp;	// kh_##name##_t -> kh_chrNFO_t
khash_t(bamNFO) *bamNFOp;
//kh_chrNFO_t *chrNFOp = kh_init(chrNFO);
/*
error: initializer element is not a compile-time constant
全局变量是保存在静态存储区的，因此在编译的时候只能用常量进行初始化，而不能用变量进行初始化。需要在执行时确定（如：在main函数里赋值）
g++编译器会先把全局变量保存到.bss段中，而且默认值为0，但是会在main函数之前添加一条赋值语句，也就是相当于局部变量进行处理了。
 */

typedef struct {
	int32_t tid;
	int32_t pos;
	int32_t endpos;
	uint8_t qual;
} __attribute__ ((__packed__)) chrRange_t;
typedef struct {
	chrRange_t HumanRange;
	chrRange_t VirusRange;
	kvec_t(uint8_t) quals;
	kvec_t(bam1_t) Reads;
} __attribute__ ((__packed__)) pierCluster_t;
// int32_t bam_endpos(const bam1_t *b); or overlap_push()

typedef struct {
	bam1_t b;
	const char * sam1txt;
	//int8_t has2ndhit;
	int32_t beg, end;
	uint8_t qual;
} __attribute__ ((__packed__)) samdat_t;

//int do_grep(const Config_t *pConfig, const khash_t(chrNFO) *chrNFOp, const khash_t(bamNFO) *bamNFOp);
int do_grep();
int do_analyse();

int getPairedSam(htsFile *fp, hts_idx_t *idx, bam1_t *b, bam1_t *d);
int checkMapQ(int8_t *ChrIsHum, bam1_t *b, bool save_tid);

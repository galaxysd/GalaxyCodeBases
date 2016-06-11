#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <stddef.h> //offsetof
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "klib/khash.h"

typedef struct {
	const char * ProjectID;
	const char * WorkDir;
	const char * RefileName;
} __attribute__ ((__packed__)) Config_t;

Config_t myConfig;

typedef struct {
	int8_t isHum;
	uint32_t ChrLen;
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

//int do_grep(const Config_t *pConfig, const khash_t(chrNFO) *chrNFOp, const khash_t(bamNFO) *bamNFOp);
int do_grep();
int do_analyse();


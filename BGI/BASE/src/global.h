#ifndef __GLOBAL__H_
#define __GLOBAL__H_
/*

   global.h        global variables in BASE.

#    Copyright (C) 2015, The University of Hong Kong.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 3
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  03110-1301, USA.

    Date   : 28th Jan 2016
    Author : Binghang Liu, bhliu@cs.hku.hk
*/

#include <pthread.h>
#include <unistd.h>  //usleep

#define MAX_BIT_LEN 8

typedef unsigned long long ULL;

typedef long long LL;

struct Contig
{
    char* seq;
    char* qual;
    char* uniq;
    ULL len;
    ULL id;
    int s_len;
    int s_depth;
    int back_end;
    int forw_end;
    int back_iter;
    int forw_iter;
    int back_rc;
    int forw_rc;
};

struct Id
{
    ULL c;
    ULL size;
};

struct EndInf
{
    ULL ur; //used read number.
    ULL c; //total simple count.
    ULL len; //total extend length;
    ULL uc; //used count.
    ULL ulen; //total used extension length, length < read_length are not used.
};

struct Base  //72=8*9.
{
    ULL saR;
    ULL saL;
    ULL saRR;
    ULL saLL;
    ULL saRRrev;
    ULL saLLrev;
    ULL plus;
    ULL minus;
    ULL depth;
};

struct Node //88=72+8+2+2+1+1+1+1
{
    unsigned short len:14;      //sequence length.
    unsigned short level:15;    //last level number.
    unsigned short is_end:1;   //no child.
    unsigned short is_extend:1;//is extend.
    unsigned short parent;  //parent node id.
    unsigned short conf:1;
    unsigned short depth;   //initial depth when creating the node, the sum depth of all its children.
    char* seq;              //node sequence.
    Base* b;                 //SA range.
};

struct Lev  //16 = 8 + 2 + 2 + 2 + 2.
{
    unsigned short depth; //depth for consensus allele.
    unsigned short base;  //consensus allele.
    unsigned short node;  //majority node id.
    unsigned short nnum;  //nodes number in current level.
    unsigned short conf;  //confliction starting level.
    unsigned short* c;    //count for 0, 1, 2, 3, 4.
};

struct ReadPair //1.
{
    unsigned char flag:1;  //used or not.
    unsigned char bi:1;    //fully used or not.
    unsigned char s:1;     //read1 or read2.
    unsigned char pflag:1; //pair read used or not.
    unsigned char curr:1;  //pair read used in current extending.
    unsigned char size:3;  //order of insert size, <8 insert size libraries.
};

struct TmpRead
{
    ULL id:36;          //read id.
    ULL size:3;         //order of insert size.
    ULL strand:1;       //read strand in this contig.
    ULL unique:6;
    ULL iter:16;
    ULL cons:1;          //reads is same to consensus(0) or not(1).
    ULL keep:1;         //set bwt flags or not.
    LL pos;             //read start position in contig.
};

struct ReadId
{
    ULL id:36;          //read id.
    ULL strand:1;       //read strand in this iteration.
    ULL level:9;        //read end level number, the level meet $.
    ULL node:15;        //read last node number.
    ULL cons:1;         //read is same to consensus(0) or not(1);
    ULL size:2;         //order of insert size.
};

struct Iter
{
    unsigned short  repeat:1;    //seed is repeat.
    unsigned short  unique:1;    //seed is unique.
    unsigned short  pe:6;
    unsigned short  len:9;   //extending length in this iteration.
    unsigned short  end:6;   //end type in this iteration.
    unsigned short  slen:9;  //length for this seed.
    unsigned short rnum;  //read number found in this iteration.
    unsigned short depth; //starting depth for this iteration.
    char* seq;            //seed for this iteration.
    Base* b;               //starting SA range for this iteration.
};

extern int depth_3_length;
extern int depth_6_length;
extern int depth_9_length;
extern int depth_half_length;

extern unsigned char charMap[256];
extern unsigned char reverseCharMap[CHAR_MAP_SIZE];
extern unsigned char complementMap[256];

extern float Upper_bound;
extern int Expect_depth;
extern int Error_depth;
extern int Solve_Hete;
extern int Solve_Conf;

extern int MINOVERLAP;
extern int Read_length;
extern int Consensus_depth;
extern double MINRATEDIF;
extern double MINLENDIF;
extern double MAXDIFRATIO;
extern int MAXDIFNUM;

extern double Similar_cutoff;

extern int MaxIterNum;
extern int MaxCtgLen;
extern int MaxReadCount;

extern Idx2BWT** idx2BWT;
extern int bwtCount;
extern BWT **bwts;
extern BWT **rev_bwts;
extern int* InsertSize;
extern int* SD;
extern int* FLAGS;

extern int StopReadNum;
extern ULL UsedReadNum;

extern int MAXDEPTH;

extern ULL pair_num;

// For multi-thread version.
//for one extension.
extern TmpRead** tmpRead; 
extern int* tmpCount;
extern Iter** g_iters;
extern int* iterCount;

//for one cycle.
extern char** g_seq;
extern ReadId** g_reads; 
extern Lev** g_levs;
extern Node** g_nodes;

//for each thread, melloc one.
extern Base** tmpbase;
extern Base*** sbase;

//for thread.
extern int* signal;
extern pthread_t * threads;
extern pthread_mutex_t mut;
extern int threadNum;

//for thread buffers.
extern Contig** ctgs; //contigs per thread.
extern int* ctg_num; //contig number per thread.
extern int buffer_size; //output 100 contigs once.
extern Id** conf_buffer;
extern int* conf_size;
extern Id* current_id;
extern int* threadEnd;
extern int g_count;

extern EndInf* endInf;
extern double read_ratio;
extern ULL ass_len;

extern int bitmatrix[64][6];
extern int id_to_sparse_word[64];

extern ULL count_pair;
extern int count_pair_thread;
extern int gthreadNum;

extern int increase;
extern int old_avg_len;
extern int is_platform;

extern int* four_iter_used_read;
extern int* four_iter_checked_read;
extern int* four_iter;

extern int RepeatMask;
extern int DEBUG ;

#endif


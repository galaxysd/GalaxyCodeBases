/*

   main.c        main file of BASE.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include "./2bwt/Timing.h"
#include "./2bwt/TextConverter.h"
#include "contigs.h"

#define APP "BASE"
#define VERSION 1.0
#define RELEASE "28/1/2016"

int DEBUG = 0;

int MaxIndel = 2;
int MINOVERLAP = 19;
int MAXOVERLAP = 81;

double MINRATEDIF = 0.005;
double MINLENDIF = 0.005;
double MAXDIFRATIO = 0.2;
int MAXDIFNUM = 5;

int MINNUM = 10;
unsigned char charMap[256];
unsigned char complementMap[256];
unsigned char reverseCharMap[CHAR_MAP_SIZE];

int Error_depth = 3;
int Consensus_depth=3;
float Upper_bound = 1.2;
int Expect_depth = 10;
int Read_length=100;
int RepeatMask = 0;

int depth_3_length = Read_length-int(double(Error_depth)*(float)Read_length/(float)Expect_depth);
int depth_6_length, depth_9_length, depth_half_length;
double match_ratio = 0.75;
double Similar_cutoff=0.6;

int MAXDEPTH=400; //65535;
int MaxIterNum = 40000;
int MaxCtgLen = 1000000; //1000000; 
int MaxReadCount = 1000000; //10000000;
int Solve_Hete = 0;
int Solve_Conf = 0;
unsigned long long UsedReadNum;

char lib_file[256];
char prefix[256];

Idx2BWT** idx2BWT;
int bwtCount;
BWT **bwts;
BWT **rev_bwts;
unsigned long long pair_num;

int* InsertSize;
int* SD;
int* FLAGS;
int StopReadNum=0;

// For multi-thread version.
//for one extension.
TmpRead** tmpRead; 
int* tmpCount;
Iter** g_iters;
int* iterCount;

//for one cycle.
char** g_seq;
ReadId** g_reads; 
Lev** g_levs;
Node** g_nodes;

//for each thread, melloc one.
Base** tmpbase;
Base*** sbase;

//for thread.
int* signal;
pthread_t * threads;
pthread_mutex_t mut = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
int threadNum=6;

//for thread buffers.
Contig** ctgs; //contigs per thread.
int* ctg_num; //contig number per thread.
int buffer_size = 300; //output 100 contigs once.
Id** conf_buffer;
int* conf_size;
Id* current_id;
int* threadEnd;
int g_count=0;

EndInf* endInf;
double read_ratio=0;
unsigned long long ass_len=0;
unsigned long long count_pair;
int count_pair_thread;
int gthreadNum;
int bitmatrix[64][6];
int id_to_sparse_word[64];
int increase;
int old_avg_len;
int is_platform = 0;

int* four_iter_used_read;
int* four_iter_checked_read;
int* four_iter;

int exists_file(const char* filename)
{
	return (access(filename, 0)== 0); //F_OK =0, 0 is used to check F_OK.
}

void fill_sparse_word()
{
    int i, j;
    unsigned int v, s;
    for(i=0; i< 64; i++)
    {
        s = 0; v = i;
        for (j = 0; j < 6; ++j) {
            s |= (v & 1) << (3 + 4 * j);
            v >>= 1;
        }
        id_to_sparse_word[i] = s;
    }
}

void get_rcseq(char* seq, char* rcseq)
{
        int i;
        for(i=0; i<strlen(seq); i++)
        {
            rcseq[i] = complementMap[ seq[ strlen(seq) - i - 1 ] ];
        }
        rcseq[i]='\0';
}

void Usage(const char* arg_0)
{
    fprintf(stderr, "\nUsage Guide:\n");
    fprintf(stderr, "%s is used to assemble pair end illumina reads.\n", APP);
    fprintf(stderr, "Copyright (c) HKU-BGI Bioinformatics Algorithms and Core Technology Research Laboratory (BAL),\nDepartment of Computer Science, University of Hong Kong\n");
    fprintf(stderr, "Contact: bhliu@cs.hku.hk\n");
    fprintf(stderr, "Version: %.2f; Compile Date: %s, Time: %s\n", VERSION, __DATE__, __TIME__);
    fprintf(stderr, "Usage:\n\t./base <command> [option] >log 2>error \n", APP);
    
    fprintf(stderr, "\t-l <string>\t2BWT library file\n");
    fprintf(stderr, "\t-o <string>\tPrefix of output files\n");
    fprintf(stderr, "\t-E <int>\tExpected depth\n");
    fprintf(stderr, "\t-L <int>\tLow depth threshold, default=%d\n", Error_depth);
    fprintf(stderr, "\t-M <int>\tMinimum overlap, default = %d\n", MINOVERLAP);
    fprintf(stderr, "\t-T <double>\tInferred unique threshold, default=%.2f\n", Upper_bound);
    fprintf(stderr, "\t-P <int>\tThread number, default=%d\n", threadNum);
    fprintf(stderr, "\t-C (optional)\tSolve branches or not, [No]\n");
    fprintf(stderr, "\t-H (optional)\tSolve heterozygosis, [No]\n");
    fprintf(stderr, "\t-h\t\tprint the usage message.\n\n");
}

void loadOptions(int argc, char **argv)
{
    int opt;
    extern char * optarg;
    char temp[128];
    char copt;
    while ( ( copt = getopt ( argc, argv, "l:L:E:T:P:o:M:CDHh" ) ) != EOF )
    {
        switch ( copt )
        {
            case 'L':
                sscanf ( optarg, "%s", temp );
                Error_depth = atoi ( temp );
                break;
            case 'E':
                sscanf ( optarg, "%s", temp);
                Expect_depth = atoi(temp);
                break;
            case 'T':
                sscanf (optarg, "%s", temp );
                Upper_bound = atof(temp);
                break;
            case 'P':
                sscanf (optarg, "%s", temp );
                threadNum = atoi(temp);
                break;
            case 'M':
                sscanf (optarg, "%s", temp );
                MINOVERLAP = atoi(temp);
                break;
            case 'l':
                sscanf (optarg, "%s", lib_file );
                break;
            case 'o':
                sscanf (optarg, "%s", prefix );
                break;
            case 'C':
                Solve_Conf = 1;
                break;
            case 'H':
                Solve_Hete = 1;
                break;
            case 'D':
                DEBUG = 1;
                break;
	    case 'h':
		Usage(argv[0]);
		exit(0);
            default:
                Usage(argv[0]);
                break;
        }
    }
    fprintf(stderr,"\nVersion %.2f: Released at %s and Compilied at %s\t%s\n", VERSION, RELEASE, __DATE__ , __TIME__ );
    fprintf(stderr,"\nMain arguments:\n");
    fprintf(stderr, "Expect depth: %d\nLow depth: %d\nInfer threshold: %.2f\nSolve branches: %d\nSolve heterozygosis: %d\n\n", 
            Expect_depth, Error_depth, Upper_bound, Solve_Conf, Solve_Hete);
}
int main(int argc, char ** argv) {

   if(argc < 2)
    {
        Usage(argv[0]);
        return 1;
    }
    loadOptions(argc, argv);

    time_t start_time, now;
    time(&start_time);
    fprintf(stderr, "[Main] invoked at %s\n", ctime(&start_time));

    fprintf(stderr, "Loading index ... \n");

    bwtCount = load_bwts(lib_file);
    HSPFillCharMap(charMap);
    HSPFillComplementMap(complementMap);
    time(&now);
    fprintf(stderr, "Loading index DONE at %s\n\n", ctime(&now));

    Read_length = idx2BWT[0]->readLength;
    depth_3_length = Read_length-int(double(Error_depth)*(float)Read_length/(float)Expect_depth);
    depth_6_length = Read_length-int(double(Error_depth*2)*(float)Read_length/(float)Expect_depth);
    depth_9_length = Read_length-int(double(Error_depth*3)*(float)Read_length/(float)Expect_depth);
    depth_half_length = Read_length - int(double(Expect_depth/2))*(double)Read_length/(double)Expect_depth;
    MAXDEPTH = 4 * Expect_depth;

    fprintf(stderr, "Assembly contig with library number %d and generating contigs with Output prefix %s\n", bwtCount, prefix);

    bwts = (BWT**)calloc(bwtCount, sizeof(BWT*));
    rev_bwts = (BWT**)calloc(bwtCount, sizeof(BWT*));
    int i;
    for(i=0; i< bwtCount; i++)
    {
        bwts[i] = idx2BWT[i]->bwt;
        rev_bwts[i] = idx2BWT[i]->rev_bwt;
    }
    fill_sparse_word();

    contig_assembly(prefix, 0);

    time(&now);
    fprintf(stderr, "Assembly DONE at %s\n", ctime(&now));

    fprintf(stderr, "Free index ... \n");
    free(bwts);
    free(rev_bwts);
    free_bwts(bwtCount);
    fprintf(stderr, "ALL DONE. THANK YOU!\n");

  return 0;
}



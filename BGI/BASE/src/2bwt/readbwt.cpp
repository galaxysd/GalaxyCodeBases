/*

   readbwt.cpp        Build index of short reads.

   This program builds index of short reads for use of assembler BASE and some other programs.

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

   Date   : 1st Jan 2014
   Author : Chi Man LIU
   Change : Developped this product to build BWT of short reads.
            This is the implementation of CX1 method.

   Date   : 1st Jan 2015
   Author : Binghang LIU
   Change : Add the implementation of CX1 method without GPU.
   			Add the loading of SOAPdenovo read libarary files, add qual_cutoff to set quality threshold.
   			Add a new output file for pair end reads with different insert sizes.
*/


#include <omp.h>
#include <pthread.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/sysinfo.h> // memory size.
#include <vector>
#include <algorithm>
#include <string>

// for sorting
#include <zlib.h>
#include "lv2_cpu_sort.h"
#include "lv2_gpu_sort.h"

using namespace std;

typedef unsigned int word;

#define NUM_CPU_THREADS 24
//#define THREADS_PER_BLOCK 256 // GPU

// do not modify below
#define BITS_PER_WORD 32 //used for quality.
#define BITS_PER_CHAR 2

#define BITS_PER_QUAL 1 //used for quality.

#define CHARS_PER_WORD 16
#define CHAR_MASK 0x3

#define PREFIX_LENGTH 8
#define SPECIAL_HANDLE_SUFFIX_LENGTH 3 // suffixes with length <= 3 will be handled specially
#define SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS 125 // 5^3
#define NUM_BUCKETS 390625
#define PREFIX_MASK 0xFFFF

char dna_map[256]; // map from ACGT to 0123
char dna_char[5] = {'A','C','G','T','$'};
static char tabs[2][1024];
word* g_packed_reads;
word* packed_quals;
word* packed_reads_p;

int words_per_qual;
int last_shift;
int last_shift_qual;
int number_cpu_threads = 24;

#define SENTINEL_INT 4 // '$'

//**********************************define structures start
struct global_data_t;
struct section_data_t {
    int section_id;
    struct global_data_t *parent;
    int64 section_start, section_end; // start and end read ID ('end' is exlusive)
    string bwt_buffer[SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS];
    int64 *offset_pos;
    int64 offset_diff_base;
};

struct suffix_fetch_data_t {
    int thread_id;
    struct global_data_t *parent;
    int bucket_start, bucket_end;
    int64 index_start;
};

struct global_data_t {
    int node_bucket_start, node_bucket_end; // buckets to be processed on this node
    int num_cpu_threads;
    int64 max_num_reads;
    int64 num_reads;
    word *packed_reads;
    uint64_t *cpu_sort_space;

    int read_len;
    int words_per_read;
    int words_per_suffix;
    int offset_num_bits; // number of bits to represent 0...read_len
    int offset_mask; // -(1 << offset_num_bits)

    int64 max_lv1_entries;
    int64 max_lv2_entries;

    int cur_special_handle_prefix;
    int64 num_lv2_to_output;
    bool need_ridt;

    int *buckets, *section_buckets[ NUM_CPU_THREADS ]; // buckets (aggregated for all read; for each section). TODO int64?
    pthread_t threads[ NUM_CPU_THREADS ];
    struct section_data_t locals[ NUM_CPU_THREADS ];
    struct suffix_fetch_data_t suffix_fetches[ NUM_CPU_THREADS - 1 ];
    int lv1_bucket_start, lv1_bucket_end; // 'end' is exclusive
    int lv2_bucket_start, lv2_bucket_end; // 'end' is exclusive
    int64 lv2_offset_start_index;
    int64 lv2_suffixes_row_width;
    int *lv1_offsets;
    word *sort_indexes; //permutation.
    word *sort_indexes_to_output;
    unsigned int *lv2_read_id;
    unsigned int *lv2_read_id_to_output;
    word *lv2_suffixes;
    word *lv2_suffixes_to_output;
    FILE *out_file;
    FILE *map_out_file;
    // for large diff
    vector<int64> v_large_diff;
    pthread_spinlock_t large_diff_lock;
};

typedef struct lib_info
{
    int sd;
    int qual;
    int max_ins;
    int min_ins;
    int avg_ins;
    int rd_len_cutoff;  //read length cutoff
    int reverse;
    int asm_flag;
    int map_len;
    int pair_num_cut;
    int rank;
    
    //whether last read is read1 in pair
    int paired;  // 0 -- single; 1 -- read1; 2 -- read2;

    //type1
    char ** a1_fname;
    char ** a2_fname;
    int num_a1_file;
    int num_a2_file;

    //type2
    char ** q1_fname;
    char ** q2_fname;
    int num_q1_file;
    int num_q2_file;

    //type3
    char ** p_fname;
    int num_p_file;  //fasta only

    //type4 &5
    char ** s_a_fname;
    int num_s_a_file;
    char ** s_q_fname;
    int num_s_q_file;
    
    char ** b_fname; //the name of the bam file
    int num_b_file; //the number of the bam file
} LIB_INFO;
LIB_INFO* lib_array;

//***************** define structure finished.

//***************** Initialize start.
//initial dna_map, here N is treated as G.
void init_dna_map() {
    for (int i='A';i<='Z';++i) dna_map[i] = 2; // G
    for (int i='a';i<='z';++i) dna_map[i] = 2; // G
    dna_map['A'] = dna_map['a'] = 0;
    dna_map['C'] = dna_map['c'] = 1;
    dna_map['G'] = dna_map['g'] = 2;
    dna_map['T'] = dna_map['t'] = 3;
    dna_map['N'] = dna_map['n'] = 2; // G
}

static boolean splitColumn ( char * line )
{
	int len = strlen ( line );
	int i = 0, j;
	int tabs_n = 0;

	while ( i < len )
	{
		if ( line[i] >= 32 && line[i] <= 126 && line[i] != '=' )
		{
			j = 0;
			while ( i < len && line[i] >= 32 && line[i] <= 126 && line[i] != '=' )
			{
				tabs[tabs_n][j++] = line[i];
				i++;
			}

			tabs[tabs_n][j] = '\0';
			tabs_n++;

			if ( tabs_n == 2 )
				return 1;
		}
		i++;
	}

	if ( tabs_n == 2 )
		return 1;
    return 0;
}

static int cmp_lib ( const void * a, const void * b )
{
	LIB_INFO * A, *B;
	A = ( LIB_INFO * ) a;
	B = ( LIB_INFO * ) b;

	if ( A->avg_ins > B->avg_ins )
		return 1;
	if ( A->avg_ins == B->avg_ins )
		return 0;
    return -1;
}

void scan_libInfo ( char * libfile, int *lib_num, int *read_len)
{
    char line[1024], ch;
    int i, j, index;
    boolean flag;
    boolean * pe;
    FILE* fp = fopen ( libfile, "r" );
    int num_libs = 0;
    int maxReadLen = 0;
    
	while ( fgets ( line, 1024, fp ) )
	{
		ch = line[5];
		line[5] = '\0';
		if ( strcmp ( line, "[LIB]" ) == 0 )
			num_libs++;

		if ( !num_libs )
		{
			line[5] = ch;
			flag = splitColumn ( line );

			if ( !flag )
				continue;

			if ( strcmp ( tabs[0], "max_rd_len" ) == 0 || strcmp ( tabs[0], "max_read_len" ) == 0)
				maxReadLen = atoi ( tabs[1] );
		}
	}

	if ( num_libs == 0 )
	{
		fprintf ( stderr, "Config file error: no [LIB] in file\n" );
		exit ( -1 );
	}
    *lib_num = num_libs;
    *read_len = maxReadLen;
    
	//count file numbers of each type
    lib_array = ( LIB_INFO * ) calloc (1, num_libs * sizeof ( LIB_INFO ) );
	pe = ( boolean * ) calloc (1, num_libs * sizeof ( boolean ) );

	for ( i = 0; i < num_libs; i++ )
	{
        lib_array[i].qual = 74;
        lib_array[i].sd = 0;
        lib_array[i].avg_ins = 0;
		lib_array[i].asm_flag = 3;
		lib_array[i].rank = 0;
		lib_array[i].pair_num_cut = 0;
		lib_array[i].rd_len_cutoff = 0;
		lib_array[i].map_len = 0;
		lib_array[i].num_s_a_file = 0;
		lib_array[i].num_s_q_file = 0;
		lib_array[i].num_p_file = 0;
		lib_array[i].num_a1_file = 0;
		lib_array[i].num_a2_file = 0;
		lib_array[i].num_q1_file = 0;
		lib_array[i].num_q2_file = 0;
		lib_array[i].num_b_file = 0;    //init
		pe[i] = false;
	}

	rewind ( fp );
	i = -1;

	while ( fgets ( line, 1024, fp ) )
	{
		ch = line[5];
		line[5] = '\0';

		if ( strcmp ( line, "[LIB]" ) == 0 )
		{
			i++;
			continue;
		}

		line[5] = ch;
		flag = splitColumn ( line );

		if ( !flag )
		{
			continue;
		}

		if ( strcmp ( tabs[0], "f1" ) == 0 )
		{
			lib_array[i].num_a1_file++;
			pe[i] = true;
		}
		else if ( strcmp ( tabs[0], "q1" ) == 0 )
		{
			lib_array[i].num_q1_file++;
			pe[i] = true;
		}
		else if ( strcmp ( tabs[0], "f2" ) == 0 )
		{
			lib_array[i].num_a2_file++;
			pe[i] = true;
		}
		else if ( strcmp ( tabs[0], "q2" ) == 0 )
		{
			lib_array[i].num_q2_file++;
			pe[i] = true;
		}
		else if ( strcmp ( tabs[0], "f" ) == 0 )
		{
			lib_array[i].num_s_a_file++;
		}
		else if ( strcmp ( tabs[0], "q" ) == 0 )
		{
			lib_array[i].num_s_q_file++;
		}
		else if ( strcmp ( tabs[0], "p" ) == 0 )
		{
			lib_array[i].num_p_file++;
			pe[i] = true;
		}
		else if ( strcmp ( tabs[0], "b" ) == 0 ) // the bam file
		{
			lib_array[i].num_b_file++;
			pe[i] = true;
		}
	}

	//allocate memory for filenames
	for ( i = 0; i < num_libs; i++ )
	{
		if ( lib_array[i].num_a2_file != lib_array[i].num_a1_file )
		{
			fprintf ( stderr, "Config file error: the number of mark \"f1\" is not the same as \"f2\"!\n" );
			exit ( -1 );
		}

		if ( lib_array[i].num_q2_file != lib_array[i].num_q1_file )
		{
			fprintf ( stderr, "Config file error: the number of mark \"q1\" is not the same as \"q2\"!\n" );
			exit ( -1 );
		}

		if ( lib_array[i].num_s_a_file )
		{
			lib_array[i].s_a_fname = ( char ** ) calloc (1, lib_array[i].num_s_a_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_s_a_file; j++ )
			{
				lib_array[i].s_a_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_s_q_file )
		{
			lib_array[i].s_q_fname = ( char ** ) calloc (1, lib_array[i].num_s_q_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_s_q_file; j++ )
			{
				lib_array[i].s_q_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_p_file )
		{
			lib_array[i].p_fname = ( char ** ) calloc (1, lib_array[i].num_p_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_p_file; j++ )
			{
				lib_array[i].p_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_a1_file )
		{
			lib_array[i].a1_fname = ( char ** ) calloc (1, lib_array[i].num_a1_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_a1_file; j++ )
			{
				lib_array[i].a1_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_a2_file )
		{
			lib_array[i].a2_fname = ( char ** ) calloc (1, lib_array[i].num_a2_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_a2_file; j++ )
			{
				lib_array[i].a2_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_q1_file )
		{
			lib_array[i].q1_fname = ( char ** ) calloc (1, lib_array[i].num_q1_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_q1_file; j++ )
			{
				lib_array[i].q1_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_q2_file )
		{
			lib_array[i].q2_fname = ( char ** ) calloc (1, lib_array[i].num_q2_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_q2_file; j++ )
			{
				lib_array[i].q2_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) );
			}
		}

		if ( lib_array[i].num_b_file )  //allot memory for bam file name
		{
			lib_array[i].b_fname = ( char ** ) calloc (1, lib_array[i].num_b_file * sizeof ( char * ) );

			for ( j = 0; j < lib_array[i].num_b_file; j++ )
				{ lib_array[i].b_fname[j] = ( char * ) calloc (1, 1024 * sizeof ( char ) ); }
		}
	}

	// get file names
	for ( i = 0; i < num_libs; i++ )
	{
		lib_array[i].num_s_a_file = 0;
		lib_array[i].num_s_q_file = 0;
		lib_array[i].num_p_file = 0;
		lib_array[i].num_a1_file = 0;
		lib_array[i].num_a2_file = 0;
		lib_array[i].num_q1_file = 0;
		lib_array[i].num_q2_file = 0;
		lib_array[i].num_b_file = 0;    //init
	}

	rewind ( fp );
	i = -1;

	while ( fgets ( line, 1024, fp ) )
	{
		ch = line[5];
		line[5] = '\0';

		if ( strcmp ( line, "[LIB]" ) == 0 )
		{
			i++;
			continue;
		}

		line[5] = ch;
		flag = splitColumn ( line );

		if ( !flag )
		{
			continue;
		}

		if ( strcmp ( tabs[0], "f1" ) == 0 )
		{
			index = lib_array[i].num_a1_file++;
			strcpy ( lib_array[i].a1_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "q1" ) == 0 )
		{
			index = lib_array[i].num_q1_file++;
			strcpy ( lib_array[i].q1_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "f2" ) == 0 )
		{
			index = lib_array[i].num_a2_file++;
			strcpy ( lib_array[i].a2_fname[index], tabs[1] );

			if ( strcmp ( lib_array[i].a2_fname[index], lib_array[i].a1_fname[index] ) == 0 )
			{
				fprintf ( stderr, "Config file error: f2 file is the same as f1 file\n" );
				fprintf ( stderr, "f1=%s\n", lib_array[i].a1_fname[index] );
				fprintf ( stderr, "f2=%s\n", lib_array[i].a2_fname[index] );
				exit ( -1 );
			}
		}
		else if ( strcmp ( tabs[0], "q2" ) == 0 )
		{
			index = lib_array[i].num_q2_file++;
			strcpy ( lib_array[i].q2_fname[index], tabs[1] );

			if ( strcmp ( lib_array[i].q2_fname[index], lib_array[i].q1_fname[index] ) == 0 )
			{
				fprintf ( stderr, "Config file error: q2 file is the same as q1 file\n" );
				fprintf ( stderr, "q1=%s\n", lib_array[i].q1_fname[index] );
				fprintf ( stderr, "q2=%s\n", lib_array[i].q2_fname[index] );
				exit ( -1 );
			}
		}
		else if ( strcmp ( tabs[0], "f" ) == 0 )
		{
			index = lib_array[i].num_s_a_file++;
			strcpy ( lib_array[i].s_a_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "q" ) == 0 )
		{
			index = lib_array[i].num_s_q_file++;
			strcpy ( lib_array[i].s_q_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "p" ) == 0 )
		{
			index = lib_array[i].num_p_file++;
			strcpy ( lib_array[i].p_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "b" ) == 0 )
		{
			//bam file
			index = lib_array[i].num_b_file++;
			strcpy ( lib_array[i].b_fname[index], tabs[1] );
		}
		else if ( strcmp ( tabs[0], "min_ins" ) == 0 )
		{
			lib_array[i].min_ins = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "max_ins" ) == 0 )
		{
			lib_array[i].max_ins = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "avg_ins" ) == 0 )
		{
			lib_array[i].avg_ins = atoi ( tabs[1] );
		}
        else if ( strcmp ( tabs[0], "sd" ) == 0)
        {
            lib_array[i].sd = atoi(tabs[1]);
        }
        else if( strcmp ( tabs[0], "qual_cutoff" ) == 0 )
        {
            lib_array[i].qual = atoi(tabs[1]);
        }
		else if ( strcmp ( tabs[0], "rd_len_cutoff" ) == 0 )
		{
			lib_array[i].rd_len_cutoff = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "reverse_seq" ) == 0 )
		{
			lib_array[i].reverse = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "asm_flags" ) == 0 )
		{
			lib_array[i].asm_flag = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "rank" ) == 0 )
		{
			lib_array[i].rank = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "pair_num_cutoff" ) == 0 )
		{
			lib_array[i].pair_num_cut = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "rd_len_cutoff" ) == 0 )
		{
			lib_array[i].rd_len_cutoff = atoi ( tabs[1] );
		}
		else if ( strcmp ( tabs[0], "map_len" ) == 0 )
		{
			lib_array[i].map_len = atoi ( tabs[1] );
		}
	}

	for ( i = 0; i < num_libs; i++ )
	{
		if ( pe[i] && lib_array[i].avg_ins == 0 )
		{
			fprintf ( stderr, "Config file error: PE reads need avg_ins in [LIB] %d\n", i + 1 );
			exit ( -1 );
		}
        
        //check sd.
        if( pe[i] && lib_array[i].sd == 0)
        {
            if(lib_array[i].max_ins > 0 || lib_array[i].min_ins > 0)
            {
                if( lib_array[i].max_ins > 0 && lib_array[i].min_ins > 0 )
                    lib_array[i].sd = lib_array[i].max_ins - lib_array[i].avg_ins > lib_array[i].min_ins - lib_array[i].avg_ins ? (lib_array[i].max_ins - lib_array[i].avg_ins)/3 : (lib_array[i].min_ins - lib_array[i].avg_ins)/3;
                else if( lib_array[i].max_ins > 0 )
                    lib_array[i].sd = (lib_array[i].max_ins - lib_array[i].avg_ins)/3;
                else
                    lib_array[i].sd = (lib_array[i].min_ins - lib_array[i].avg_ins)/3;
            }else{
                lib_array[i].sd = lib_array[i].avg_ins * 0.05;
            }
        }
        
	}

	fclose ( fp );
	qsort ( &lib_array[0], num_libs, sizeof ( LIB_INFO ), cmp_lib );
}

void free_libs (LIB_INFO* lib_array, int num_libs)
{
	if ( !lib_array )
	{
		return;
	}

	int i, j;
	fprintf ( stderr, "LIB(s) information:\n" );

	for ( i = 0; i < num_libs; i++ )
	{
		fprintf ( stderr, " [LIB] %d, avg_ins %d, reverse %d.\n", i, lib_array[i].avg_ins, lib_array[i].reverse );

		if ( lib_array[i].num_s_a_file )
		{
			//printf("%d single fasta files\n",lib_array[i].num_s_a_file);
			for ( j = 0; j < lib_array[i].num_s_a_file; j++ )
			{
				free ( ( void * ) lib_array[i].s_a_fname[j] );
			}

			free ( ( void * ) lib_array[i].s_a_fname );
		}

		if ( lib_array[i].num_s_q_file )
		{
			//printf("%d single fastq files\n",lib_array[i].num_s_q_file);
			for ( j = 0; j < lib_array[i].num_s_q_file; j++ )
			{
				free ( ( void * ) lib_array[i].s_q_fname[j] );
			}

			free ( ( void * ) lib_array[i].s_q_fname );
		}

		if ( lib_array[i].num_p_file )
		{
			//printf("%d paired fasta files\n",lib_array[i].num_p_file);
			for ( j = 0; j < lib_array[i].num_p_file; j++ )
			{
				free ( ( void * ) lib_array[i].p_fname[j] );
			}

			free ( ( void * ) lib_array[i].p_fname );
		}

		if ( lib_array[i].num_a1_file )
		{
			//printf("%d read1 fasta files\n",lib_array[i].num_a1_file);
			for ( j = 0; j < lib_array[i].num_a1_file; j++ )
			{
				free ( ( void * ) lib_array[i].a1_fname[j] );
			}

			free ( ( void * ) lib_array[i].a1_fname );
		}

		if ( lib_array[i].num_a2_file )
		{
			//printf("%d read2 fasta files\n",lib_array[i].num_a2_file);
			for ( j = 0; j < lib_array[i].num_a2_file; j++ )
			{
				free ( ( void * ) lib_array[i].a2_fname[j] );
			}

			free ( ( void * ) lib_array[i].a2_fname );
		}

		if ( lib_array[i].num_q1_file )
		{
			//printf("%d read1 fastq files\n",lib_array[i].num_q1_file);
			for ( j = 0; j < lib_array[i].num_q1_file; j++ )
			{
				free ( ( void * ) lib_array[i].q1_fname[j] );
			}

			free ( ( void * ) lib_array[i].q1_fname );
		}

		if ( lib_array[i].num_q2_file )
		{
			//printf("%d read2 fastq files\n",lib_array[i].num_q2_file);
			for ( j = 0; j < lib_array[i].num_q2_file; j++ )
			{
				free ( ( void * ) lib_array[i].q2_fname[j] );
			}

			free ( ( void * ) lib_array[i].q2_fname );
		}

		if ( lib_array[i].num_b_file )
		{
			//free the bam file name
			//printf("%d bam files\n",lib_array[i].num_b_file);
			for ( j = 0; j < lib_array[i].num_b_file; j++ )
				{ free ( ( void * ) lib_array[i].b_fname[j] ); }

			free ( ( void * ) lib_array[i].b_fname );
		}
	}

	num_libs = 0;
	free ( ( void * ) lib_array );
}
//***************** Initialize finished.

//***************** for debugging only start 
void debug_prefix( FILE* file, word prefix ) {
    for (int i = PREFIX_LENGTH - 1; i >= 0; --i) {
        int c = prefix % 5;
        prefix /= 5;
        fprintf( file, "%c", dna_char[c] );
    }
    fprintf( file, "^r" );
}

void log(const char* format, ...) {  //add const to char* can remove warning.
    va_list args;
    va_start( args, format );
    vfprintf( stderr, format, args );
    va_end( args );
    fflush( stderr );
}

// for debugging only
void debug_word( FILE* file, word w ) {
    for (int i = CHARS_PER_WORD - 1; i >= 0; --i) {
        int c = (w>>(BITS_PER_CHAR*i)) & CHAR_MASK;
        fprintf( file, "%c", dna_char[c] );
    }
}

// for debugging only
void dump_special_prefix(int key) {
    for (int i = 0; i < SPECIAL_HANDLE_SUFFIX_LENGTH; ++i) {
        log("%c", dna_char[key / (SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS / 5)]);
        key = key * 5 % SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS;
    }
    log("\n");
}

// for debug only
void dump_read(word *read_p, int words_per_read, int read_len) {
    for (int i = 0; i < words_per_read; ++i) {
        for (int j = 0; j < CHARS_PER_WORD; ++j) {
            log("%c", dna_char[(read_p[i] >> (CHARS_PER_WORD - 1 - j) * BITS_PER_CHAR) & 3]);
            if (i * CHARS_PER_WORD + j == read_len - 1) {
                log("_");
            } 
        }
    }
    log("\n");
}

void debug_suffix( word* start, struct global_data_t *_data ) {
    struct global_data_t &globals = *_data;
    for (int j=0; j<globals.words_per_read; ++j) {
        word w = *(start + j * globals.lv2_suffixes_row_width);
        for (int i = CHARS_PER_WORD - 1; i >= 0; --i) {
            int c = (w>>(BITS_PER_CHAR*i)) & CHAR_MASK;
            fprintf(stderr, "%c", dna_char[c]);
        }
    }
}
//***************** for debugging only end 

inline void* MallocAndCheck(size_t size_in_byte,
                            const char *malloc_from_which_file = __FILE__,
                            int malloc_from_which_line = __LINE__) {
    void *ptr = malloc(size_in_byte);
    if (ptr == NULL && size_in_byte != 0) {
        log( "[ERROR] Ran out of memory while applying %llubytes\n", (unsigned long long)size_in_byte);
        log( "In file: %s, line %d\n", malloc_from_which_file, malloc_from_which_line);
        log( "There may be errors as follows:\n");
        log( "1) Not enough memory.\n");
        log( "2) The ARRAY may be overrode.\n");
        log( "3) The wild pointers.\n");
        exit(-1);
    }

    return ptr;
}

inline bool is_valid_special_prefix(int prefix_key) {
    // TODO hard
    if (prefix_key % 5 != 4) { return false; }
    for (int i = 0, start = 0; i < SPECIAL_HANDLE_SUFFIX_LENGTH; ++i) {
        if (prefix_key % 5 != 4) {
            start = 1;
        } else if (start) {
            return false;
        }
        prefix_key /= 5;
    }
    return prefix_key == 0;
}

void* preprocess_fill_buckets( void *_data ) {
    struct section_data_t &locals = *( (struct section_data_t*)_data );
    struct global_data_t &globals = *(locals.parent);
    int section_id = locals.section_id;
    int *buckets = globals.section_buckets[ section_id ];

    int read_len = globals.read_len;
    int words_per_read = globals.words_per_read;
    int w_shift = BITS_PER_WORD - BITS_PER_CHAR;

    for ( int64 r = locals.section_start; r < locals.section_end; ++r )  {
        word *w_p = globals.packed_reads + ( words_per_read * r );
        word key = 0;
        word w = *(w_p++);
        for ( int i = 0; i < PREFIX_LENGTH - 1; ++i ) {
            key = key * 5 + ( w >> w_shift ); // replaced by this
            w <<= BITS_PER_CHAR;
        }
        for ( int i = PREFIX_LENGTH - 1; i < read_len; ++i ) {
            if ( i % CHARS_PER_WORD == 0 ) // TODO optimize
                w = *(w_p++);
            key = ( key * 5 + ( w >> w_shift ) ) % NUM_BUCKETS; // replaced by this
            w <<= BITS_PER_CHAR;
            buckets[key]++;
        }
        int special_key = (key * 5 + SENTINEL_INT) % (SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS * 5);

        for ( int i = PREFIX_LENGTH - 1; i >= SPECIAL_HANDLE_SUFFIX_LENGTH; --i ) {
            key = ( key * 5 + SENTINEL_INT ) % NUM_BUCKETS; // replaced by this
            buckets[key]++;
        }

        // special handle buckets
        for ( int i = 0; i < SPECIAL_HANDLE_SUFFIX_LENGTH; ++i ) {
            int index = special_key % SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS;
            assert(is_valid_special_prefix(index));
            assert(special_key / SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS != 4);
            locals.bwt_buffer[index].push_back(dna_char[ special_key / SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS ]);
            special_key = (special_key * 5 + SENTINEL_INT) % (SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS * 5);
        }
    }

    return NULL;
}

void* lv1_fetch_offsets( void *_data) {
    struct section_data_t &locals = *( (struct section_data_t*)_data );
    struct global_data_t &globals = *(locals.parent);
  
    const int bucket_start = globals.lv1_bucket_start;
    const int bucket_end = globals.lv1_bucket_end;

    const int section_start = locals.section_start;
    const int section_end = locals.section_end;

    int read_len = globals.read_len;
    int words_per_read = globals.words_per_read;
    int w_shift = BITS_PER_WORD - BITS_PER_CHAR;
    int64 *offset_pos = locals.offset_pos;

    int id_shift = globals.offset_num_bits;
    int64 *recent_offsets = (int64*) malloc( NUM_BUCKETS * sizeof( int64 ) ); // for calculating differences
    int64 recent_base = ( locals.section_start << id_shift ); // offset of first read first char
    for ( int b = bucket_start; b < bucket_end; ++b ) {
        recent_offsets[b] = recent_base;
    }
    locals.offset_diff_base = recent_base;

    // scan all reads and record offsets
    for ( int64 r = section_start; r < section_end; ++r )  {
        word *w_p = globals.packed_reads + ( words_per_read * r );
        word key = 0;
        word w = *(w_p++);
        //        int prev_char = SENTINEL_INT;
        for ( int i = 0; i < PREFIX_LENGTH - 1; ++i ) {
            //            key = ((key<<BITS_PER_CHAR) | (w>>w_shift));
            key = key * 5 + ( w >> w_shift ); // replaced by this
            w <<= BITS_PER_CHAR;
        }
        for ( int i = PREFIX_LENGTH - 1; i < read_len; ++i ) {
            if ( i % CHARS_PER_WORD == 0 ) // TODO optimize
                w = *(w_p++);
            //            key = (((key<<BITS_PER_CHAR)&PREFIX_MASK)|(w>>w_shift));
            key = ( key * 5 + ( w >> w_shift ) ) % NUM_BUCKETS; // replaced by this
            w <<= BITS_PER_CHAR;
            // check if key within range
            // TODO 0x80000000 hardcode
            // this condition is equivalent to bucket_start <= key < bucket_end
            if ( (( key - bucket_start ) ^ ( key - bucket_end )) & 0x80000000 ) {
                // read ID = r, within-read offset = i - PREFIX_LENGTH + 1
                int64 true_offset = ( r << id_shift ) | ( i - PREFIX_LENGTH + 1 );
                int64 diff_offset = true_offset - recent_offsets[ key ];

                if ( diff_offset < 0 ) { // TODO remove this assertion
                    log("ERROR: offset diff negative (%lld %lld)! ABORT.\n", recent_offsets[key], true_offset);
                    exit(1);
                }
                if ( diff_offset > 2147483647 ) {
                    pthread_spin_lock(&globals.large_diff_lock);
                    diff_offset = -1 - (int64)globals.v_large_diff.size();
                    globals.v_large_diff.push_back(true_offset);
                    pthread_spin_unlock(&globals.large_diff_lock);
                }
                globals.lv1_offsets[ offset_pos[key]++ ] = diff_offset;
                recent_offsets[ key ] = true_offset;
            }
        }
        for ( int i = PREFIX_LENGTH - 1; i >= SPECIAL_HANDLE_SUFFIX_LENGTH; --i ) {
            key = ( key * 5 + SENTINEL_INT ) % NUM_BUCKETS; // replaced by this            
            // check if key within range
            // TODO 0x80000000 hardcode
            // this condition is equivalent to bucket_start <= key < bucket_end
            if ( (( key - bucket_start ) ^ ( key - bucket_end )) & 0x80000000 ) {
                // read ID = r, within-read offset = read_len - i
                int64 true_offset = ( r << id_shift ) | ( read_len - i );
                int64 diff_offset = true_offset - recent_offsets[ key ];

                if ( diff_offset < 0 ) { // TODO remove this assertion
                    log("ERROR: offset diff negative (%lld %lld)! ABORT.\n", recent_offsets[key], true_offset);
                    exit(1);
                }
                if ( diff_offset > 2147483647 ) {
                    pthread_spin_lock(&globals.large_diff_lock);
                    diff_offset = -1 - (int64)globals.v_large_diff.size();
                    globals.v_large_diff.push_back(true_offset);
                    pthread_spin_unlock(&globals.large_diff_lock);
                }
                globals.lv1_offsets[ offset_pos[key]++ ] = diff_offset;
                recent_offsets[ key ] = true_offset;
            }            
        }

    }

    free( recent_offsets );

    return NULL;
}

void* lv2_fetch_suffixes( void *_data ){
    struct suffix_fetch_data_t &fetch = *( (struct suffix_fetch_data_t*)_data );
    struct global_data_t &globals = *(fetch.parent);

    int read_len = globals.read_len;
    int offset_num_bits = 1;
    while ((1 << offset_num_bits) - 1 < read_len) {
        ++offset_num_bits;
    }
    int words_per_read = globals.words_per_read;
    int words_per_suffix = globals.words_per_suffix;
    const int bucket_start = fetch.bucket_start;
    const int bucket_end = fetch.bucket_end;
    word *packed_reads = globals.packed_reads;

    int id_shift = globals.offset_num_bits; // offsets are in the form [read_id][read_offset]; id_shift = # bits for read_offset
    int read_offset_mask = globals.offset_mask;

    word *suffix_p = globals.lv2_suffixes + fetch.index_start;
    int num_suffixes = fetch.index_start;
    
    int *offset_p = globals.lv1_offsets + globals.lv2_offset_start_index + fetch.index_start;
    for ( int b = bucket_start; b < bucket_end; ++b ) {

        for ( int section = 0; section < globals.num_cpu_threads; ++section ) {

            int64 offset = globals.locals[section].offset_diff_base; // reset per local bucket
            int64 num_in_section_bucket = globals.section_buckets[section][b];
            while ( num_in_section_bucket-- ) {
                if (*offset_p < 0) {
                    offset = globals.v_large_diff[-(*offset_p++) - 1];
                } else {
                    offset += *(offset_p++); // remember we are using difference encoding
                }
                // offset in the form [read_id][read_offset]
                int64 read_id = offset >> id_shift;
                int read_offset = offset & read_offset_mask;
                word *read_p = packed_reads + ( words_per_read * read_id );

                // find previous character (the char to be stored in BWT)
                int prev_char;
                if ( !read_offset ) { // special case
                    prev_char = SENTINEL_INT;
                    globals.lv2_read_id[ num_suffixes ] = read_id;
                } else {
                    int prev_which_word = ( read_offset - 1 ) / CHARS_PER_WORD;
                    int prev_word_offset = ( read_offset - 1 ) % CHARS_PER_WORD;
                    word w = *(read_p + prev_which_word);
                    prev_char = ( w >> ( BITS_PER_WORD - BITS_PER_CHAR - prev_word_offset * 2 ) ) & CHAR_MASK;
                }

                // copy words of the suffix to the suffix pool
                int which_word = read_offset / CHARS_PER_WORD;
                int word_offset = read_offset % CHARS_PER_WORD;
                word *dest = suffix_p;
                word *src = read_p + which_word;
                if ( !word_offset ) { // special case (word aligned), easy
                    while ( which_word < words_per_read ) {
                        *dest = *src; // write out
                        dest += globals.lv2_suffixes_row_width;
                        src++;
                        which_word++;
                    }
                } else { // not word-aligned
                    int bit_shift = read_offset * 2;
                    word s = *src;
                    word d = s << bit_shift;
                    which_word++;
                    while ( which_word < words_per_read ) {
                        s = *(++src);
                        d |= s >> (BITS_PER_WORD - bit_shift);
                        *dest = d; // write out
                        dest += globals.lv2_suffixes_row_width;
                        d = s << bit_shift;
                        which_word++;
                    }
                    // write last word
                    *dest = d | ( 0xFFFFFFFF >> ( BITS_PER_WORD - bit_shift ) ); // bit padding by '1's. TODO hard
                }
                *(suffix_p + (globals.lv2_suffixes_row_width * (words_per_suffix - 1) ) ) &= ~globals.offset_mask; // TODO hard
                *(suffix_p + (globals.lv2_suffixes_row_width * (words_per_suffix - 1) ) ) |= read_offset;
                suffix_p++;
                globals.sort_indexes[ num_suffixes ] = ( prev_char << 29 ) | num_suffixes; // TODO hard
                num_suffixes++;
            }
        }
    }
        
    return NULL;
}

int get_base5_prefix(word item) {
    item >>= (CHARS_PER_WORD - SPECIAL_HANDLE_SUFFIX_LENGTH) * BITS_PER_CHAR;
    static int base4_mask = (1 << SPECIAL_HANDLE_SUFFIX_LENGTH * BITS_PER_CHAR) - 1;
    int ret = 0;
    for (int i = 0; i < SPECIAL_HANDLE_SUFFIX_LENGTH; ++i) {
        ret = ret * 5 + (item >> (SPECIAL_HANDLE_SUFFIX_LENGTH - 1) * BITS_PER_CHAR);
        item = item * 4 & base4_mask;
    }
    return ret;
}

inline unsigned int mirror(unsigned int v) {
    // swap consecutive pairs
    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    // swap nibbles ... 
    v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    // swap bytes
    v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    // swap 2-byte long pairs
    v = ( v >> 16             ) | ( v               << 16);
    return v;
}

inline void reverse_one_read(word *read_p, int words_per_read, int last_shift) {
    // shift all words
    if (last_shift > 0) {
        for (int i = words_per_read - 1; i > 0; --i) {
            read_p[i] >>= last_shift;
            read_p[i] |= read_p[i - 1] << (BITS_PER_WORD - last_shift);
        }
        read_p[0] >>= last_shift;
        read_p[0] |= -(1 << (BITS_PER_WORD - last_shift));
    }
    // mirror all words
    for (int i = 0; i < words_per_read; ++i) {
        read_p[i] = mirror(read_p[i]);
    }
}

void reverse_all_reads(global_data_t &globals) {
    omp_set_num_threads(globals.num_cpu_threads);
    int last_shift = globals.read_len * BITS_PER_CHAR % BITS_PER_WORD;
    if (last_shift != 0) last_shift = BITS_PER_WORD - last_shift;
    
#pragma omp for
    for (int64 i = 0; i < globals.num_reads; ++i) {
        reverse_one_read(globals.packed_reads + globals.words_per_read * i, globals.words_per_read, last_shift);
    }

    int64 mid = globals.num_reads / 2;

#pragma omp for
    for (int64 i = 0; i < mid; ++i) {
        int64 i2 = globals.num_reads - 1 - i;
        for (int j = 0; j < globals.words_per_read; ++j) {
            swap(globals.packed_reads[i * globals.words_per_read + j], globals.packed_reads[(i2 + 1) * globals.words_per_read - 1 - j]);
        }
    }

    if (globals.num_reads % 2 == 1) {
        for (int j = 0; j < globals.words_per_read / 2; ++j) {
            swap(globals.packed_reads[mid * globals.words_per_read + j], globals.packed_reads[(mid + 1) * globals.words_per_read - 1 - j]);
        }
    }
}

void init_global_data(global_data_t &globals, int read_len, int64 max_host_mem, int number_cpu_threads) {
    globals.read_len = read_len;
    globals.words_per_read = ( read_len + ( CHARS_PER_WORD - 1 ) ) / CHARS_PER_WORD;
    globals.num_reads = 0;
    globals.max_num_reads = max_host_mem * 0.8 / (globals.words_per_read * sizeof(word));
    globals.num_cpu_threads = number_cpu_threads;
    globals.max_lv2_entries = 0;

    globals.offset_num_bits = 1;
    while ((1 << globals.offset_num_bits) - 1 < read_len) {
        ++globals.offset_num_bits;
    }
    globals.offset_mask = (1 << globals.offset_num_bits) - 1;
    globals.words_per_suffix = (read_len * BITS_PER_CHAR + globals.offset_num_bits + BITS_PER_WORD - 1) / BITS_PER_WORD;

    log("Max host mem: %lld, max number of reads can be loaded: %lld\n", max_host_mem, globals.max_num_reads);
    log("Words per read: %d, words per suffix: %d\n", globals.words_per_read, globals.words_per_suffix);
}

void set_global_bucket(global_data_t &globals)
{
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        globals.locals[t].parent = &globals;
        globals.locals[t].section_id = t;
    }
    for ( int t = 0; t < globals.num_cpu_threads - 1; ++t ) {
        globals.suffix_fetches[t].parent = &globals;
        globals.suffix_fetches[t].thread_id = t;
    }

    globals.buckets = NULL;
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        globals.section_buckets[t] = NULL;
        globals.locals[t].offset_pos = NULL;
    }

    // realloc
    int64 mem_packed_reads = ( globals.num_reads + 32 ) * globals.words_per_read * sizeof(word); // +32 for safety
    globals.packed_reads = (word*) realloc( globals.packed_reads, mem_packed_reads );
    if ( ! globals.packed_reads ) {
        fprintf( stderr, "realloc error, quit\n" );
        exit(1);
    }
    log("Mem allocated for packed reads: %lld bytes\n", mem_packed_reads );

    { // distribute reads among threads
        int64 each = globals.num_reads / globals.num_cpu_threads;
        for ( int t = 0; t < globals.num_cpu_threads - 1; ++t ) {
            globals.locals[t].section_start = t * each;
            globals.locals[t].section_end = (t + 1) * each;
        }
        globals.locals[ globals.num_cpu_threads - 1 ].section_start = (globals.num_cpu_threads-1) * each;
        globals.locals[ globals.num_cpu_threads - 1 ].section_end = globals.num_reads;
    }
}


void set_global_after_reading_initial(global_data_t &globals, int64 max_host_mem, int64 free_gpu_mem) {
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        globals.locals[t].parent = &globals;
        globals.locals[t].section_id = t;
    }
    for ( int t = 0; t < globals.num_cpu_threads - 1; ++t ) {
        globals.suffix_fetches[t].parent = &globals;
        globals.suffix_fetches[t].thread_id = t;
    }

    globals.buckets = NULL;
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        globals.section_buckets[t] = NULL;
        globals.locals[t].offset_pos = NULL;
    }

    // realloc
    int64 mem_packed_reads = ( globals.num_reads + 32 ) * globals.words_per_read * sizeof(word); // +32 for safety
    globals.packed_reads = (word*) realloc( globals.packed_reads, mem_packed_reads );
    if ( ! globals.packed_reads ) {
        fprintf( stderr, "realloc error, quit\n" );
        exit(1);
    }
    log("Mem allocated for packed reads: %lld bytes\n", mem_packed_reads );

    { // distribute reads among threads
        int64 each = globals.num_reads / globals.num_cpu_threads;
        for ( int t = 0; t < globals.num_cpu_threads - 1; ++t ) {
            globals.locals[t].section_start = t * each;
            globals.locals[t].section_end = (t + 1) * each;
        }
        globals.locals[ globals.num_cpu_threads - 1 ].section_start = (globals.num_cpu_threads-1) * each;
        globals.locals[ globals.num_cpu_threads - 1 ].section_end = globals.num_reads;
    }

    // preprocessing - compute memory limits

    free_gpu_mem -= 1073741824; // remain 1G for b40c
    const int GPU_BYTES_PER_ENTRY = 16; // key, value, and x2 for radix sort buffer
    int64 max_lv2_entries = free_gpu_mem / GPU_BYTES_PER_ENTRY;

    if (max_lv2_entries >= (1 << 29)) {
        max_lv2_entries = (1 << 29) - 1; // TODO hard
    }
    globals.max_lv2_entries = max_lv2_entries;
    log( "Free GPU MEM = %lld\n", free_gpu_mem);

    int64 lv2_bytes_per_entry = ( globals.words_per_suffix + 1 + 1) * sizeof(word); // read_id_table, permutation
    int64 mem_lv2 = lv2_bytes_per_entry * max_lv2_entries * 2; // *2 for double buffering

    log( "max host mem: %lld, mem_packed_reads %lld, mem_lv2 %lld, global.num_reads %lld", max_host_mem, mem_packed_reads, mem_lv2, globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2);
    int64 avail_host_mem = max_host_mem - mem_packed_reads - mem_lv2 - globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2;
    
    if(avail_host_mem < 1024*1024)
    {
        log( "Please set max_host_mem larger than: %lld\n", mem_packed_reads + mem_lv2 + globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2 + 1024*1024);
        exit(2);
    }
    const int LV1_BYTES_PER_ENTRY = 4; // 32-bit adjacent-difference representation
    int64 max_lv1_entries = avail_host_mem / LV1_BYTES_PER_ENTRY;
    globals.max_lv1_entries = max_lv1_entries;

    log( "Avail host MEM = %lld\n", avail_host_mem );
    log( "Max # of lv.1 entries = %lld\n", max_lv1_entries );
    log( "Max # of lv.2 entries = %lld\n", max_lv2_entries );
    log( "\n" );

#ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_entries, __FILE__, __LINE__); // as CPU memory is used to simulate GPU
    assert(globals.cpu_sort_space != NULL);
#endif

    // allocate memory for lv.1 and lv.2 arrays
    globals.lv1_offsets = (int*) malloc( max_lv1_entries * sizeof(int) );
    globals.lv2_suffixes = (word*) malloc( max_lv2_entries * globals.words_per_suffix * sizeof(word) );
    globals.lv2_suffixes_row_width = max_lv2_entries;
    globals.lv2_read_id = (unsigned int*) malloc( max_lv2_entries * sizeof(unsigned int) );
    globals.sort_indexes = (word*) malloc( max_lv2_entries * sizeof(word) );

    globals.lv2_suffixes_to_output = (word*) malloc( max_lv2_entries * globals.words_per_suffix * sizeof(word) );
    globals.lv2_read_id_to_output = (unsigned int*) malloc( max_lv2_entries * sizeof(unsigned int) );
    globals.sort_indexes_to_output = (word*) malloc( max_lv2_entries * sizeof(word) );

    pthread_spin_init(&globals.large_diff_lock, 0);
}

void set_global_after_lv0(global_data_t &globals, int64 max_host_mem, int64 free_gpu_mem, uint64_t max_bucket_size)
{
    int64 max_lv2_entries;

#ifdef DISABLE_GPU
    if(max_bucket_size*1.2 < 2*1024*1024)
        max_lv2_entries = 2*1024*1024;
    else
        max_lv2_entries = max_bucket_size * 1.2;

    globals.max_lv2_entries = max_lv2_entries;
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.max_lv2_entries, __FILE__, __LINE__); // as CPU memory is used to simulate GPU
    assert(globals.cpu_sort_space != NULL);
#else
    const int GPU_BYTES_PER_ENTRY = 16; // key, value, and x2 for radix sort buffer
    max_lv2_entries = (free_gpu_mem - 1073741824) / GPU_BYTES_PER_ENTRY;
    if (max_lv2_entries >= (1 << 29)) {
        max_lv2_entries = (1 << 29) - 1;
    }
    if(max_lv2_entries < (int64_t)max_bucket_size)
    {
        log ("Too many items for GPU sorting %lld > %lld! Please try larger GPU memory or use CPU version instead.", (int64_t)max_bucket_size, max_lv2_entries);
        exit(1);
    }
    globals.max_lv2_entries = max_lv2_entries;
#endif
    int64 lv2_bytes_per_entry = ( globals.words_per_suffix + 1 + 1) * sizeof(word); // read_id_table, permutation
    int64 mem_lv2 = lv2_bytes_per_entry * max_lv2_entries * 2; // *2 for double buffering

    int64 mem_packed_reads = ( globals.num_reads + 32 ) * globals.words_per_read * sizeof(word); // +32 for safety

    log( "max host mem: %lld, mem_packed_reads %lld, mem_lv2 %lld, global.num_reads %lld", max_host_mem, mem_packed_reads, mem_lv2, globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2);
    int64 avail_host_mem = max_host_mem - mem_packed_reads - mem_lv2 - globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2;

    if(avail_host_mem < 1024*1024)
    {
        log( "Please set max_host_mem larger than: %lld\n", mem_packed_reads + mem_lv2 + globals.num_reads * sizeof(char) * SPECIAL_HANDLE_SUFFIX_LENGTH * 2 + 1024*1024);
        exit(2);
    }

    const int LV1_BYTES_PER_ENTRY = 4; // 32-bit adjacent-difference representation
    int64 max_lv1_entries = avail_host_mem / LV1_BYTES_PER_ENTRY;
    globals.max_lv1_entries = max_lv1_entries;

    log( "Avail host MEM = %lld\n", avail_host_mem );
    log( "Max # of lv.1 entries = %lld\n", max_lv1_entries );
    log( "Max # of lv.2 entries = %lld\n", max_lv2_entries );
    log( "\n" );

    globals.lv1_offsets = (int*) malloc( max_lv1_entries * sizeof(int) );
    globals.lv2_suffixes = (word*) malloc( max_lv2_entries * globals.words_per_suffix * sizeof(word) );
    globals.lv2_suffixes_row_width = max_lv2_entries;
    globals.lv2_read_id = (unsigned int*) malloc( max_lv2_entries * sizeof(unsigned int) );
    globals.sort_indexes = (word*) malloc( max_lv2_entries * sizeof(word) );

    globals.lv2_suffixes_to_output = (word*) malloc( max_lv2_entries * globals.words_per_suffix * sizeof(word) );
    globals.lv2_read_id_to_output = (unsigned int*) malloc( max_lv2_entries * sizeof(unsigned int) );
    globals.sort_indexes_to_output = (word*) malloc( max_lv2_entries * sizeof(word) );

    pthread_spin_init(&globals.large_diff_lock, 0);
}

uint64_t lv0_fill_buckets(global_data_t &globals) {
    // preprocessing - fill buckets
    if (globals.buckets == NULL) 
        globals.buckets = (int*) malloc( NUM_BUCKETS * sizeof(int) );
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        if (globals.section_buckets[t] == NULL) 
            globals.section_buckets[t] = (int*) malloc( NUM_BUCKETS * sizeof(int) );
        memset(globals.section_buckets[t], 0, sizeof(int) * NUM_BUCKETS);
    }
    log( "Preprocessing: filling buckets... " );
    for ( int t = 1; t < globals.num_cpu_threads; ++t ) {
        pthread_create( &globals.threads[t], NULL, preprocess_fill_buckets, &globals.locals[t] );
    }
    preprocess_fill_buckets( &globals.locals[0] );
    for ( int t = 1; t < globals.num_cpu_threads; ++t ) {
        pthread_join( globals.threads[t], NULL );
    }
    log( " Done\n" );

    // sum up local buckets
    memset( globals.buckets, 0, NUM_BUCKETS * sizeof(int) );
    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        for ( int b = 0; b < NUM_BUCKETS; ++b )
            globals.buckets[b] += globals.section_buckets[t][b];
        if (globals.locals[t].offset_pos == NULL)
            globals.locals[t].offset_pos = (int64*) malloc( NUM_BUCKETS * sizeof( int64 ) );
        memset(globals.locals[t].offset_pos, 0, sizeof(int64) * NUM_BUCKETS);
    }

    int64_t max_buckets=0;
    for (int i = 0; i < NUM_BUCKETS; ++i) 
    {
        if(globals.max_lv2_entries == 0)
        {
            if(globals.buckets[i] > max_buckets)
                max_buckets = globals.buckets[i];
        }else{
            if (globals.buckets[i] > globals.max_lv2_entries) {
                fprintf(stderr, "Too many items for GPU sorting! Bucket: %d, number of items: %d\n", i, globals.buckets[i]);
                exit(1);
            }
        }
    }
    if(globals.max_lv2_entries == 0)
        log ("Maximum Bucket size: %lld\n", max_buckets);

    return max_buckets;
}

void* lv2_output_thread(void *_data) {
    global_data_t &globals = *((global_data_t*)_data);
    log( "Writing output... \n" );
    for (int i = 0; i < globals.num_lv2_to_output; ++i) {
        int c = globals.sort_indexes_to_output[i] >> 29; // TODO hard
        int cur_prefix = get_base5_prefix( globals.lv2_suffixes_to_output[globals.sort_indexes_to_output[i] & ((1 << 29) - 1)]);
        assert(cur_prefix < SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS);
        while (cur_prefix > globals.cur_special_handle_prefix) {
            dump_special_prefix(globals.cur_special_handle_prefix);
            for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
                fwrite( globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].c_str(), sizeof(char),
                        globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].size(), globals.out_file );
                globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].clear();
            }

            ++globals.cur_special_handle_prefix;
            while (globals.cur_special_handle_prefix < SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS) {
                if (is_valid_special_prefix(globals.cur_special_handle_prefix)) {
                    break;
                }
                ++globals.cur_special_handle_prefix;
            }
        }
        fprintf(globals.out_file, "%c", dna_char[ c ] );
        if ( c == SENTINEL_INT && globals.need_ridt ) {
            // TODO offsets are accumulated....
            unsigned int id = globals.lv2_read_id_to_output[ globals.sort_indexes_to_output[i] & 0x1FFFFFFF ]; // TODO hard
            fwrite( &id, sizeof(unsigned int), 1, globals.map_out_file );
        }
    }
    return NULL;
}

void build_bwt(global_data_t &globals, const char *output_prefix, bool need_ridt = false) {

    for (int64 i = globals.num_reads - 10; i < globals.num_reads; ++i) {
        dump_read(globals.packed_reads + i * globals.words_per_read, globals.words_per_read, globals.read_len);
    }

    pthread_t output_thread;
    bool output_started = false;

    if(need_ridt)
        globals.out_file = fopen64( (string(output_prefix) + ".f.bwt.ascii").c_str(), "w" );
    else
        globals.out_file = fopen64( (string(output_prefix) + ".r.bwt.ascii").c_str(), "w" );
    globals.need_ridt = need_ridt;

    if (need_ridt) {
        globals.map_out_file = fopen64( (string(output_prefix) + ".ridt").c_str(), "wb" );
        fwrite( &globals.num_reads, sizeof(unsigned int), 1, globals.map_out_file );
    }
    // level 1 loop - pick consecutive buckets
    int lv1_round = 0;
    globals.cur_special_handle_prefix = 0;
    while (!is_valid_special_prefix(globals.cur_special_handle_prefix)) {
        ++globals.cur_special_handle_prefix;
    }
    globals.lv1_bucket_start = 0;
    while ( globals.lv1_bucket_start < NUM_BUCKETS ) {
        lv1_round++;
        globals.v_large_diff.clear();
        int64 num_lv1 = globals.buckets[ globals.lv1_bucket_start ];
        int temp_end = globals.lv1_bucket_start + 1; // remember that 'end' is exclusive
        while ( temp_end < NUM_BUCKETS && num_lv1 + globals.buckets[ temp_end ] <= globals.max_lv1_entries ) {
            num_lv1 += globals.buckets[ temp_end++ ];
        }
        globals.lv1_bucket_end = temp_end; // remember that 'end' is exclusive

        log( "Lv.1 Round %d: buckets %d to %d, num items = %lld\n", lv1_round,
             globals.lv1_bucket_start, globals.lv1_bucket_end-1, num_lv1 );
        // lv 1 bucket start & end found! fetch the offsets now

        // compute bucket positions for the offset arrays
        { // t = 0
            int64 *pos = globals.locals[0].offset_pos;
            pos[globals.lv1_bucket_start] = 0;
            for ( int i = globals.lv1_bucket_start+1; i < globals.lv1_bucket_end; ++i ) {
                pos[i] = pos[i-1] + globals.buckets[i-1];
            }
        }

        for ( int t = 1; t < globals.num_cpu_threads; ++t ) { // t > 0
            int64 *pos = globals.locals[t].offset_pos + globals.lv1_bucket_start;
            int64 *pre_pos = globals.locals[t-1].offset_pos + globals.lv1_bucket_start;
            int *pre_buckets = globals.section_buckets[t-1] + globals.lv1_bucket_start;
            for ( int i = globals.lv1_bucket_start; i < globals.lv1_bucket_end; ++i ) {
                *(pos++) = *(pre_pos++) + *(pre_buckets++);
            }
        }

        log( "  Fetching offsets... " );
        for ( int t = 1; t < globals.num_cpu_threads; ++t ) {
            pthread_create( &globals.threads[t], NULL, lv1_fetch_offsets, &globals.locals[t] );
        }
        lv1_fetch_offsets( &globals.locals[0] );
        for ( int t = 1; t < globals.num_cpu_threads; ++t ) {
            pthread_join( globals.threads[t], NULL );
        }
        log( "Done. Number of large diff: %lld\n", globals.v_large_diff.size());

        // level 2 loop
        int lv2_round = 0;
        globals.lv2_offset_start_index = 0;
        globals.lv2_bucket_start = globals.lv1_bucket_start;
        while ( globals.lv2_bucket_start < globals.lv1_bucket_end ) {
            lv2_round++;
            int64 num_lv2 = globals.buckets[ globals.lv2_bucket_start ];
            if (num_lv2 > globals.max_lv2_entries) { // TROUBLESOME CASES
                fprintf( stderr, "ERROR: Single bucket too large for lv.2: %lld (", num_lv2 );
                debug_prefix( stderr, globals.lv2_bucket_start );
                fprintf( stderr, ")\n");
                globals.lv2_bucket_start++;
                globals.lv2_offset_start_index += num_lv2;
                continue;
            }
            int temp_end2 = globals.lv2_bucket_start + 1; // remember that 'end' is exclusive
            while ( temp_end2 < globals.lv1_bucket_end && num_lv2 + globals.buckets[ temp_end2 ] <= globals.max_lv2_entries ) {
                num_lv2 += globals.buckets[ temp_end2++ ];
            }
            globals.lv2_bucket_end = temp_end2; // remember that 'end' is exclusive
            
            if (num_lv2 == 0) { // TROUBLESOME CASES
                fprintf( stderr, "ERROR: Single bucket too large for lv.2: %d (", globals.buckets[temp_end2] );
                debug_prefix( stderr, temp_end2 );
                fprintf( stderr, ")\n");
                globals.lv2_bucket_start = temp_end2 + 1;
                globals.lv2_offset_start_index += globals.buckets[temp_end2];
                continue;
            }
            
            {        
                log("  Lv.2 Round %d: buckets %d to %d, num items = %lld\n", lv2_round,
                    globals.lv2_bucket_start, globals.lv2_bucket_end-1, num_lv2 );
                // lv 2 bucket start & end found! fetch the suffixes now

                // distribute to threads
                int64 each = num_lv2 / ( globals.num_cpu_threads - 1 );
                each = !each ? 1 : each; // min is 1
                int b_start = globals.lv2_bucket_start;
                int64 cumulate_num_items = 0;
                for (int t = 0; t < (globals.num_cpu_threads-1) - 1; ++t) {
                    globals.suffix_fetches[t].bucket_start = b_start;
                    globals.suffix_fetches[t].index_start = cumulate_num_items;
                    int64 num_items = 0;
                    while ( num_items < each && b_start < globals.lv2_bucket_end ) {
                        num_items += globals.buckets[ b_start++ ];
                    }
                    globals.suffix_fetches[t].bucket_end = b_start;
                    cumulate_num_items += num_items;
                }
                globals.suffix_fetches[ (globals.num_cpu_threads-1) - 1 ].bucket_start = b_start;
                globals.suffix_fetches[ (globals.num_cpu_threads-1) - 1 ].index_start = cumulate_num_items;
                globals.suffix_fetches[ (globals.num_cpu_threads-1) -1 ].bucket_end = globals.lv2_bucket_end;

                memset( globals.lv2_suffixes, 0xFF, globals.lv2_suffixes_row_width * globals.words_per_suffix * sizeof(word) ); // important: must fill all '1' bits
                log( "   -> Fetching suffixes... " );

                for ( int t = 0; t < globals.num_cpu_threads-1; ++t ) {
                    pthread_create( &globals.threads[t], NULL, lv2_fetch_suffixes, &globals.suffix_fetches[t] );
                }
                for ( int t = 0; t < globals.num_cpu_threads-1; ++t ) {
                    pthread_join( globals.threads[t], NULL );
                }
            
                // suffixes fetched! now sort them with GPU
                /////////////////////////////////// GPU SORTING ////////////////////////////////
                log( "Sorting suffixes with CPU ... " );
#ifdef DISABLE_GPU
				omp_set_num_threads(globals.num_cpu_threads - 1);
				lv2_cpu_sort(globals.lv2_suffixes, globals.sort_indexes, globals.cpu_sort_space, globals.words_per_suffix, globals.lv2_suffixes_row_width, num_lv2);
				omp_set_num_threads(globals.num_cpu_threads);
#else
				lv2_gpu_sort(globals.lv2_suffixes, globals.sort_indexes, globals.words_per_suffix, globals.lv2_suffixes_row_width, num_lv2);
#endif
                /////// output here ///////
                if (output_started) {
                    pthread_join(output_thread, NULL);
                }
                swap(globals.sort_indexes, globals.sort_indexes_to_output);
                swap(globals.lv2_suffixes, globals.lv2_suffixes_to_output);
                swap(globals.lv2_read_id, globals.lv2_read_id_to_output);
                globals.num_lv2_to_output = num_lv2;
                pthread_create(&output_thread, NULL, lv2_output_thread, (void*)&globals);
                output_started = true;
            }

            globals.lv2_bucket_start = globals.lv2_bucket_end;
            globals.lv2_offset_start_index += num_lv2;
            log( "Round ends.\n" );
        }

        globals.lv1_bucket_start = globals.lv1_bucket_end;
        log( "Lv.1 Round %d ends\n\n", lv1_round );
    }
    // END level 1 loop

    if (output_started) {
        pthread_join(output_thread, NULL);
    }
    // print last BWT segment (those chars corresponding to '$')
    while (globals.cur_special_handle_prefix < SPECIAL_HANDLE_SUFFIX_NUM_BUCKETS) {
        if (is_valid_special_prefix(globals.cur_special_handle_prefix)) {
            dump_special_prefix(globals.cur_special_handle_prefix);
            for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
                fwrite( globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].c_str(), sizeof(char),
                        globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].size(), globals.out_file );
                globals.locals[t].bwt_buffer[globals.cur_special_handle_prefix].clear();
            }
        }
        ++globals.cur_special_handle_prefix;
    }

    fclose( globals.out_file );
    if (need_ridt)
        fclose( globals.map_out_file );
}

void read_input_file( struct global_data_t &globals, char* filename ) {
    /////////////////// READ INPUT (FASTQ) /////////////////////
    int read_len = globals.read_len;
    int words_per_read = globals.words_per_read;
    word *packed_reads = NULL;
    int64 num_reads = 0;

#define READ_BUFFER_SIZE 4096
    gzFile in_file = gzopen( filename, "r" );
    int num_bytes = 0;
    char buffer[ READ_BUFFER_SIZE + 10 ];
    packed_reads = (word*) malloc( globals.max_num_reads * words_per_read * sizeof(word) );
    // memset( packed_reads, 0xFF, globals.max_num_reads * words_per_read * sizeof(word) );
    word *packed_reads_p = packed_reads;
    int last_shift = read_len % CHARS_PER_WORD ? ( CHARS_PER_WORD - read_len % CHARS_PER_WORD ) * BITS_PER_CHAR : 0;

#define REFILL_BUFFER( _in_file, _buffer, _num_bytes )                  \
    do {                                                                \
        (_num_bytes) = gzread((_in_file),(_buffer),READ_BUFFER_SIZE);   \
        (_buffer)[(_num_bytes)] = 0;                                    \
    } while (0)

    REFILL_BUFFER( in_file, buffer, num_bytes );
    char *p = buffer;
    word *read_start_ptr;
    int is_FASTQ = ( buffer[0] == '@' );
    char read_terminator = is_FASTQ ? '+' : '>';
    while ( num_bytes ) {
        read_start_ptr = packed_reads_p;
        // consume read name @xxxxxxx
        p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
        while (!p) { 
            REFILL_BUFFER( in_file, buffer, num_bytes );
            if ( !num_bytes ) break; // end of file
            p = (char*) memchr( buffer, '\n', num_bytes );
        }
        if ( !num_bytes ) break; // end of file
        if (! *(++p) ) {
            REFILL_BUFFER( in_file, buffer, num_bytes );
            p = buffer;
        }
        // consume read sequence ACGT
        word w = 0;
        int index = 0;
        while ( *p != read_terminator ) { // FASTQ: +  FASTA:  >
            if ( *p >= 'A' ) {
                if ( index % CHARS_PER_WORD == 0 ) {
                    if ( index ) {
                        *packed_reads_p = w;
                        packed_reads_p++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_CHAR) | dna_map[(unsigned char)( *p) ]; // order: MSB to LSB
                index++;
            }
            p++;
            if (!(*p)) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break;
                p = buffer;
            }
        }
        // write last word
        *packed_reads_p = w << last_shift;
        *packed_reads_p |= ( 0xFFFFFFFF >> ( BITS_PER_WORD - last_shift ) ); // TODO hard. add padding
        packed_reads_p++;

        if ( is_FASTQ ) {
            // consume + sign
            p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
            while (!p) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                p = (char*) memchr( buffer, '\n', num_bytes );
            }
            p++;
            // consume base quality
            p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
            while (!p) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break; // end of file
                p = (char*) memchr( buffer, '\n', num_bytes );
            }
            p++; 
        } // end FASTQ only

        if (index == globals.read_len) {
            ++num_reads;
        } else { // trouble...
            log( "ERROR: read length of read #%lld is %d\n", num_reads, index);
            packed_reads_p = read_start_ptr;
        }

        if ( !num_bytes ) break; // end of file
    }

    gzclose( in_file );
        
    globals.packed_reads = packed_reads;
    globals.num_reads = num_reads;

    ////////////////////// END READING INPUT //////////////////////
}

void read_single_file(FILE *fq,  char* filename,  int64* read_num, int read_len, int threshold) 
{
    /////////////////// READ INPUT (FASTQ) /////////////////////
    int64 num_reads = 0;

#define READ_BUFFER_SIZE 4096
    gzFile in_file = gzopen( filename, "r" );
    //#define GZBUFFER_SIZE 65536        
    //    gzbuffer( infile, GZBUFFER_SIZE ); // needs zlib 1.2.4
    int num_bytes = 0;
    char buffer[ READ_BUFFER_SIZE + 10 ];
    
    //memset( packed_quals, 0, 1 * words_per_qual * sizeof(word) );
    
#define REFILL_BUFFER( _in_file, _buffer, _num_bytes )                  \
    do {                                                                \
        (_num_bytes) = gzread((_in_file),(_buffer),READ_BUFFER_SIZE);   \
        (_buffer)[(_num_bytes)] = 0;                                    \
    } while (0)

    REFILL_BUFFER( in_file, buffer, num_bytes );
    char *p = buffer;
    word *read_start_ptr;
    int is_FASTQ = ( buffer[0] == '@' );
    char read_terminator = is_FASTQ ? '+' : '>';
    while ( num_bytes ) {
        read_start_ptr = packed_reads_p;
        // consume read name @xxxxxxx
        p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
        while (!p) { 
            REFILL_BUFFER( in_file, buffer, num_bytes );
            if ( !num_bytes ) break; // end of file
            p = (char*) memchr( buffer, '\n', num_bytes );
        }
        if ( !num_bytes ) break; // end of file
        if (! *(++p) ) {
            REFILL_BUFFER( in_file, buffer, num_bytes );
            p = buffer;
        }
        // consume read sequence ACGT
        word w = 0;
        int index = 0;
        while ( *p != read_terminator ) { // FASTQ: +  FASTA:  >
            if ( *p >= 'A' ) {
                if ( index % CHARS_PER_WORD == 0 ) {
                    if ( index ) {
                        *packed_reads_p = w;
                        packed_reads_p++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_CHAR) | dna_map[(unsigned char)( *p )]; // order: MSB to LSB
                index++;
            }
            p++;
            if (!(*p)) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break;
                p = buffer;
            }
        }
        // write last word
        *packed_reads_p = w << last_shift;
        *packed_reads_p |= ( 0xFFFFFFFF >> ( BITS_PER_WORD - last_shift ) ); // TODO hard. add padding
        packed_reads_p++;

        if ( is_FASTQ ) {
            word *packed_quals_p = packed_quals;
            // consume + sign
            p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
            while (!p) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                p = (char*) memchr( buffer, '\n', num_bytes );
            }
            p++;
            if (!(*p)) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break;
                p = buffer;
            }
            // consume base quality
            w = 0;
            index = 0;
            while( *p != '@' && *p != '\n')
            {
                //fprintf(stderr, ";%c,%d", *p, w);
                if(index % BITS_PER_WORD == 0){
                    if(index){
                        *packed_quals_p = w;
                        packed_quals_p++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_QUAL) | (*p >= threshold);
                index++;
                p++;
                if (!(*p)) {
                    REFILL_BUFFER( in_file, buffer, num_bytes );
                    if ( !num_bytes ) break;
                    p = buffer;
                }
            }
            *packed_quals_p = w << last_shift_qual;
            //*packed_quals_p |= (0xFFFFFFFF >> ( BITS_PER_QUAL - last_shift_qual ));
            packed_quals_p ++;
            fwrite( packed_quals, sizeof(word), 1 * words_per_qual, fq);
            //memset( packed_quals, 0, 1 * words_per_qual * sizeof(word) );
        } // end FASTQ only

        if (index == read_len) {
            ++num_reads;
        } else { // trouble...
            log( "SINGLE READ ERROR: read length of read #%lld is %d\n", num_reads, index);
            packed_reads_p = read_start_ptr;
        }
            
        if ( !num_bytes ) break; // end of file
    }

    gzclose( in_file );
        
    (*read_num) = (*read_num) + num_reads;
    ////////////////////// END READING INPUT //////////////////////
}

void read_pair_file(FILE *fq, char* file1, char* file2,  int64* read_num, int read_len, int threshold)
{
    /////////////////// READ INPUT (FASTQ) /////////////////////
    int64 num_reads = 0;

#define READ_BUFFER_SIZE 4096
    gzFile in_file = gzopen( file1, "r" );
    gzFile in_file2 = gzopen( file2, "r" );
    //#define GZBUFFER_SIZE 65536        
    //    gzbuffer( infile, GZBUFFER_SIZE ); // needs zlib 1.2.4
    int num_bytes = 0, num_bytes2 = 0;
    char buffer[ READ_BUFFER_SIZE + 10 ];
    char buffer2[ READ_BUFFER_SIZE + 10 ];
    
   // memset( packed_quals, 0, 1 * words_per_qual * sizeof(word) );
#define REFILL_BUFFER( _in_file, _buffer, _num_bytes )                  \
    do {                                                                \
        (_num_bytes) = gzread((_in_file),(_buffer),READ_BUFFER_SIZE);   \
        (_buffer)[(_num_bytes)] = 0;                                    \
    } while (0)

    REFILL_BUFFER( in_file, buffer, num_bytes );
    REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
    char *p = buffer;
    char *p2 = buffer2;
    word *read_start_ptr;
    int is_FASTQ = ( buffer[0] == '@' );
    char read_terminator = is_FASTQ ? '+' : '>';
    while ( num_bytes && num_bytes2) {
        read_start_ptr = packed_reads_p;
        // consume read name @xxxxxxx
        p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
        p2 = (char*) memchr( p2, '\n', num_bytes2 - (p2-buffer2) );
        /*while (!p || !p2) { 
            REFILL_BUFFER( in_file, buffer, num_bytes );
            REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
            if ( !num_bytes || !num_bytes2) break; // end of file
            p = (char*) memchr( buffer, '\n', num_bytes );
            p2 = (char*) memchr( buffer2, '\n', num_bytes2 );
        }*/
        while(!p)
        {
            REFILL_BUFFER( in_file, buffer, num_bytes );
            if ( !num_bytes)
                break;
            p = (char*) memchr( buffer, '\n', num_bytes );
        }
        while(!p2)
        {
            REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
            if( !num_bytes2)
                break;
            p2 = (char*) memchr( buffer2, '\n', num_bytes2 );
        }
        if ( !num_bytes || !num_bytes2) break; // end of file
        // consume read sequence ACGT
        word w = 0;
        int index = 0;
        //fprintf(stderr, "R1:");
        //if(*p == '@' || *p == '>')
        {
            while(*p != '\n')
            {
                p++;
                if (!(*p)) {
                    REFILL_BUFFER( in_file, buffer, num_bytes );
                    if ( !num_bytes ) break;
                    p= buffer;
                }
            }
            p++;
            if (!(*p)) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break;
                p= buffer;
            }

        }
        while ( *p != read_terminator ) { // FASTQ: +  FASTA:  >
            //fprintf(stderr, "%c", *p);
            if ( *p >= 'A') {
                if ( index % CHARS_PER_WORD == 0 ) {
                    if ( index ) {
                        *packed_reads_p = w;
                        packed_reads_p++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_CHAR) | dna_map[(unsigned char)(*p) ]; // order: MSB to LSB
                index++;
            }
            p++;
            if (!(*p)) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                if ( !num_bytes ) break;
                p = buffer;
            }
        }
        
        if (index != read_len) {
            log( "ERROR: read length of read1 #%lld is %d\n", num_reads, index);
            p -= index;
            int c=0;
            while(c < index)
            {
                log( " %d,%c", c, *p);
                p++;
                c++;
            }
            log( "\n");
            exit(1);
        }
        //fprintf(stderr, "\nR2:");
        // write last word
        *packed_reads_p = w << last_shift;
        *packed_reads_p |= ( 0xFFFFFFFF >> ( BITS_PER_WORD - last_shift ) ); // TODO hard. add padding
        packed_reads_p++;

        w = 0;
        index = 0;
        //if(*(p2) == '@' || *(p2) == '>')
        {
            while(*p2 != '\n')
            {
                p2++;
                if (!(*p2)) {
                    REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                    if ( !num_bytes2 ) break;
                    p2 = buffer2;
                }
            }
            p2++;
            
            if (!(*p2)) {
                REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                if ( !num_bytes2 ) break;
                p2 = buffer2;
            }
        }
       
        while ( *p2 != read_terminator ) { // FASTQ: +  FASTA:  >
            //fprintf(stderr, "%c", *p2);
            if ( *p2 >= 'A') {
                if ( index % CHARS_PER_WORD == 0 ) {
                    if ( index ) {
                        *packed_reads_p = w;
                        packed_reads_p++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_CHAR) | dna_map[(unsigned char) (*p2) ]; // order: MSB to LSB
                index++;
            }
            p2++;
            if (!(*p2)) {
                REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                if ( !num_bytes2 ) break;
                p2 = buffer2;
            }
        }
        //fprintf(stderr, "\nq1:");
        // write last word
        *packed_reads_p = w << last_shift;
        *packed_reads_p |= ( 0xFFFFFFFF >> ( BITS_PER_WORD - last_shift ) ); // TODO hard. add padding
        packed_reads_p++;
        
        if (index != read_len) {
            log( "READ ERROR: read length of read2 #%lld is %d\n", num_reads, index);
            p2 -= index;
            int c=0;
            while(c < index)
            {
                log( " %d,%c", c, *p2);
                p2++;
                c++;
            }
            log( "\n");
            exit(1);
        }
        if ( is_FASTQ ) {
            word *packed_quals_p1 = packed_quals;
            word *packed_quals_p2 = packed_quals;
            // consume + sign
            p = (char*) memchr( p, '\n', num_bytes - (p-buffer) );
            while (!p) {
                REFILL_BUFFER( in_file, buffer, num_bytes );
                p = (char*) memchr( buffer, '\n', num_bytes );
            }
            //p++;
            while( *p != '\n')
            {
               p++;
               if(!(*p)){
                   REFILL_BUFFER( in_file, buffer, num_bytes );
                   if ( !num_bytes ) break;
                   p = buffer;
               }
            }
            p++;
            
            if(!(*p)){
               REFILL_BUFFER( in_file, buffer, num_bytes );
               if ( !num_bytes ) break;
               p = buffer;
            }
            // consume base quality
            w = 0;
            index = 0;
            while( *p != '\n')
            {
                //fprintf(stderr, "%c", *p);
                if(index % BITS_PER_WORD == 0){
                    if(index){
                        *packed_quals_p1 = w;
                        packed_quals_p1++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_QUAL) | (*p >= threshold);
                index++;
                p++;
                if (!(*p)) {
                    REFILL_BUFFER( in_file, buffer, num_bytes );
                    if ( !num_bytes ) break;
                    p = buffer;
                }
            }
            p++;
            //fprintf(stderr, "\nq2:");
            *packed_quals_p1 = w << last_shift_qual;
            //*packed_quals_p1 |= (0xFFFFFFFF >> ( BITS_PER_QUAL - last_shift_qual ));
            packed_quals_p1 ++;
            fwrite( packed_quals, sizeof(word), 1 * words_per_qual, fq);
            //memset( packed_quals, 0, 1 * words_per_qual * sizeof(word) );
            // consume + sign
            p2 = (char*) memchr( p2, '\n', num_bytes2 - (p2-buffer2) );
            while (!p2) {
                REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                p2 = (char*) memchr( buffer2, '\n', num_bytes2 );
            }
            //p2++;
            while( *p2 != '\n')
            {
                p2++;
                if (!(*p2)) {
                   REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                    if ( !num_bytes2 ) break;
                    p2 = buffer2;
                }
            }
            p2++;
             
            if (!(*p2)) {
                REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                if ( !num_bytes2 ) break;
                p2 = buffer2;
            }
            // consume base quality
            w = 0;
            index = 0;
            while( *p2 != '\n')
            {
                //fprintf(stderr, "%c", *p2);
                if(index % BITS_PER_WORD == 0){
                    if(index){
                        *packed_quals_p2 = w;
                        packed_quals_p2++;
                        w = 0;
                    }
                }
                w = (w << BITS_PER_QUAL) | (*p2 >= threshold);
                index++;
                p2++;
                if (!(*p2)) {
                    REFILL_BUFFER( in_file2, buffer2, num_bytes2 );
                    if ( !num_bytes2 ) break;
                    p2 = buffer2;
                }
            }
            p2++;
            //fprintf(stderr, "\n");
            *packed_quals_p2 = w << last_shift_qual;
            //*packed_quals_p2 |= (0xFFFFFFFF >> ( BITS_PER_QUAL - last_shift_qual ));
            packed_quals_p2 ++;
            fwrite( packed_quals, sizeof(word), 1 * words_per_qual, fq);
            //memset( packed_quals, 0, 1 * words_per_qual * sizeof(word) );
        } // end FASTQ only
        if (index == read_len) {
            num_reads += 2;
        } else { // trouble...
            log( "READ ERROR: read length of read #%lld is %d\n", num_reads, index);
            p2 -= index;
            int c=0;
            while(c < index)
            {
                log( " %d,%c", c, *p2);
                p2++;
                c++;
            }
            log( "\n");
            exit(1);
            packed_reads_p = read_start_ptr;
        }
            
        if ( !num_bytes ) break; // end of file
    }

    gzclose( in_file );

    (*read_num) = (*read_num) + num_reads;
    ////////////////////// END READING INPUT //////////////////////
}

void lib2read( struct global_data_t &globals, char* lib_file, char* prefix)
{
    int num_libs=0, i, j, read_length = globals.read_len;
    int64 start_num=0, g_num_reads=0;
    
    g_packed_reads = (word*) malloc( globals.max_num_reads * globals.words_per_read * sizeof(word) );
    //memset( g_packed_reads, 0, globals.max_num_reads * globals.words_per_read * sizeof(word) );
    last_shift = globals.read_len % CHARS_PER_WORD ? ( CHARS_PER_WORD - globals.read_len % CHARS_PER_WORD ) * BITS_PER_CHAR : 0;
    packed_reads_p = g_packed_reads;
    g_num_reads=0;
    
    words_per_qual = ( globals.read_len + ( BITS_PER_WORD - 1 ) ) / BITS_PER_WORD;
    packed_quals = (word*) malloc( 1 * words_per_qual * sizeof(word) );
    last_shift_qual = globals.read_len % BITS_PER_WORD ? ( BITS_PER_WORD - globals.read_len % BITS_PER_WORD) * BITS_PER_QUAL : 0;
    
    log("Starting load libraries...\n");
    scan_libInfo(lib_file, &num_libs, &read_length);

    if(read_length != globals.read_len)
    {
        log( "Ignore max read length in read library file\n");
        read_length = globals.read_len;
    }
    
    char name[255];
    sprintf(name, "%s.qual.bit", prefix);
    FILE *fq = fopen(name, "w"); //output base quality as 0/1.
    
    sprintf(name, "%s.index.peGrade", prefix);
    FILE *fg = fopen(name, "w"); //output insert size grade.
   
    log( "Totally there are %d libs\n", num_libs); 
    for( i=0; i< num_libs; i++ )
    {
        if(lib_array[i].asm_flag != 1 && lib_array[i].asm_flag != 3)
            continue;
        if(lib_array[i].avg_ins == 0)
            continue;
        log( "lib %d has %d paired fastqs\n", i, lib_array[i].num_q1_file);
        start_num = g_num_reads; 
        for( j = 0; j < lib_array[i].num_q1_file; j++ )
        {
            log( "Loading file: %s ", lib_array[i].q1_fname[j]);
            read_pair_file(fq, lib_array[i].q1_fname[j], lib_array[i].q2_fname[j], &g_num_reads, globals.read_len, lib_array[i].qual);
            log( ", Using %lld reads\n", g_num_reads);
        }    
        for( j = 0; j < lib_array[i].num_a1_file; j++ )
        {
            log( "Loading file: %s ", lib_array[i].a1_fname[j]);
            read_pair_file(fq, lib_array[i].a1_fname[j], lib_array[i].a2_fname[j], &g_num_reads, globals.read_len, lib_array[i].qual);
            log( ", Using %lld reads\n", g_num_reads);
        }    
        for( j = 0; j < lib_array[i].num_p_file; j++)
        {
            log( "Loading file: %s ", lib_array[i].p_fname[j]);
            read_single_file(fq, lib_array[i].p_fname[j], &g_num_reads, globals.read_len, lib_array[i].qual);
            log( ", Using %lld reads\n", g_num_reads);
        }

        if(g_num_reads > start_num)
            fprintf(fg, "%lld %lld %d %d\n", start_num, g_num_reads, lib_array[i].avg_ins, lib_array[i].sd);
    }
    
    log( "Starting Single-End libraries to reads...\n");
    for( i=0; i< num_libs; i++ )
    {
        if(lib_array[i].asm_flag != 1 && lib_array[i].asm_flag != 3)
            continue;
        log( "lib %d has %d single fastqs\n", i, lib_array[i].num_s_q_file);
        for( j = 0; j < lib_array[i].num_s_q_file; j++)
	{
            read_single_file(fq, lib_array[i].s_q_fname[j], &g_num_reads, globals.read_len, lib_array[i].qual);
        }   
        for( j = 0; j < lib_array[i].num_s_a_file; j++)
            read_single_file(fq, lib_array[i].s_a_fname[j], &g_num_reads, globals.read_len, lib_array[i].qual);
    }
 
    globals.packed_reads = g_packed_reads;
    globals.num_reads = g_num_reads;
    
    fclose(fq);
    fclose(fg);
    free_libs(lib_array, num_libs);    
    free( packed_quals );
    log("Loading library finished.\n");
}

int main( int argc, char** argv ) {
    struct sysinfo s_info;
    sysinfo(&s_info);

    if ( argc <3 ) {
#ifdef DISABLE_GPU
        log( "Usage:\n  %s read_library read_len output_prefix [thread_num] [max_host_mem] 2> bwt.log\n", argv[0] );
        log( "Input:\n 1) read_library has the same format to that of SOAPdenovo, except to set qual_cutoff to define high quality threshold\n 2) read length must be set and reads only with fixed length could be used.\n");
#else
        log( "Usage:\n  %s read_library read_len output_prefix [thread_num] [max_host_mem] [GPU_mem] 2> bwt.log\n", argv[0] );
        log( "Input:\n 1) read_library has the same format to that of SOAPdenovo, except to set qual_cutoff to define high quality threshold\n 2) read length must be set and reads only with fixed length could be used.\n");
#endif
        log( "Output:\n output_prefix.f.bwt.ascii, output_prefix.r.bwt.ascii, output_prefix.ridt\n" );
        log("\nMachine Information:\n"
               "RAM: total %lld / free %lld\n"
               "Swap: total %lld / free %lld\n"
               "Number of CPU = %d\n",
               s_info.totalram, s_info.freeram,
               s_info.totalswap, s_info.freeswap,
               sysconf(_SC_NPROCESSORS_ONLN)
              );
#ifdef DISABLE_GPU
        exit(1);
#else
	get_gpu_mem();
	exit(1);
#endif
    }

    init_dna_map();
    
    int read_len = atoi(argv[2]);
    char *output_prefix = argv[3];
    uint64_t gpu_mem = 0; // = atoll(argv[4]);
    uint64_t max_host_mem = s_info.totalram;
    int cpu_count = sysconf(_SC_NPROCESSORS_ONLN);

    if(argc > 4)
    {
        number_cpu_threads = atoi(argv[4]);
        if(number_cpu_threads > cpu_count || number_cpu_threads == 0 || number_cpu_threads > NUM_CPU_THREADS)
            number_cpu_threads = cpu_count < NUM_CPU_THREADS ? cpu_count : NUM_CPU_THREADS;
        log("Using thread number: %d\n", number_cpu_threads);
    }else{
        number_cpu_threads = cpu_count < NUM_CPU_THREADS ? cpu_count : NUM_CPU_THREADS;
        log("Using thread number: %d\n", number_cpu_threads);
    }

    if(argc > 5)
    {
        uint64_t set_host_mem = atoll(argv[5]);
        if(set_host_mem > max_host_mem || set_host_mem == 0)
        {
            log("Use maximum host memory %lld\n", max_host_mem);
        }else{
            max_host_mem = set_host_mem;
        }
    }else{
        log("Using total RAM %lld\n", max_host_mem);
    }

#ifdef DISABLE_GPU
    gpu_mem = 0;
#else
    size_t free_gpu_mem = get_gpu_mem();
    if(gpu_mem < 100 || gpu_mem > free_gpu_mem)
    {
        gpu_mem = free_gpu_mem;
    }
#endif

    if(argc > 6)
    {
#ifdef DISABLE_GPU
        gpu_mem = atoll(argv[6]);
#else
        uint64_t set_gpu_mem = atoll(argv[6]);
        if(set_gpu_mem < 100 || set_gpu_mem > gpu_mem)
        {
            log("Use Maximum GPU memory %lld.\n", gpu_mem);
        }else{
            gpu_mem = set_gpu_mem;
            log("Use Set GPU memory %lld.\n", gpu_mem);
        }
#endif
    }else{
#ifdef DISABLE_GPU
        log("Use maximum bucket size.\n");
#else
        log("Use Maximum GPU memory %lld.\n", gpu_mem);
#endif
    }

    struct global_data_t globals;
    
    init_global_data(globals, read_len, max_host_mem, number_cpu_threads);

    // read input
    lib2read( globals, argv[1], output_prefix);
    log("Number of reads: %lld\n", globals.num_reads);
    
    //building bwt.
    //set_global_after_reading(globals, max_host_mem);
    set_global_bucket(globals);
    //level 1 fill buckets.
    uint64_t max_bucket_size = lv0_fill_buckets(globals);
    
    //set level2 memeory.
    set_global_after_lv0(globals, max_host_mem, gpu_mem, max_bucket_size);

    build_bwt(globals, output_prefix, true);

    log("Reverse all reads..");
    reverse_all_reads(globals);
    log("Done\n");

    lv0_fill_buckets(globals);
    build_bwt(globals, output_prefix, false);

    log("All finished.\n");   
    // free memory

    free( globals.packed_reads );
    free( globals.buckets );
    free( globals.lv1_offsets );    
    free( globals.lv2_suffixes );
    free( globals.lv2_read_id );
    free( globals.sort_indexes );

#ifdef DISABLE_GPU
    free(globals.cpu_sort_space);
#endif

    pthread_spin_destroy(&globals.large_diff_lock);

    for ( int t = 0; t < globals.num_cpu_threads; ++t ) {
        free( globals.section_buckets[t] );
        free( globals.locals[t].offset_pos );
    }

    return 0;

}

/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#define _GNU_SOURCE
#include "thread.h"
#include "list.h"
#include "string.h"
#include <regex.h>
#include <unistd.h>

typedef struct {
	char *name;
	char *prefix;
	int ins_size;
	int pair_idx;
	int is_fq;
} SeqFile;

define_list(sfv, SeqFile);

static inline void parse_seqfile(SeqFile *sf, char *name){
	regex_t pat;
	regmatch_t mats[5];
	int ret;
	regcomp(&pat, "(.*?)\\.ins([0-9]+)\\.s([12])\\.f([aq])", REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	ret = regexec(&pat, name, 5, mats, 0);
	if(ret == REG_NOMATCH){
		fprintf(stdout, " -- Wrong file name \"%s\" in %s -- %s:%d --\n", name, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
	}
	sf->name = catstr(1, name);
	strcpy(sf->name, name);
	sf->prefix = malloc(mats[1].rm_eo - mats[1].rm_so + 1);
	memcpy(sf->prefix, sf->name + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so);
	sf->prefix[mats[1].rm_eo - mats[1].rm_so] = 0;
	sf->ins_size = strtol(sf->name + mats[2].rm_so, NULL, 10);
	sf->pair_idx = strtol(sf->name + mats[3].rm_so, NULL, 10);
	sf->is_fq    = (sf->name[mats[4].rm_so] == 'q');
	regfree(&pat);
}

static inline int cmp_seqfile(const void *e1, const void *e2){
	SeqFile *s1, *s2;
	int cmp;
	s1 = (SeqFile*)e1;
	s2 = (SeqFile*)e2;
	if(s1->ins_size == s2->ins_size){
		cmp = strcmp(s1->prefix, s2->prefix);
		if(cmp == 0){
			if(s1->pair_idx == s2->pair_idx) return 0;
			else if(s1->pair_idx < s2->pair_idx) return -1;
			else return 1;
		} else return cmp;
	} else if(s1->ins_size < s2->ins_size) return 1;
	else return -1;
}

typedef struct {
	int      n_cpu;
	uint32_t limit;
	uint32_t min_lib, max_lib;
	uint32_t min_ins, max_ins;
	uint32_t low_cpx;
	uint32_t rd_len;
	uint32_t mode;
	uint32_t gap_cutoff;
	uint32_t n_seed[2];
	uint32_t kmer_size[2];
	uint32_t min_ol[2];
	float    min_sm[2];
	uint32_t max_mm[2];
	char     *outf;
} ATOptions;

static inline int parse_options(ATOptions *opt, int mode, int argc, char **argv){
	int c;
	opt->n_cpu   = 1;
	opt->limit   = 400;
	opt->min_lib = 0;
	opt->min_ins = 500;
	opt->max_ins = 700;
	opt->low_cpx = 12;
	opt->rd_len  = 0;
	opt->mode    = mode;
	opt->gap_cutoff   = 100;
	opt->n_seed[0]    = 4;
	opt->n_seed[1]    = 4;
	opt->kmer_size[0] = 15;
	opt->kmer_size[1] = 4;
	opt->min_ol[0]    = 30;
	opt->min_ol[1]    = 8;
	opt->min_sm[0]    = 0.97;
	opt->min_sm[1]    = 0.90;
	opt->max_mm[0]    = 4;
	opt->max_mm[1]    = 6;
	opt->outf = NULL;
	while((c = getopt(argc, argv, "ht:c:X:Y:x:y:f:r:12n:k:l:s:m:g:o:")) != -1){
		switch(c){
			case 'h': return -1;
			case 't': opt->n_cpu = atoi(optarg); break;
			case 'c': opt->limit = atoi(optarg); break;
			case 'X': opt->min_lib = atoi(optarg); break;
			case 'x': opt->min_ins = atoi(optarg); break;
			case 'y': opt->max_ins = atoi(optarg); break;
			case 'f': opt->low_cpx = atoi(optarg); break;
			case 'r': opt->rd_len  = atoi(optarg); break;
			case '1': opt->mode = 0; break;
			case '2': opt->mode = 1; break;
			case 'n': opt->n_seed[opt->mode] = atoi(optarg); break;
			case 'k': opt->kmer_size[opt->mode] = atoi(optarg); break;
			case 'l': opt->min_ol[opt->mode] = atoi(optarg); break;
			case 's': opt->min_sm[opt->mode] = atof(optarg); break;
			case 'm': opt->max_mm[opt->mode] = atoi(optarg); break;
			case 'g': opt->gap_cutoff = atoi(optarg); break;
			case 'o': opt->outf = optarg; break;
			default: printf("Unknown option '-%c'\n", c); return -1;
		}
	}
	return optind;
}


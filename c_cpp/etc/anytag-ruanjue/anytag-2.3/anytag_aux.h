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
 
#ifndef __ANYTAG_AUX_RJ_H
#define __ANYTAG_AUX_RJ_H

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
	int var_size;
	int pair_idx;
	int is_fq;
} SeqFile;

define_list(sfv, SeqFile);

static inline void parse_seqfile(SeqFile *sf, char *name, int var_size){
	regex_t pat;
	regmatch_t mats[7];
	int ret;
	regcomp(&pat, "(.*?)\\.ins([0-9]+)(\\.var([0-9]+))?\\.s([12])\\.f([aq])", REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	ret = regexec(&pat, name, 7, mats, 0);
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
	if(mats[4].rm_so < mats[4].rm_eo) sf->var_size = strtol(sf->name + mats[4].rm_so, NULL, 10);
	else sf->var_size = var_size;
	sf->pair_idx = strtol(sf->name + mats[5].rm_so, NULL, 10);
	sf->is_fq    = (sf->name[mats[6].rm_so] == 'q');
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
	uint32_t var_size;
	uint32_t low_cpx;
	uint32_t rd_len;
	uint32_t rd_clip;
	uint32_t n_part;
	uint32_t inc_ar;
	uint32_t n_full;
	int      mode;
	uint32_t n_seed[2];
	uint32_t kmer_size[2];
	uint32_t min_ol[3];
	float    min_sm[3];
	uint32_t max_mm[3];
	uint32_t max_hits[3];
	uint32_t allow_gap[3];
	uint32_t gap_cutoff[2];
	uint32_t flags[3];
	char     *load_aln;
} ATOptions;

static inline void print_options(ATOptions *opt, FILE *out){
	fprintf(out, "Options:\n");
	fprintf(out, "-t %d\n", opt->n_cpu);
	fprintf(out, "-c %d\n", opt->limit);
	if(opt->min_lib) fprintf(out, "-X %d\n", opt->min_lib);
	fprintf(out, "-x %d\n", opt->min_ins);
	fprintf(out, "-y %d\n", opt->max_ins);
	fprintf(out, "-d %d\n", opt->var_size);
	fprintf(out, "-r %d\n", opt->rd_len);
	fprintf(out, "-R %d\n", opt->rd_clip);
	fprintf(out, "-p %d\n", opt->n_part);
	fprintf(out, "-f %d\n", opt->low_cpx);
	fprintf(out, "-S %d\n", opt->n_full);
	fprintf(out, "-1\n");
	fprintf(out, "-C %d\n", opt->max_hits[0]);
	fprintf(out, "-n %d\n", opt->n_seed[0]);
	if(opt->inc_ar) fprintf(out, "-a\n");
	fprintf(out, "-k %d\n", opt->kmer_size[0]);
	fprintf(out, "-l %d\n", opt->min_ol[0]);
	fprintf(out, "-s %f\n", opt->min_sm[0]);
	fprintf(out, "-m %d\n", opt->max_mm[0]);
	fprintf(out, "-2\n");
	fprintf(out, "-C %d\n", opt->max_hits[1]);
	fprintf(out, "-n %d\n", opt->n_seed[1]);
	fprintf(out, "-k %d\n", opt->kmer_size[1]);
	fprintf(out, "-l %d\n", opt->min_ol[1]);
	fprintf(out, "-s %f\n", opt->min_sm[1]);
	fprintf(out, "-m %d\n", opt->max_mm[1]);
	fprintf(out, "-3\n");
	fprintf(out, "-C %d\n", opt->max_hits[2]);
	fprintf(out, "-l %d\n", opt->min_ol[2]);
	fprintf(out, "-s %f\n", opt->min_sm[2]);
	fprintf(out, "-m %d\n", opt->max_mm[2]);
	if(opt->allow_gap[2] == 2) fprintf(out, "-G\n");
	if(opt->load_aln) fprintf(out, "-i %s\n", opt->load_aln);
	fflush(out);
}

static inline int parse_options(ATOptions *opt, int mode, int argc, char **argv){
	int c;
	opt->n_cpu   = 8;
	opt->limit   = 400;
	opt->min_lib = 0;
	opt->min_ins = 400;
	opt->max_ins = 800;
	opt->var_size     = 100;
	opt->low_cpx = 12;
	opt->rd_len  = 0;
	opt->rd_clip = 0;
	opt->n_part  = 1;
	opt->inc_ar  = 0;
	opt->n_full  = 10000;
	opt->allow_gap[0] = 0;
	opt->allow_gap[1] = 0;
	opt->allow_gap[2] = 0;
	opt->gap_cutoff[0]= 0;
	opt->gap_cutoff[1]= 0;
	opt->n_seed[0]    = 4;
	opt->n_seed[1]    = 4;
	opt->kmer_size[0] = 15;
	opt->kmer_size[1] = 5;
	opt->min_ol[0]    = 30;
	opt->min_ol[1]    = 10;
	opt->min_ol[2]    = 50;
	opt->min_sm[0]    = 0.97;
	opt->min_sm[1]    = 0.90;
	opt->min_sm[2]    = 0.97;
	opt->max_mm[0]    = 4;
	opt->max_mm[1]    = 6;
	opt->max_mm[2]    = 6;
	opt->max_hits[0]  = 500;
	opt->max_hits[1]  = 5000;
	opt->max_hits[2]  = 1000;
	opt->flags[0]     = 0;
	opt->flags[1]     = 0;
	opt->flags[2]     = 0;
	opt->load_aln     = NULL;
	while((c = getopt(argc, argv, "ht:c:X:Y:x:y:d:f:r:R:p:S:123i:ab:n:k:l:s:m:C:G")) != -1){
		switch(c){
			case 'h': return -1;
			case 'i': opt->load_aln = optarg; break;
			case 't': opt->n_cpu = atoi(optarg); break;
			case 'c': opt->limit = atoi(optarg); break;
			case 'X': opt->min_lib = atoi(optarg); break;
			case 'x': opt->min_ins = atoi(optarg); break;
			case 'y': opt->max_ins = atoi(optarg); break;
			case 'd': opt->var_size = atoi(optarg); break;
			case 'f': opt->low_cpx = atoi(optarg); break;
			case 'r': opt->rd_len  = atoi(optarg); break;
			case 'R': opt->rd_clip = atoi(optarg); break;
			case 'p': opt->n_part  = atoi(optarg); break;
			case 'S': opt->n_full = atoll(optarg); break;
			case 'a': opt->inc_ar = 1; break;
			case 'G': opt->allow_gap[2] = 2; break;
			case '1': mode = 0; break;
			case '2': mode = 1; break;
			case '3': mode = 2; break;
			case 'C': opt->max_hits[mode]= atoi(optarg); break;
			case 'b': opt->flags[mode] = atoi(optarg); break;
			case 'n': opt->n_seed[mode] = atoi(optarg); break;
			case 'k': opt->kmer_size[mode] = atoi(optarg); break;
			case 'l': opt->min_ol[mode] = atoi(optarg); break;
			case 's': opt->min_sm[mode] = atof(optarg); break;
			case 'm': opt->max_mm[mode] = atoi(optarg); break;
			default: printf("Unknown option '-%c'\n", c); return -1;
		}
	}
	return optind;
}

typedef struct {
	uint32_t n_lib;
	u32list  *lib_offs;
	u32list  *min_ins, *max_ins;
} LibInserts;

static inline LibInserts* guess_lib_inserts(sfv *sfs){
	LibInserts *libs;
	SeqFile *sf;
	uint32_t i;
	libs = malloc(sizeof(LibInserts));
	libs->n_lib = count_sfv(sfs) / 2;
	libs->lib_offs = init_u32list(libs->n_lib);
	libs->min_ins  = init_u32list(libs->n_lib);
	libs->max_ins  = init_u32list(libs->n_lib);
	for(i=0;i<libs->n_lib;i++){
		sf = ref_sfv(sfs, i * 2 + 0);
		if(sf->ins_size > sf->var_size) push_u32list(libs->min_ins, sf->ins_size - sf->var_size);
		else push_u32list(libs->min_ins, 0);
		push_u32list(libs->max_ins, sf->ins_size + sf->var_size);
	}
	return libs;
}

static inline uint32_t get_read_lib(LibInserts *libs, uint32_t rid){
	uint32_t i;
	for(i=1;i<libs->n_lib;i++){ if(rid < get_u32list(libs->lib_offs, i)) break; }
	return i - 1;
}

static inline void free_lib_inserts(LibInserts *libs){
	free_u32list(libs->lib_offs);
	free_u32list(libs->min_ins);
	free_u32list(libs->max_ins);
	free(libs);
}

#endif

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

static inline void parse_seqfile(SeqFile *sf, char *name, float var){
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
	else sf->var_size = sf->ins_size * var;
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
	uint32_t genome_size;
	uint32_t limit;
	uint32_t min_ar;
	uint32_t inc_ar;
	uint32_t low_cpx;
	uint32_t qc_ksize;
	uint32_t rd_trim;
	uint32_t rd_clip;
	uint32_t n_full;
	uint32_t asm_flags;
	uint8_t  _asm_flags[5];
	int      mode;
	uint32_t n_seed;
	uint32_t seed_skips;
	uint32_t kmer_size[2];
	uint32_t min_ol[2];
	float    min_sm[2];
	uint32_t flags[2];
	int      re_cns;
	int      debug;
	uint32_t debug_beg;
	char     *scaff;
} ATOptions;

static inline void print_options(ATOptions *opt, FILE *out){
	uint32_t i;
	fprintf(out, "Options:\n");
	fprintf(out, "-g %u\n", opt->genome_size);
	fprintf(out, "-t %d\n", opt->n_cpu);
	if(opt->min_ar) fprintf(out, "-X %d\n", opt->min_ar);
	if(opt->inc_ar) fprintf(out, "-A\n");
	fprintf(out, "-m %u\n", opt->asm_flags);
	fprintf(out, "-r %d\n", opt->rd_trim);
	//fprintf(out, "-R %d\n", opt->rd_clip);
	fprintf(out, "-f %d\n", opt->low_cpx);
	fprintf(out, "-q %d\n", opt->qc_ksize);
	fprintf(out, "-S %d\n", opt->n_full);
	fprintf(out, "-1\n");
	fprintf(out, "-n %d\n", opt->n_seed);
	for(i=0;i<opt->n_seed;i++){
		if((opt->seed_skips >> i) & 0x01){
			fprintf(out, "-N %d\n", i);
		}
	}
	fprintf(out, "-k %d\n", opt->kmer_size[0]);
	fprintf(out, "-l %d\n", opt->min_ol[0]);
	fprintf(out, "-s %f\n", opt->min_sm[0]);
	fprintf(out, "-2\n");
	fprintf(out, "-k %d\n", opt->kmer_size[1]);
	fprintf(out, "-l %d\n", opt->min_ol[1]);
	fprintf(out, "-s %f\n", opt->min_sm[1]);
	if(opt->re_cns) fprintf(out, "-c\n");
	if(opt->debug == 1){
	} else if(opt->debug == 2){
		fprintf(out, "-D %u\n", opt->debug_beg);
	} else if(opt->debug == 3){
		fprintf(out, "-d\n");
	}
	if(opt->scaff){
		fprintf(out, "-G %s\n", opt->scaff);
	}
	fflush(out);
}

static inline int parse_options(ATOptions *opt, int mode, int argc, char **argv){
	int i, c;
	opt->n_cpu   = 8;
	opt->genome_size = 0;
	opt->limit   = 0;
	opt->min_ar  = 0;
	opt->inc_ar  = 0;
	opt->low_cpx = 7;
	opt->rd_trim = 0;
	opt->rd_clip = 0;
	opt->qc_ksize = 0;
	opt->n_full  = 10000;
	opt->n_seed  = 4;
	opt->seed_skips   = 0;
	opt->kmer_size[0] = 15;
	opt->kmer_size[1] = 9;
	opt->min_ol[0]    = 40;
	opt->min_ol[1]    = 20;
	opt->min_sm[0]    = 0.95;
	opt->min_sm[1]    = 0.95;
	opt->flags[0]     = 0;
	opt->flags[1]     = 0;
	opt->re_cns = 0;
	opt->debug = 0;
	opt->debug_beg = 0;
	opt->scaff = NULL;
	while((c = getopt(argc, argv, "hg:t:X:m:f:q:r:R:AS:12b:n:N:k:l:s:cdD:G:")) != -1){
		switch(c){
			case 'h': return -1;
			case 'g': opt->genome_size = atol(optarg); break;
			case 't': opt->n_cpu = atoi(optarg); break;
			case 'X': opt->min_ar = atoi(optarg); break;
			case 'f': opt->low_cpx = atoi(optarg); break;
			case 'q': opt->qc_ksize = atoi(optarg); break;
			case 'r': opt->rd_trim = atoi(optarg); break;
			case 'R': opt->rd_clip = atoi(optarg); break;
			case 'm': opt->asm_flags = atoi(optarg); break;
			case 'S': opt->n_full = atoll(optarg); break;
			case 'A': opt->inc_ar = 1; break;
			case 'n': opt->n_seed = atoi(optarg); break;
			case 'N': opt->seed_skips |= 1U << (atoi(optarg) - 1); break;
			case '1': mode = 0; break;
			case '2': mode = 1; break;
			case 'b': opt->flags[mode] = atoi(optarg); break;
			case 'k': opt->kmer_size[mode] = atoi(optarg); break;
			case 'l': opt->min_ol[mode] = atoi(optarg); break;
			case 's': opt->min_sm[mode] = atof(optarg); break;
			case 'c': opt->re_cns = 1; break;
			case 'D': if(optarg[0] == '-') opt->debug = 3; else { opt->debug = 2; opt->debug_beg = atoll(optarg); } break;
			case 'd': opt->debug = 1; break;
			case 'G': opt->scaff = optarg; break;
			default: printf("Unknown option '-%c'\n", c); return -1;
		}
	}
	if(((opt->asm_flags / 1000) % 10) == 1){
		opt->asm_flags = (opt->asm_flags / 100) * 100 + 211;
		opt->n_seed = 2;
		opt->seed_skips = 0;
		opt->min_ol[0] = 60;
		opt->min_ol[0] = 40;
		opt->min_sm[0] = 1.0;
	}
	c = 1;
	for(i=0;i<5;i++){
		opt->_asm_flags[i] = (opt->asm_flags / c) % 10;
		c *= 10;
	}
	return optind;
}

typedef struct {
	uint32_t n_lib;
	u32list  *lib_offs;
	u32list  *lib_cnts;
	u32list  *min_ins, *max_ins;
	u32list  *lib_ins, *lib_vars;
} LibInserts;

static inline LibInserts* guess_lib_inserts(sfv *sfs){
	LibInserts *libs;
	SeqFile *sf;
	uint32_t i;
	libs = malloc(sizeof(LibInserts));
	libs->n_lib = count_sfv(sfs) / 2;
	libs->lib_offs = init_u32list(libs->n_lib);
	libs->lib_cnts = init_u32list(libs->n_lib);
	libs->min_ins  = init_u32list(libs->n_lib);
	libs->max_ins  = init_u32list(libs->n_lib);
	libs->lib_ins  = init_u32list(libs->n_lib);
	libs->lib_vars = init_u32list(libs->n_lib);
	for(i=0;i<libs->n_lib;i++){
		sf = ref_sfv(sfs, i * 2 + 0);
		if(sf->ins_size > sf->var_size) push_u32list(libs->min_ins, sf->ins_size - sf->var_size);
		else push_u32list(libs->min_ins, 0);
		push_u32list(libs->max_ins, sf->ins_size + sf->var_size);
		push_u32list(libs->lib_ins, sf->ins_size);
		push_u32list(libs->lib_vars, sf->var_size);
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
	free_u32list(libs->lib_cnts);
	free_u32list(libs->min_ins);
	free_u32list(libs->max_ins);
	free_u32list(libs->lib_ins);
	free_u32list(libs->lib_vars);
	free(libs);
}

#endif

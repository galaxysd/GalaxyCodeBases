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
 
#include "sr_aln.h"
#include "local_assembly.h"
#include "file_reader.h"
#include "thread.h"
#include "anytag_aux.h"

int usage_asm();

thread_begin_def(masm);
LGraph *g;
FILE *cnsf, *samf, *msaf;
int asm_ret;
thread_end_def(masm);

thread_begin_func(masm);
thread_begin_loop(masm);
align_lgraph(masm->g);
masm->asm_ret = layout_lgraph(masm->g);
if(masm->asm_ret){
	consensus_lgraph(masm->g);
	realign_lgraph(masm->g);
	//thread_begin_syn(masm);
	//output_lgraph(masm->g, masm->cnsf, masm->samf, masm->msaf);
	//thread_end_syn(masm);
}
thread_end_loop(masm);
thread_end_func(masm);

int main_asm(int argc, char **argv){
	ATOptions opt;
	FileReader *fr;
	FILE *cnsf, *samf, *msaf;
	uint32_t n, n_seq, dir;
	int c;
	char flag, *seq1, *seq2, *prefix, cmd[256];
	uint32_t pid, len1, len2;
	thread_preprocess(masm);
	if((c = parse_options(&opt, 1, argc, argv)) == -1) return usage_asm();
	if(c + 1 >= argc) return usage_asm();
	if((fr = fopen_filereader(argv[c])) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", argv[c], __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	prefix = argv[c+1];
	sprintf(cmd, "%s.fasta", prefix);
	if((cnsf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.sam.gz", prefix);
	if((samf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.msa.gz", prefix);
	if((msaf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	thread_begin_init(masm, (int)opt.n_cpu);
	masm->g = init_lgraph(opt.kmer_size[1], opt.rd_len, opt.min_ol[1], opt.min_sm[1], opt.max_mm[1], opt.gap_cutoff, opt.min_ins, opt.max_ins);
	masm->cnsf = cnsf;
	masm->samf = samf;
	masm->msaf = msaf;
	masm->asm_ret = 0;
	thread_end_init(masm);
	n = 0;
	n_seq = 0;
	thread_begin_operate(masm, 0);
	while((c = fread_table(fr)) != -1){
		if(c < 5) continue;
		flag = get_col_str(fr, 0)[0];
		pid  = atol(get_col_str(fr, 1));
		dir  = (get_col_str(fr, 2)[0] == '-');
		seq1 = get_col_str(fr, 3);
		len1 = get_col_len(fr, 3);
		seq2 = get_col_str(fr, 4);
		len2 = get_col_len(fr, 4);
		if(flag == 'T'){
			thread_wake(masm);
			thread_waitfor_one_idle(masm);
			n = thread_index(masm);
			if(masm->asm_ret) output_lgraph(masm->g, masm->cnsf, masm->samf, masm->msaf);
			masm->asm_ret = 0;
			reset_lgraph(masm->g);
			if(dir){
				push_lgraph(masm->g, 2 * pid + 1, seq1, len1, 0);
				push_lgraph(masm->g, 2 * pid + 0, seq2,  len2, 1);
			} else {
				push_lgraph(masm->g, 2 * pid + 0, seq1, len1, 0);
				push_lgraph(masm->g, 2 * pid + 1, seq2, len2, 1);
			}
			n_seq = 2;
		} else if(flag == 'S'){
			if(n_seq > opt.limit) continue;
			n_seq += 2;
			if(dir){
				push_lgraph(masm->g, 2 * pid + 1, seq1, len1, 0);
				push_lgraph(masm->g, 2 * pid + 0, seq2, len2, 1);
			} else {
				push_lgraph(masm->g, 2 * pid + 0, seq1, len1, 0);
				push_lgraph(masm->g, 2 * pid + 1, seq2, len2, 1);
			}
		}
	}
	fclose_filereader(fr);
	thread_wake(masm);
	thread_waitfor_all_idle(masm);
	thread_begin_close(masm);
	if(masm->asm_ret) output_lgraph(masm->g, masm->cnsf, masm->samf, masm->msaf);
	masm->asm_ret = 0;
	free_lgraph(masm->g);
	thread_end_close(masm);
	fclose(cnsf);
	pclose(samf);
	pclose(msaf);
	return 0;
}


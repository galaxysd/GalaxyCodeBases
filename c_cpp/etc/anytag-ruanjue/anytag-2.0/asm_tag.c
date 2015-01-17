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
 
#include "simp_asm.h"
#include "file_reader.h"
#include <unistd.h>

void asm_tag(SimpAssembler *sa, char *tag, FILE *out){
	SimpContigInfo *ctg;
	SimpSeqInfo *rd;
	uint32_t i, ctg_id;
	simple_assemble(sa);
	begin_iter_simpasm(sa);
	ctg_id = 1;
	while((ctg = iter_simpasm(sa))){
		fprintf(out, ">%s_%03u len=%u n_rd=%u rds=", tag, ctg_id++, ctg->len, (unsigned)count_u32list(ctg->sids));
		for(i=0;i<count_u32list(ctg->sids);i++){
			rd = ref_seqv(sa->rds, get_u32list(ctg->sids, i));
			if(i) fprintf(out, ",");
			fprintf(out, "%u:%u:%c:%u", rd->seqid, rd->len, "+-"[rd->ctg_dir], rd->ctg_off);
		}
		fprintf(out, "\n");
		for(i=0;i<ctg->len;i++){
			if(i && i % 100 == 0) fprintf(out, "\n");
			fprintf(out, "%c", ctg->seq->string[i]);
		}
		fprintf(out, "\n");
	}
}

int usage(){
	printf(
"Usage: asm_tag [option]\n"
"Options:\n"
" -k <int>     Kmer size 3 - 16, [4]\n"
" -l <int>     Read length, 0: auto [0]\n"
" -m <int>     Min overlap [9]\n"
" -s <int>     Similarity [0.9]\n"
" -x <int>     Max mismatches when aligning reads [4]\n"
" -g           Allow gaps when aligning reads [no]\n"
" -t <int>     Num of thread [1]\n"
		  );
	return 1;
}

int main(int argc, char **argv){
	SimpAssembler *sa;
	FileReader *fr;
	FILE *out;
	char *seq1, *seq2, *dir, tag[64];
	uint32_t seqid, ksize, rd_len, min_ol, max_mm, rank1, rank2, gap, n_cpu;
	float min_sm;
	int n_col, c;
	ksize  = 4;
	rd_len = 0;
	min_ol = 9;
	min_sm = 0.9;
	max_mm = 4;
	gap = 0;
	n_cpu = 1;
	while((c = getopt(argc, argv, "hgk:l:m:s:x:t:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'k': ksize = atoi(optarg); break;
			case 'l': rd_len = atoi(optarg); break;
			case 'm': min_ol = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'x': max_mm = atoi(optarg); break;
			case 'g': gap = 1; break;
			case 't': n_cpu = atoi(optarg); break;
			default: return usage();
		}
	}
	sa = init_simpasm(n_cpu, ksize, rd_len, 1, min_ol, min_sm, max_mm, gap);
	fr = stdin_filereader();
	out = stdout;
	tag[0] = 0;
	while((n_col = fread_table(fr)) != -1){
		if(fr->line->string[0] == '#') continue;
		if(n_col < 6) continue;
		if(strcmp(tag, get_col_str(fr, 0))){
			if(tag[0]) asm_tag(sa, tag, out);
			strcpy(tag, get_col_str(fr, 0));
			reset_simpasm(sa);
		}
		seqid = atoll(get_col_str(fr, 1));
		dir  = get_col_str(fr, 2);
		seq1 = get_col_str(fr, 4);
		seq2 = get_col_str(fr, 5);
		if(dir[0] == 'F'){
			rank1 = 2;
			rank2 = 3;
			if(dir[1] == 'C') reverse_dna(seq1, get_col_len(fr, 4));
			else reverse_dna(seq2, get_col_len(fr, 5));
		} else {
			rank1 = 1;
			rank2 = 2;
			if(dir[1] == 'C') reverse_dna(seq2, get_col_len(fr, 5));
			else reverse_dna(seq1, get_col_len(fr, 4));
		}
		push_simpasm(sa, seqid * 2 + 0, seq1, get_col_len(fr, 4), rank1);
		push_simpasm(sa, seqid * 2 + 1, seq2, get_col_len(fr, 5), rank2);
		//printf(">%u_%u\n%s\n", seqid * 2 + 0, rank1, seq1);
		//printf(">%u_%u\n%s\n", seqid * 2 + 1, rank2, seq2);
		//if(rank1 == 2) push_simpasm(sa, seqid * 2 + 0, seq1, get_col_len(fr, 4), rank1);
		//if(rank2 == 2) push_simpasm(sa, seqid * 2 + 1, seq2, get_col_len(fr, 5), rank2);
	}
	fclose_filereader(fr);
	if(tag[0]) asm_tag(sa, tag, out);
	free_simpasm(sa);
	return 0;
}

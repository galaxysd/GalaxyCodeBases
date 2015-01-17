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
#include "file_reader.h"
#include <unistd.h>

int usage(){
	printf(
			"Usage: sra [Options]\n"
			"Options:\n"
			" -k <int>    Seed size, 3 - 16 bp, [8]\n"
			" -l <int>    Min overlap [21]\n"
			" -s <float>  Min similiarity [0.9]\n"
			" -m <int>    Max mismatch [5]\n"
			" -f <int>    Min number of tri-kmers in one sequences, filtering low complexity [10]\n"
			" -g          Gap alignment (Smith-Waterman), [no]\n"
			" -S          Reads must be aligned on the same strand [no]\n"
			"\n"
			"Example: sra < rds.fa >rds.ols\n"
		  );
	return 1;
}


int main(int argc, char **argv){
	SR_SeqDB *sdb;
	SR_AlnAux *aux;
	FileReader *fr;
	Sequence *seq;
	uint32_t i, j;
	int c, kmer_size, rd_len, min_ol, max_mm, s, gap, n_cpu, min_cpx, format;
	float min_sm;
	kmer_size = 8;
	rd_len = 0;
	min_ol = 21;
	min_sm = 0.9;
	max_mm = 5;
	s = 3;
	gap = 0;
	n_cpu = 1;
	min_cpx = 10;
	format = 0;
	while((c = getopt(argc, argv, "hg01Sk:r:l:s:f:m:t:")) != -1){
		switch(c){
			case 'k': kmer_size = atoi(optarg); break;
			case 'r': rd_len = atoi(optarg); break;
			case 'l': min_ol = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'S': s = 1; break;
			case 'm': max_mm = atoi(optarg); break;
			case 'g': gap = 1; break;
			case '0': format = 0; break;
			case '1': format = 1; break;
			case 't': n_cpu = atoi(optarg); break;
			case 'f': min_cpx = atoi(optarg); break;
			default : return usage();
		}
	}
	sdb = sr_init_sdb(kmer_size, 4, rd_len, min_cpx);
	fr = stdin_filereader();
	seq = NULL;
	while(fread_fasta(&seq, fr)){
		if(seq->seq.size > MAX_RD_LEN) continue;
		sr_push_sdb(sdb, seq->seq.string, seq->seq.size);
	}
	fclose_filereader(fr);
	sr_ready_sdb(sdb);
	for(i=0;i<sdb->n_idx;i++) sr_index_sdb(sdb, i);
	aux = sr_init_aux(s, min_ol, max_mm, min_sm, gap, 1, 1024);
	sr_fit_aux2sdb(aux, sdb);
	for(i=0;i<sdb->n_rd;i++){
		sr_clear_aux(aux);
		sr_align_sdb(sdb, i, aux);
		for(j=0;j<count_sr_hitv(aux->hits);j++){
			sr_output_hit(sdb, ref_sr_hitv(aux->hits, j), stdout);
		}
	}
	sr_free_sdb(sdb);
	sr_free_aux(aux);
	return 0;
}

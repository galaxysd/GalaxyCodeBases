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
			" -0          Output as 'rid1 strand1 rid2 strand2 cigar seq1 seq2' [yes]\n"
			" -1          Output AMOS overlap_t format [no]\n"
			" -t <int>    Number of thread [1]\n"
			"\n"
			"Example: sra < rds.fa >rds.ols\n"
		  );
	return 1;
}

static const char __amos_adj_types[] = "NIOA";

void output_amos_overlap_message(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	char adj;
	int ahg, bhg;
	adj = __amos_adj_types[(((int)hit->dir1) << 1) | ((int)hit->dir2)];
	ahg = hit->off;
	bhg = ((int)(sr_rdseq_length(sdb, hit->rid2) + hit->off)) - (int)(sr_rdseq_length(sdb, hit->rid1));
	fprintf(out, "{OVL\n");
	fprintf(out, "rds:%d,%d\n", hit->rid1 + 1, hit->rid2 + 1);
	fprintf(out, "adj:%c\n", adj);
	fprintf(out, "scr:0\n");
	fprintf(out, "ahg:%d\n", ahg);
	fprintf(out, "bhg:%d\n", bhg);
	fprintf(out, "}\n");
}

int main(int argc, char **argv){
	SR_SeqDB *sdb;
	FileReader *fr;
	Sequence *seq;
	int c, kmer_size, rd_len, min_ol, max_mm, s, gap, n_cpu, min_cpx, format;
	float min_sm;
	kmer_size = 8;
	rd_len = 0;
	min_ol = 21;
	min_sm = 0.9;
	max_mm = 5;
	s = 2;
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
	sdb = sr_init_sdb("-", n_cpu, kmer_size, rd_len);
	sr_set_align_parameters(sdb, s, min_ol, min_sm, max_mm, gap);
	sr_set_filter_parameters(sdb, min_cpx, 1024);
	if(format == 1) sr_set_output_func(sdb, output_amos_overlap_message);
	fr = stdin_filereader();
	seq = NULL;
	while(fread_fasta(&seq, fr)){
		if(seq->seq.size > MAX_RD_LEN) continue;
		sr_push_sdb(sdb, seq->seq.string, seq->seq.size, 1);
	}
	fclose_filereader(fr);
	sr_ready_sdb(sdb);
	sr_aln_sdb(sdb);
	sr_free_sdb(sdb);
	return 0;
}

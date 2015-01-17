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
#include "sr_aln.h"
#include "error_correction.h"
#include "file_reader.h"
#include "thread.h"
#include "list.h"
#include <unistd.h>

#include "memstream.h"

int usage(){
	printf(
			"Usage: sralign [Options] <fq(/a)_file1> [<fq(/a)_file2> ...]\n"
			"Options:\n"
			" -t <int>    Number of threads [1]\n"
			" -0          Fasta input [auto]\n"
			" -1          Fastq input [auto]\n"
			" -2          Output alignment [yes]\n"
			" -k <int>    Seed size, 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap [50]\n"
			" -s <float>  Min similiarity [0.9]\n"
			" -m <int>    Max mismatch [5]\n"
			" -c <int>    Max number of alignments per reads [50]\n"
			"Error Correction\n"
			" -3          Do read correction, and output corrected sequence\n"
			"             doesn't support <-g>\n"
			" -b <int>    If a base occurs more than <int>, don't correct, [4]\n"
			" -p <float>  If the frequency of a base is more than <float>, don't correct, [0.25]\n"
			"Advance options\n"
			" -r <int>    Trim all reads into <-r> bp, 0: no trim [0]\n"
			" -R <int>    Clip first <-R> bases, useful bases is <-R> ~ <-r or all> [0]\n"
			" -n <int>    Number of seeds [4]\n"
			" -f <int>    Min number of tri-kmers in one sequences, filtering low complexity [10]\n"
			" -g          Gap alignment (Smith-Waterman), [no]\n"
			" -S          Reads must be aligned on the same strand [no]\n"
			"\n"
			"Examples:\n"
			" sralign -t 8 test.fq # output test.fq.ovl\n"
			" sralign -t 8 test1.fa test2.fq # output test1.fa.ovl test2.fq.ovl\n"
			" sralign -t 16 -3 test1.fa test2.fa # output test1.fa.corr test2.fa.corr\n"
		  );
	return 1;
}

thread_beg_def(aln);
SR_SeqDB *sdb;
SR_AlnAux *aux;
int outype;
uint32_t max_cnt;
float max_freq;
FILE *out;
int task;
uint32_t *out_idx;
uint32_t batch_size;
uint32_t beg, end;
thread_end_def(aln);

thread_beg_func(aln);
SR_SeqDB *sdb;
SR_AlnAux *aux;
ReadCorrector *corr;
corr_t *cr;
uint32_t batch, beg, end, i, j, tidx, ncpu, n_rd;
FILE *cache;
size_t cache_len;
char *cache_buf;
int64_t t0, t1, t2;
sdb = aln->sdb;
aux = aln->aux;
tidx = aln->t_idx;
ncpu = aln->n_cpu;
cache_buf = NULL;
cache_len = 0;
cache = open_memstream(&cache_buf, &cache_len);
if(aln->outype){
	corr = init_corrector(sdb, aln->max_cnt, aln->max_freq);
} else {
	corr = NULL;
}
thread_beg_loop(aln);
if(aln->task == 1){
	sr_mask_low_complexity(sdb, ncpu, tidx);
} else if(aln->task == 2){
	for(i=tidx;i<sdb->n_idx;i+=ncpu){
		fprintf(stdout, "Indexing %u\n", i); fflush(stdout);
		sr_index_sdb(sdb, i);
		fprintf(stdout, "Indexed %u\n",  i); fflush(stdout);
	}
} else if(aln->task == 3){
	n_rd = aln->end - aln->beg;
	t0 = microtime();
	for(batch=tidx;batch*aln->batch_size<n_rd;batch+=ncpu){
		beg = batch * aln->batch_size;
		end = beg + aln->batch_size;
		if(end > n_rd) end = n_rd;
		t1 = microtime();
		for(i=beg+aln->beg;i<end+aln->beg;i++){
			sr_clear_aux(aux);
			sr_align_sdb(sdb, i, aux);
			if(aln->outype){
				do_corrector(corr, i, aux->hits);
				fprintf(cache, ">%u n_corr=%u", i - aln->beg, (unsigned)corr->corrs->size);
				for(j=0;j<corr->corrs->size;j++){
					cr = ref_corrv(corr->corrs, j);
					fprintf(cache, " %d:%c(%d)->%c(%d)", cr->pos, bit_base_table[cr->base[0]], cr->error_cnt, bit_base_table[cr->base[1]], cr->major_cnt);
				}
				fputc('\n', cache);
				for(j=0;j<aux->rd_len;j++){
					fputc(bit_base_table[corr->rd_bases[j]], cache);
				}
				fputc('\n', cache);
			} else {
				for(j=0;j<aux->hits->size;j++) sr_output_hit(sdb, ref_sr_hitv(aux->hits, j), cache);
			}
		}
		while(1){
			if(*aln->out_idx != tidx) microsleep(1); else break;
		}
		fflush(cache); fwrite(cache_buf, ftell(cache), 1, aln->out); fseek(cache, 0, SEEK_SET);
		t2 = microtime();
		fprintf(stdout, "[%03d] %u reads, %0.2f(%0.2f) secs. \n", tidx, end, (t2 - t0) / 1000000.0, (t2 - t1) / 1000000.0);
		fflush(stdout);
		*aln->out_idx = ((*aln->out_idx) + 1) % ncpu;
	}
	{
		while(1){
			if(*aln->out_idx != tidx) microsleep(1); else break;
		}
		fflush(cache); fwrite(cache_buf, ftell(cache), 1, aln->out); fseek(cache, 0, SEEK_SET);
		*aln->out_idx = ((*aln->out_idx) + 1) % ncpu;
	}
} else {
	microsleep(1);
}
thread_end_loop(aln);
fclose(cache);
if(cache_buf) free(cache_buf);
if(corr) free_corrector(corr);
thread_end_func(aln);

int main(int argc, char **argv){
	SR_SeqDB *sdb;
	FileReader *fr;
	FILE *out, *d_inp, *d_out;
	Sequence *seq;
	int c, kmer_size, n_seed, rd_len, rd_clip, limit, min_ol, max_mm, s, gap;
	int n_cpu, min_cpx, format, outype, max_cnt, is_fq, debug_mode;
	uint32_t out_idx, rd_num;
	u32list *rd_nums;
	float min_sm, max_freq;
	char *str;
	thread_preprocess(aln);
	limit = 50;
	kmer_size = 15;
	n_seed = 4;
	rd_len = 0;
	rd_clip = 0;
	min_ol = 50;
	min_sm = 0.9;
	max_mm = 5;
	s = 3;
	gap = 0;
	n_cpu = 1;
	min_cpx = 10;
	format = 2;
	outype = 0;
	max_cnt = 4;
	max_freq = 0.25;
	debug_mode = 0;
	while((c = getopt(argc, argv, "hg0123SdDk:r:R:l:s:f:m:n:t:c:b:p:")) != -1){
		switch(c){
			case 'k': kmer_size = atoi(optarg); break;
			case 'n': n_seed = atoi(optarg); break;
			case 'r': rd_len = atoi(optarg); break;
			case 'R': rd_clip = atoi(optarg); break;
			case 'l': min_ol = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'S': s = 1; break;
			case 'm': max_mm = atoi(optarg); break;
			case 'g': gap = 1; break;
			case '0': format = 0; break;
			case '1': format = 1; break;
			case '2': outype = 0; break;
			case '3': outype = 1; break;
			case 't': n_cpu = atoi(optarg); break;
			case 'f': min_cpx = atoi(optarg); break;
			case 'c': limit = atoi(optarg); break;
			case 'b': max_cnt = atoi(optarg); break;
			case 'p': max_freq = atof(optarg); break;
			case 'd': debug_mode = 1; break;
			case 'D': debug_mode = 2; break;
			default : return usage();
		}
	}
	if(outype == 1 && gap == 1){
		fprintf(stdout, "We don't implement consensus module for gap alignment in current version\n");
		return usage();
	}
	if(optind >= argc) return usage();
	if(debug_mode == 2){
		d_inp = fopen("_sralign_dump_file", "r");
		sdb = sr_load_sdb(d_inp);
		rd_nums = load_u32list(d_inp);
		fclose(d_inp);
	} else {
		sdb = sr_init_sdb(kmer_size, n_seed, rd_len, min_cpx);
		rd_nums = init_u32list(12);
		for(c=optind;c<argc;c++){
			if((fr = fopen_filereader(argv[c])) == NULL){
				fprintf(stdout, "Cannot open %s\n", argv[c]); fflush(stdout); abort();
			}
			if(format == 2){
				switch(guess_seq_file_type(fr)){
					case 1: is_fq = 0; break;
					case 2: is_fq = 1; break;
					default: fprintf(stdout, "Unknown file type: %s\n", argv[c]); fflush(stdout); abort();
				}
			} else is_fq = format;
			seq = NULL;
			rd_num = 0;
			while(is_fq? fread_fastq(&seq, fr) : fread_fasta(&seq, fr)){
				if(seq->seq.size <= rd_clip) continue;
				sr_push_sdb(sdb, seq->seq.string + rd_clip, seq->seq.size - rd_clip);
				rd_num ++;
			}
			fclose_filereader(fr);
			push_u32list(rd_nums, rd_num);
			fprintf(stdout, "Load %u reads from %s\n", rd_num, argv[c]); fflush(stdout);
		}
		sr_ready_sdb(sdb);
	}
	out_idx = 0;
	thread_beg_init(aln, n_cpu);
	aln->sdb = sdb;
	aln->aux = sr_init_aux(s, MAX_RD_LEN, min_ol, max_mm, min_sm, gap, 0, limit);
	sr_fit_aux2sdb(aln->aux, aln->sdb);
	aln->outype = outype;
	aln->max_cnt = max_cnt;
	aln->max_freq = max_freq;
	aln->task = 0;
	aln->out  = NULL;
	aln->batch_size = 1024;
	aln->out_idx = &out_idx;
	aln->beg = 0;
	aln->end = 0;
	thread_end_init(aln);
	if(debug_mode != 2){
		thread_apply_all(aln, aln->task = 1);
		thread_apply_all(aln, aln->task = 2);
		fprintf(stdout, " Done\n"); fflush(stdout);
	}
	if(debug_mode == 1){
		d_out = fopen("_sralign_dump_file", "w");
		sr_dump_sdb(sdb, d_out);
		dump_u32list(rd_nums, d_out);
		fclose(d_out);
	}
	for(c=0;c<(int)rd_nums->size;c++){
		if(outype){
			str = malloc(strlen(argv[c + optind]) + strlen(".corr") + 1);
			sprintf(str, "%s.corr", argv[c + optind]);
		} else {
			str = malloc(strlen(argv[c + optind]) + strlen(".ovl") + 1);
			sprintf(str, "%s.ovl", argv[c + optind]);
		}
		out = fopen(str, "w");
		fprintf(stdout, "Processing %s\n", argv[optind + c]); fflush(stdout);
		out_idx = 0;
		thread_apply_all(aln, (aln->out = out, aln->task = 3, aln->beg = aln->end, aln->end += get_u32list(rd_nums, c)));
		fprintf(stdout, "Done\n"); fflush(stdout);
		fclose(out);
		free(str);
	}
	thread_beg_close(aln);
	sr_free_aux(aln->aux);
	thread_end_close(aln);
	sr_free_sdb(sdb);
	free_u32list(rd_nums);
	return 0;
}

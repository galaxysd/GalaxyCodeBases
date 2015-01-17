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
 

#include <unistd.h>
#include <signal.h>
#include "tag.h"

static uint32_t __n_seqs = 0;
static char * progress_log = NULL;

void report_progress(int signum){
	FILE *log;
	if(progress_log == NULL){
		fprintf(stderr, "[SIG:%d] count_tag(%d) processed %u sequences\n", signum, getpid(), __n_seqs);
		fflush(stderr);
	} else {
		if((log = fopen(progress_log, "a"))){
			fprintf(log, "[SIG:%d] count_tag(%d) processed %u sequences\n", signum, getpid(), __n_seqs);
			fflush(log);
			fclose(log);
		}
	}
}

void count_tag_core(kmerhash *hash, uint8_t *seq, uint32_t seqlen, uint32_t kmer_size, uint32_t n_part, uint32_t i_part){
	Kmer KMER, *kmer;
	uint32_t i, j;
	int exists, dir;
	KMER.cnt = 0;
	for(i=0;i+kmer_size<seqlen+1;i++){
		KMER.k1 = 0;
		KMER.k2 = 0;
		dir = min_value_kmer_dir(seq + i, kmer_size);
		if(dir){
			for(j=0;j<kmer_size;j++){
				if(j < 32){
					KMER.k1 |= ((uint64_t)((~seq[seqlen-(i+j+1)])&0x03)) << (j << 1);
				} else {
					KMER.k2 |= ((uint64_t)((~seq[seqlen-(i+j+1)])&0x03)) << ((j - 32) << 1);
				}
			}
		} else {
			for(j=0;j<kmer_size;j++){
				if(j < 32){
					KMER.k1 |= ((uint64_t)((~seq[seqlen-(i+j+1)])&0x03)) << (j << 1);
				} else {
					KMER.k2 |= ((uint64_t)((~seq[seqlen-(i+j+1)])&0x03)) << ((j - 32) << 1);
				}
			}
		}
		if(n_part > 1 && kmer2code(KMER) % n_part != i_part) continue;
		kmer = prepare_kmerhash(hash, KMER, &exists);
		if(exists){
			if(kmer->cnt < MAX_TAG_CNT) kmer->cnt ++;
		} else {
			kmer->cnt = 1;
			kmer->k1 = KMER.k1;
			kmer->k2 = KMER.k2;
		}
	}
}

void count_read(kmerhash *hash, uint8_t *seq, char *chs, uint32_t seqlen, uint32_t kmer_size, uint32_t n_part, uint32_t i_part){
	uint32_t i;
	uint8_t v;
	for(i=0;i<seqlen;i++){
		v = base_bit_table[(int)chs[i]];
		if(v == 4) v = lrand48() & 0x03;
		seq[i] = v;
	}
	count_tag_core(hash, seq, seqlen, kmer_size, n_part, i_part);
}

void output_tags(kmerhash *hash, uint32_t kmer_size, uint32_t low, uint32_t high, FILE *out){
	Kmer *kmer;
	uint32_t i;
	char tag[100];
	reset_iter_kmerhash(hash);
	while((kmer = ref_iter_kmerhash(hash))){
		if(kmer->cnt < low || kmer->cnt > high) continue;
		for(i=0;i<kmer_size;i++){
			tag[i] = bit_base_table[((i < 32)? (kmer->k1 >> (i << 1)) : (kmer->k2 >> ((i - 32) << 1))) & 0x03];
		}
		tag[kmer_size] = 0;
		fprintf(out, "%s\t%u\n", tag, kmer->cnt);
	}
	fflush(out);
}

int usage(){
	printf(
"Usage: count_tag [options]\n"
"Options:\n"
" -i <string> Input file, supports multiple '-i' [NULL]\n"
" -o <string> Output file [stdout]\n"
" -f <string> FQ|FA [FQ]\n"
" -k <int>    tag size, 17 bp - %u bp, suggests odd number [%u]\n"
" -l <int>    ignore tags count less than this value, 1 - %u [1]\n"
" -h <int>    ignore tags count more than this value, 2 - %u [%u]\n"
" -x <int>    part index, for multiple-processes <-x>/<-y>, [1]\n"
" -y <int>    Number of part, for multiple-processes,       [1]\n"
"\n"
"Multiple-process:\n"
" count_tag -y 8 -x 1 ...\n"
" count_tag -y 8 -x 2 ...\n"
" ...\n"
"Progress report:\n"
" After run count_tag, use 'kill -s SIGUSR1 <pid>' to invoke count_tag to report progress\n"
"\n", MAX_KMER_SIZE, 31, MAX_TAG_CNT - 1, MAX_TAG_CNT, MAX_TAG_CNT
	);
	return 1;
}

define_list(flist, char*);

int main(int argc, char **argv){
	kmerhash *hash;
	FileReader *fr;
	Sequence *seq;
	flist *infs;
	FILE *out;
	char *ouf;
	uint32_t f, kmer_size, low, high, x, y;
	uint8_t *buf;
	int c;
	infs = init_flist(2);
	ouf  = NULL;
	f = 1;
	kmer_size = 31;
	x = y = 0;
	low = 1;
	high = MAX_TAG_CNT;
	while((c = getopt(argc, argv, "i:o:f:k:l:h:x:y:")) != -1){
		switch(c){
			case 'i': push_flist(infs, optarg); break;
			case 'o': ouf = optarg; break;
			case 'f': f = (strcasecmp(optarg, "FA") == 0)? 0 : 1; break;
			case 'k': kmer_size = atoi(optarg); break;
			case 'l': low = atoi(optarg); break;
			case 'h': high = atoi(optarg); break;
			case 'x': x = atoi(optarg); break;
			case 'y': y = atoi(optarg); break;
		}
	}
	if(count_flist(infs) == 0) return usage();
	if(y){
		if(y == 1 || x == 0 || x > y) return usage();
	} else if(x) return usage();
	else x = y = 1;
	if((fr = fopen_m_filereader(count_flist(infs), as_array_flist(infs))) == NULL){
		fprintf(stderr, " -- Cannot open input file(s) in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if(ouf == NULL) out = stdout;
	else {
		if((out = fopen(ouf, "w")) == NULL){
			fprintf(stderr, " -- Cannot write '%s' in %s -- %s:%d --\n", ouf, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		} else {
			progress_log = catstr(2, ouf, ".progress.log");
		}
	}
	hash = init_kmerhash(1023);
	seq = NULL;
	buf = malloc(MAX_RD_LEN);
	signal(SIGUSR1, report_progress);
	while(f? fread_fastq_adv(&seq, fr, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq, fr, FASTA_FLAG_NO_NAME)){
		if(seq->seq.size > MAX_RD_LEN) seq->seq.size = MAX_RD_LEN;
		count_read(hash, buf, seq->seq.string, seq->seq.size, kmer_size, y, x - 1);
		__n_seqs ++;
	}
	fclose_filereader(fr);
	output_tags(hash, kmer_size, low, high, out);
	if(ouf) fclose(out);
	free_kmerhash(hash);
	free(buf);
	if(progress_log) free(progress_log);
	return 0;
}


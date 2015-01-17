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

typedef struct {
	uint64_t k1, k2:54, cnt:10, idx:10, offset:54;
} Xmer;

define_hashset(xmerhash, Xmer, kmer_hashcode, kmer_equals);

typedef struct {
	xmerhash *hash;
	uint64_t *seqs;
	uint32_t p1len, p2len;
	uint32_t kmer_size;
	uint64_t n_seq;
	uint64_t *ids;
} SeqTagLib;

SeqTagLib* load_tags(FileReader *fr, uint32_t fix_p1_len, uint32_t fix_p2_len){
	SeqTagLib *lib;
	Xmer X;
	uint32_t i;
	char *seq;
	lib = malloc(sizeof(SeqTagLib));
	lib->p1len = fix_p1_len;
	lib->p2len = fix_p2_len;
	lib->hash = init_xmerhash(1023);
	X.offset = 0;
	X.idx = 0;
	while(fread_table(fr) != -1){
		seq = get_col_str(fr, 0);
		lib->kmer_size = get_col_len(fr, 0);
		X.cnt = atoi(get_col_str(fr, 1));
		X.k1 = 0;
		X.k2 = 0;
		for(i=0;i<lib->kmer_size;i++){
			if(i < 32){
				X.k1 |= ((uint64_t)base_bit_table[(int)seq[i]]) << (i << 1);
			} else {
				X.k2 |= ((uint64_t)base_bit_table[(int)seq[i]]) << ((i - 32) << 1);
			}
		}
		put_xmerhash(lib->hash, X);
		X.offset += X.cnt;
	}
	lib->n_seq = X.offset;
	lib->seqs  = calloc((lib->n_seq * (lib->p1len + lib->p2len) + 31) / 32, sizeof(uint64_t));
	lib->ids   = calloc(lib->n_seq, sizeof(uint64_t));
	return lib;
}

static uint64_t __n_scan = 0;
static uint64_t __n_load = 0;
static uint64_t __n_seqs = 0;
static char *progress_log = NULL;

void report_progress(int signum){
	FILE *log;
	if(progress_log == NULL){
		fprintf(stderr, "[SIG:%d] cluster_tag(%d) %llu/%llu/%llu\n", signum, getpid(), (unsigned long long)__n_load, (unsigned long long)__n_scan, (unsigned long long)__n_seqs);
		fflush(stderr);
	} else if((log = fopen(progress_log, "a"))){
		fprintf(log, "[SIG:%d] cluster_tag(%d) %llu/%llu/%llu\n", signum, getpid(), (unsigned long long)__n_load, (unsigned long long)__n_scan, (unsigned long long)__n_seqs);
		fclose(log);
	}
}

void fill_bitseqs(SeqTagLib *lib, uint64_t off, uint8_t *seq1, uint32_t len1, uint8_t *seq2, uint32_t len2){
	uint64_t a, b, j;
	for(j=0;j<len1;j++){
		a = (off + j) >> 5;
		b = ((~(off + j)) & 0x1FLLU) << 1;
		lib->seqs[a] |= ((uint64_t)seq1[j]) << b;
	}
	off += lib->p1len;
	for(j=0;j<len2;j++){
		a = (off + j) >> 5;
		b = ((~(off + j)) & 0x1FLLU) << 1;
		lib->seqs[a] |= ((uint64_t)seq2[j]) << b;
	}
}

void scan_tag_core(SeqTagLib *lib, uint64_t id, uint8_t *seq1, uint32_t len1, uint8_t *seq2, uint32_t len2, int pdir){
	Xmer X, *x;
	uint64_t off, id_flags;
	uint32_t i, j, len;
	uint8_t *seq;
	int rdir;
	X.cnt = 0;
	X.idx = 0;
	X.offset = 0;
	id_flags  = pdir? (1LLU << 63) : 0;
	seq = pdir? seq2 : seq1;
	len = pdir? len2 : len1;
	for(i=0;i+lib->kmer_size<=len;i++){
		X.k1 = 0;
		X.k2 = 0;
		rdir = min_value_kmer_dir(seq + i, lib->kmer_size);
		for(j=0;j<lib->kmer_size;j++){
			if(j < 32){
				X.k1 |= ((uint64_t)(rdir? ((~seq[len-i-j-1])&0x03) : seq[i+j])) << (j << 1);
			} else {
				X.k2 |= ((uint64_t)(rdir? ((~seq[len-i-j-1])&0x03) : seq[i+j])) << ((j - 32) << 1);
			}
		}
		x = get_xmerhash(lib->hash, X);
		if(x == NULL || x->idx == x->cnt) continue;
		lib->ids[x->offset + x->idx] = len1 | (len2 << 10) | (i << 20) | (id << 30) | id_flags | (rdir? (1LLU << 62) : 0);
		off = (x->offset + x->idx) * (lib->p1len + lib->p2len);
		fill_bitseqs(lib, off, seq1, len1, seq2, len2);
		x->idx ++;
		__n_load ++;
	}
}

void scan_tags(SeqTagLib *lib, FileReader *fwd, FileReader *rev, int is_fq, int pdirs){
	Sequence *seq1, *seq2;
	uint64_t id;
	uint32_t i, len1, len2;
	uint8_t *buf1, *buf2, v;
	__n_seqs = lib->n_seq;
	id = 0;
	seq1 = seq2 = NULL;
	buf1 = malloc(lib->p1len);
	buf2 = malloc(lib->p2len);
	while(is_fq? (fread_fastq_adv(&seq1, fwd, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) && 
		fread_fastq_adv(&seq2, rev, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL)) :
		(fread_fasta_adv(&seq1, fwd, FASTA_FLAG_NO_NAME) && fread_fasta_adv(&seq2, rev, FASTA_FLAG_NO_NAME))){
		id ++;
		__n_scan ++;
		len1 = seq1->seq.size;
		if(len1 > lib->p1len) len1 = lib->p1len;
		for(i=0;i<len1;i++){
			v = base_bit_table[(int)seq1->seq.string[i]];
			if(v == 4) v = lrand48() & 0x03;
			buf1[i] = v;
		}
		len2 = seq2->seq.size;
		if(len2 > lib->p2len) len2 = lib->p2len;
		for(i=0;i<len2;i++){
			v = base_bit_table[(int)seq2->seq.string[i]];
			if(v == 4) v = lrand48() & 0x03;
			buf2[i] = v;
		}
		if(pdirs & 0x01){ scan_tag_core(lib, id, buf1, len1, buf2, len2, 0); }
		if(pdirs & 0x02){ scan_tag_core(lib, id, buf1, len1, buf2, len2, 1); }
	}
	free(buf1);
	free(buf2);
	if(seq1) free_sequence(seq1);
	if(seq2) free_sequence(seq2);
}

void output_tags(SeqTagLib *lib, FILE *out){
	Xmer *x;
	char *tag, *seq1, *seq2;
	uint64_t off, id;
	uint32_t i, len1, len2;
	tag = malloc(lib->kmer_size + 1);
	seq1 = malloc(lib->p1len + 1);
	seq2 = malloc(lib->p2len + 1);
	reset_iter_xmerhash(lib->hash);
	while((x = ref_iter_xmerhash(lib->hash))){
		if(x->idx == 0) continue;
		for(i=0;i<lib->kmer_size;i++){
			if(i < 32){
				tag[i] = bit_base_table[(x->k1 >> (i << 1)) & 0x03];
			} else {
				tag[i] = bit_base_table[(x->k2 >> ((i - 32) << 1)) & 0x03];
			}
		}
		tag[lib->kmer_size] = 0;
		for(i=0;i<x->idx;i++){
			id = lib->ids[x->offset + i];
			len1 = id & 0x3FFU;
			len2 = (id >> 10) & 0x3FFU;
			off = (x->offset + i) * (lib->p1len + lib->p2len);
			bits2seq(seq1, lib->seqs, off, len1);
			bits2seq(seq2, lib->seqs, off + lib->p1len, len2);
			fprintf(out, "%s\t%u\t%c%c\t%u\t%s\t%s\n", tag, (uint32_t)(id >> 30), "FC"[(id >> 63)&0x01], "FC"[(id>>62)&0x01], (uint32_t)(id >> 20) & 0x3FFU, seq1, seq2);
		}
	}
	free(tag);
	free(seq1);
	free(seq2);
}

int usage(){
	printf(
"Usage: cluster_tag [options]\n"
"Options:\n"
" -i <string> Tag file [stdin]\n"
" -a <string> Input file, pair1 [NULL]\n"
" -b <string> Input file, pair2 [NULL]\n"
" -f <string> FQ | FA [FQ]\n"
" -1 <int>    Max read length in pair1 [100]\n"
" -2 <int>    Max read length in pair2 [100]\n"
" -o <string> Output file [stdout]\n"
" -p <int>    Flags of pair, 1: pair1, 2: pair2, 3: both [3]\n"
	);
	return 1;
}

define_list(fslist, char*);

int main(int argc, char **argv){
	SeqTagLib *lib;
	FileReader *tag, *fwd, *rev;
	FILE *out;
	fslist *fa, *fb;
	char *inf, *ouf;
	int is_fq, pdirs, p1, p2, c;
	fa = init_fslist(2);
	fb = init_fslist(2);
	inf = NULL;
	ouf = NULL;
	is_fq = 1;
	pdirs = 3;
	p1 = p2 = 100;
	while((c = getopt(argc, argv, "i:a:b:1:2:f:o:p:")) != -1){
		switch(c){
			case 'i': inf = optarg; break;
			case 'a': push_fslist(fa, optarg); break;
			case 'b': push_fslist(fb, optarg); break;
			case 'o': ouf = optarg; break;
			case '1': p1 = atoi(optarg); break;
			case '2': p2 = atoi(optarg); break;
			case 'f': is_fq = (strcasecmp(optarg, "FA") == 0)? 0 : 1; break;
			case 'p': pdirs = atoi(optarg); break;
		}
	}
	if(count_fslist(fa) == 0 || count_fslist(fa) != count_fslist(fb)){
		free_fslist(fa);
		free_fslist(fb);
		return usage();
	}
	if(inf == NULL) tag = stdin_filereader();
	else if((tag = fopen_filereader(inf)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", inf, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if(ouf == NULL) out = stdout;
	else {
		if((out = fopen(ouf, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", ouf, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		} else progress_log = catstr(2, ouf, ".pregress.log");
	}
	if((fwd = fopen_m_filereader(count_fslist(fa), as_array_fslist(fa))) == NULL){
		fprintf(stderr, " -- Cannot open pair files in '-a' in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	if((rev = fopen_m_filereader(count_fslist(fb), as_array_fslist(fb))) == NULL){
		fprintf(stderr, " -- Cannot open pair files in '-b' in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	free_fslist(fa);
	free_fslist(fb);
	signal(SIGUSR1, report_progress);
	lib = load_tags(tag, p1, p2);
	fclose_filereader(tag);
	scan_tags(lib, fwd, rev, is_fq, pdirs);
	fclose_filereader(fwd);
	fclose_filereader(rev);
	output_tags(lib, out);
	if(ouf){ fclose(out); free(progress_log); }
	return 0;
}

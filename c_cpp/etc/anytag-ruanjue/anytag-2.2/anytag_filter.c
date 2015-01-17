
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
#include <unistd.h>
#include "anytag_aux.h"
#include "sr_aln.h"
#include "file_reader.h"
#include "thread.h"

const char *version = "1.0";

int validate_pseudo_sanger(SR_SeqDB *sdb, SR_AlnAux *aux, uint32_t seqid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen, BitVec *flags[3], char *str[2], FILE *log){
	SR_AlnHit *hit, *hit1, *hit2;
	uint32_t i, j, s, last_rid, b1, b2, e1, e2;
	clear_sr_hitv(aux->hits);
	sr_align_long_sdb(sdb, aux, seqid, seqs, seqoff, seqlen, 3, 0);
	if(log){
		bits2seq(str[0], seqs, seqoff, seqlen);
		fprintf(log, "+%010u\t%s\n", seqid, str[0]);
	}
	if(aux->hits->size == 0){
		return 0;
	}
	last_rid = ref_sr_hitv(aux->hits, 0)->rid2;
	zeros_bitvec(flags[0]);
	//pairing
	for(i=1;i<aux->hits->size;i++){
		hit = ref_sr_hitv(aux->hits, i);
		if(last_rid & 0x01){
			if(hit->rid2 == (last_rid & (~0x1U))){
				one_bitvec(flags[0], i - 1);
				one_bitvec(flags[0], i);
			}
		}
		last_rid = hit->rid2;
	}
	if(log){
		for(i=0;i<aux->hits->size;i++){
			if(get_bitvec(flags[0], i) == 0) continue;
			hit = ref_sr_hitv(aux->hits, i);
			if(hit->dir1) reverse_cigars(hit->cigars, hit->n_cigar);
			if(hit->dir1 ^ hit->dir2){
				bits2revseq(str[0], ref_u64list(sdb->rd_seqs, 0), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
			} else {
				bits2seq(str[0], ref_u64list(sdb->rd_seqs, 0), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
			}
			fprintf(log, "-%010u\t", hit->rid2);
			cigars_seq2aln(str[1], hit->cigars, hit->n_cigar, 1, str[0]);
			fprintf(log, "%s\t", str[1]);
			cigars2string(hit->cigars, hit->n_cigar, str[1]);
			fprintf(log, "%c%c\t%s\n", "+-"[hit->dir1], "+-"[hit->dir2], str[1]);
		}
	}
	if(log) zeros_bitvec(flags[1]);
	zeros_bitvec(flags[2]);
	j = 0xFFFFFFFFU;
	for(i=0;i<aux->hits->size;i++){
		if(get_bitvec(flags[0], i) == 0) continue;
		if(j == 0xFFFFFFFFU){ j = i; continue; }
		hit1 = ref_sr_hitv(aux->hits, j);
		hit2 = ref_sr_hitv(aux->hits, i);
		j = 0xFFFFFFFFU;
		if(hit1->dir1){ e1 = seqlen - hit1->off; b1 = e1 - sr_rdseq_length(sdb, hit1->rid2); }
		else { b1 = hit1->off; e1 = b1 + sr_rdseq_length(sdb, hit1->rid2); }
		if(hit2->dir1){ e2 = seqlen - hit2->off; b2 = e2 - sr_rdseq_length(sdb, hit2->rid2); }
		else { b2 = hit2->off; e2 = b2 + sr_rdseq_length(sdb, hit2->rid2); }
		if(log){
			for(s=b1;s<e1;s++) one_bitvec(flags[1], s);
			for(s=b2;s<e2;s++) one_bitvec(flags[1], s);
		}
		if(b1 < b2){
			for(s=b1;s<e2;s++) one_bitvec(flags[2], s);
		} else {
			for(s=b2;s<e1;s++) one_bitvec(flags[2], s);
		}
	}
	if(log){
		fprintf(log, "base coverage\t");
		for(i=0;i<seqlen;i++) fprintf(log, "%d", (int)get_bitvec(flags[1], i));
		fprintf(log, "\n");
		fprintf(log, "pair coverage\t");
		for(i=0;i<seqlen;i++) fprintf(log, "%d", (int)get_bitvec(flags[2], i));
		fprintf(log, "\n");
		fprintf(log, "//\n");
	}
	for(i=sdb->max_rd_len;i+sdb->max_rd_len<seqlen;i++){
		if(get_bitvec(flags[2], i) == 0) return 0;
	}
	return 1;
}

thread_beg_def(mval);
SR_SeqDB *sdb;
uint32_t n_seq;
uint32_t idoff;
u64list  *seqs;
u64list  *rd_offs;
u32list  *rd_lens;
String   *text;
u32list  *text_offs;
char     *cache[2];
size_t   size[2];
FILE     *out[2];
int      verbose;
int      task;
thread_end_def(mval);

thread_beg_func(mval);
uint32_t i, seqlen;
uint64_t seqoff;
SR_SeqDB *sdb;
SR_AlnAux *aux;
BitVec *flags[3];
char *str[2];
sdb = mval->sdb;
mval->cache[0] = NULL;
mval->cache[1] = NULL;
mval->size[0]  = 0;
mval->size[1]  = 0;
mval->out[0]   = open_memstream(&mval->cache[0], &mval->size[0]);
mval->out[1]   = open_memstream(&mval->cache[1], &mval->size[1]);
str[0]         = malloc(MAX_RD_LEN * 2);
str[1]         = malloc(MAX_RD_LEN * 2);
flags[0]       = init_bitvec(1024 * 1024);
flags[1]       = init_bitvec(1024 * 1024);
flags[2]       = init_bitvec(1024 * 1024);
aux = sr_init_aux();
thread_beg_loop(mval);
if(mval->task == 1){
	for(i=mval->t_idx;i<sdb->n_idx;i+=mval->n_cpu) sr_index_sdb(sdb, i);
	break;
}
for(i=0;i<mval->n_seq;i++){
	seqoff = get_u64list(mval->rd_offs, i);
	seqlen = get_u32list(mval->rd_lens, i);
	if(validate_pseudo_sanger(sdb, aux, mval->idoff + i, ref_u64list(mval->seqs, 0), seqoff, seqlen, flags, str, mval->verbose? mval->out[1] : NULL)){
		fprintf(mval->out[0], "%s", mval->text->string + get_u32list(mval->text_offs, i));
	}
}
thread_end_loop(mval);
sr_free_aux(aux);
fclose(mval->out[0]);
fclose(mval->out[1]);
if(mval->cache[0]) free(mval->cache[0]);
if(mval->cache[1]) free(mval->cache[1]);
free(str[0]);
free(str[1]);
free_bitvec(flags[0]);
free_bitvec(flags[1]);
free_bitvec(flags[2]);
thread_end_func(mval);

int usage(){
	fprintf(stdout,
"ps_filter -- A tool for validating pseudo-sanger sequence\n"
"Version: %s\n"
"Author: Jue Ruan <ruanjue@gmail.com> and Zechen Chong <chongzechen@gmail.com>\n"
"Description: First, ps_filter map paired-end reads with very small insert sizes (such as 200bp) to pseudo-sanger reads\n"
"             If a pseudo-sanger sequence has regions (excepts two ends) cannot be covered by paired-end inserts, it will\n"
"             be filtered out. Please use as many small insert paired-end reads as you can, but DON't use libraries with \n"
"             insert size bigger than half of pseudo-sanger sequences' mean size\n"
"Usage: ps_filter [options] <pseudo_sanger_file> <XXX.insYYY[.varZZZ].s1.fq/fa[.gz]> <XXX.insYYY[.varZZZ].s2.fq/fa[.gz]> ...\n"
"Output: stdout\n"
"Options:\n"
" -t <int>    Number of threads [8]\n"
" -k <int>    Kmer size [15]\n"
" -l <int>    Min overlap [50]\n"
" -s <float>  Min similarity [0.95]\n"
" -m <int>    Max mismatches [8]\n"
" -v          Print detailed into stderr\n"
"\n", version
);
	return 1;
}

int main(int argc, char **argv){
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	FileReader *s1, *s2;
	Sequence *seq1, *seq2;
	char *inf;
	uint64_t seqoff;
	uint32_t seqid, seqlen, i, idx, size;
	uint32_t ksize, nseed, ncpu, min_ovl, max_mm, verbose, block;
	float min_sim;
	int ii, c;
	thread_preprocess(mval);
	ncpu  = 8;
	ksize = 15;
	nseed = 4;
	min_ovl = 50;
	min_sim = 0.95;
	max_mm  = 8;
	verbose = 0;
	block = 100;
	while((c = getopt(argc, argv, "ht:k:n:l:s:m:v")) != -1){
		switch(c){
			case 't': ncpu = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'l': min_ovl = atoi(optarg); break;
			case 's': min_sim = atof(optarg); break;
			case 'm': max_mm = atoi(optarg); break;
			case 'v': verbose = 1; break;
			default: return usage();
		}
	}
	c = optind;
	if(argc < 3 + c) return usage();
	inf = argv[c];
	sfs = init_sfv(4);
	for(ii=c+1;ii<argc;ii++){
		sf1 = next_ref_sfv(sfs);
		parse_seqfile(sf1, argv[ii], 50);
	}
	qsort(ref_sfv(sfs, 0), count_sfv(sfs), sizeof(SeqFile), cmp_seqfile);
	for(i=0;2*i<count_sfv(sfs);i++){
		if(strcasecmp(ref_sfv(sfs, 2 * i)->prefix, ref_sfv(sfs, 2 * i + 1)->prefix)
			|| ref_sfv(sfs, 2 * i)->ins_size != ref_sfv(sfs, 2 * i + 1)->ins_size
			|| ref_sfv(sfs, 2 * i)->pair_idx >= ref_sfv(sfs, 2 * i + 1)->pair_idx){
			fprintf(stderr, " -- Wrong pair \"%s\" and \"%s\" in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	sdb = sr_init_sdb(ksize, nseed, 0);
	sr_set_align_parameters(sdb, 3, min_ovl, min_sim, max_mm, 0);
	sr_set_filter_parameters(sdb, 10, 4086, 0);
	for(i=0;2*i<count_sfv(sfs);i++){
		seq1 = seq2 = NULL;
		sf1 = ref_sfv(sfs, 2 * i);
		sf2 = ref_sfv(sfs, 2 * i + 1);
		if((s1 = fopen_filereader(sf1->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf1->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		if((s2 = fopen_filereader(sf2->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf2->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		while(1){
			if(!((sf1->is_fq)? fread_fastq_adv(&seq1, s1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq1, s1, FASTA_FLAG_NO_NAME))) break;
			if(!((sf2->is_fq)? fread_fastq_adv(&seq2, s2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq2, s2, FASTA_FLAG_NO_NAME))) break;
			sr_push_sdb(sdb, seq1->seq.string, seq1->seq.size);
			sr_push_sdb(sdb, seq2->seq.string, seq2->seq.size);
		}
		if(seq1) free_sequence(seq1);
		if(seq2) free_sequence(seq2);
		fclose_filereader(s1);
		fclose_filereader(s2);
	}
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	sr_ready_sdb(sdb);
	thread_beg_init(mval, (int)ncpu);
	mval->sdb = sdb;
	mval->verbose = verbose;
	mval->n_seq = 0;
	mval->seqs  = init_u64list(1024);
	mval->rd_offs = init_u64list(block);
	mval->rd_lens = init_u32list(block);
	mval->text = init_string(1024);
	mval->text_offs = init_u32list(block);
	mval->task = 0;
	thread_end_init(mval);
	thread_beg_iter(mval);
	mval->task = 1;
	thread_wake(mval);
	thread_end_iter(mval);
	thread_waitfor_all_idle(mval);
	if((s1 = fopen_filereader(inf)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", inf, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	seq1 = NULL;
	seqid = 0;
	idx = 0;
	while(1){
		thread_beg_operate(mval, idx);
		{
			thread_waitfor_idle(mval);
			fflush(mval->out[0]);
			size = ftell(mval->out[0]);
			if(size){ fwrite(mval->cache[0], size, 1, stdout); fflush(stdout); }
			fseek(mval->out[0], 0, SEEK_SET);
			size = ftell(mval->out[1]);
			fflush(mval->out[1]);
			if(size){ fwrite(mval->cache[1], size, 1, stderr); fflush(stderr); }
			fseek(mval->out[1], 0, SEEK_SET);
		}
		clear_u64list(mval->rd_offs);
		clear_u32list(mval->rd_lens);
		clear_string(mval->text);
		clear_u32list(mval->text_offs);
		seqoff = 0;
		mval->idoff = seqid;
		for(mval->n_seq=0;mval->n_seq<block;mval->n_seq++){
			if(fread_fasta(&seq1, s1) == 0) break;
			seqid ++;
			seqlen = seq1->seq.size;
			encap_u64list(mval->seqs, (seqoff + seq1->seq.size + 31) / 32);
			seq2bits(ref_u64list(mval->seqs, 0), seqoff, seq1->seq.string, seqlen);
			push_u64list(mval->rd_offs, seqoff);
			push_u32list(mval->rd_lens, seqlen);
			seqoff += seqlen;
			push_u32list(mval->text_offs, mval->text->size);
			add_char_string(mval->text, '>');
			append_string(mval->text, seq1->name.string, seq1->name.size);
			if(seq1->comment.size){
				add_char_string(mval->text, ' ');
				append_string(mval->text, seq1->comment.string, seq1->comment.size);
			}
			add_char_string(mval->text, '\n');
			for(i=0;i+100<=seqlen;i+=100){
				append_string(mval->text, seq1->seq.string + i, 100);
				add_char_string(mval->text, '\n');
			}
			if(i < seqlen){
				append_string(mval->text, seq1->seq.string + i, seqlen - i);
				add_char_string(mval->text, '\n');
			}
			add_char_string(mval->text, '\0');
		}
		mval->task = 2;
		thread_wake(mval);
		idx = (idx + 1) % ncpu;
		if(mval->n_seq < block) break;
	}
	fclose_filereader(s1);
	for(i=0;i<ncpu;i++){
		thread_beg_operate(mval, ((idx + i) % ncpu));
		{
			thread_waitfor_idle(mval);
			size = ftell(mval->out[0]);
			if(size){ fwrite(mval->cache[0], size, 1, stdout); fflush(stdout); }
			fseek(mval->out[0], 0, SEEK_SET);
			size = ftell(mval->out[1]);
			if(size){ fwrite(mval->cache[1], size, 1, stderr); fflush(stderr); }
			fseek(mval->out[1], 0, SEEK_SET);
		}
	}
	thread_beg_close(mval);
	free_u64list(mval->seqs);
	free_u64list(mval->rd_offs);
	free_u32list(mval->rd_lens);
	free_string(mval->text);
	free_u32list(mval->text_offs);
	thread_end_close(mval);
	sr_free_sdb(sdb);
	return 0;
}

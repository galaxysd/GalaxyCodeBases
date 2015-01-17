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
#include "file_reader.h"
#include "sr_aln.h"
#include "thread.h"
#include "anytag_aux.h"

int usage_aln();

thread_begin_def(maln);
SR_SeqDB *sdb;
uint32_t rid_beg, rid_end;
int task;
FILE *out;
thread_end_def(maln);

thread_begin_func(maln);
SR_SeqDB *sdb;
SR_AlnAux *aux;
SR_AlnHit *hit;
char seq[256];
FILE *cache;
char *buf;
size_t buf_size;
uint64_t rd_off;
uint32_t i, j, rd_len, t_idx, n_cpu, l, r, dir;
sdb = maln->sdb;
t_idx = maln->t_idx;
n_cpu = maln->n_cpu;
aux = sr_init_aux();
buf = NULL;
buf_size = 0;
cache = open_memstream(&buf, &buf_size);
thread_begin_loop(maln);
if(maln->task == 1){ // Index
	for(i=t_idx;i<sdb->n_idx;i+=n_cpu){
		sr_index_sdb(sdb, i);
	}
} else if(maln->task == 2){ // Align
	for(i=maln->rid_beg+t_idx;i<maln->rid_end;i+=n_cpu){
		clear_sr_hitv(aux->hits);
		sr_align_sdb(sdb, 2 * i, aux);
		l = count_sr_hitv(aux->hits);
		sr_align_sdb(sdb, 2 * i + 1, aux);
		r = count_sr_hitv(aux->hits) - l;
		if(count_sr_hitv(aux->hits) == 0) continue;
		rd_off = sr_rdseq_offset(sdb, 2 * i);
		rd_len = sr_rdseq_length(sdb, 2 * i);
		bits2seq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
		fprintf(cache, "T\t%u\t+\t%s", i, seq);
		rd_off = sr_rdseq_offset(sdb, 2 * i + 1);
		rd_len = sr_rdseq_length(sdb, 2 * i + 1);
		bits2revseq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
		fprintf(cache, "\t%s\t%u\t%u\n", seq, l, r);
		for(j=0;j<count_sr_hitv(aux->hits);j++){
			hit = ref_sr_hitv(aux->hits, j);
			if(hit->rid2 == aux->rid) continue;
			dir = hit->dir1 ^ hit->dir2;
			dir ^= hit->rid1 & 0x01;
			dir ^= hit->rid2 & 0x01;
			fprintf(cache, "S\t%u\t%c", hit->rid2, "+-"[dir]);
			if(dir){
				rd_off = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 1);
				rd_len = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 1);
				bits2seq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
				fprintf(cache, "\t%s", seq);
				rd_off = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 0);
				rd_len = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 0);
				bits2revseq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
				fprintf(cache, "\t%s", seq);
			} else {
				rd_off = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 0);
				rd_len = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 0);
				bits2seq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
				fprintf(cache, "\t%s", seq);
				rd_off = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 1);
				rd_len = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 1);
				bits2revseq(seq, ref_u64list(sdb->rd_seqs, 0), rd_off, rd_len);
				fprintf(cache, "\t%s", seq);
			}
			fprintf(cache, "\t%u\n", hit->off);
			fflush(cache);
		}
		if(ftell(cache) >= 16 * 1024 * 1024){
			thread_begin_syn(maln);
			fwrite(buf, ftell(cache), 1, maln->out);
			fseek(cache, 0, SEEK_SET);
			thread_end_syn(maln);
		}
	}
	thread_begin_syn(maln);
	fwrite(buf, ftell(cache), 1, maln->out);
	fseek(cache, 0, SEEK_SET);
	thread_end_syn(maln);
} else {
	microsleep(1);
}
thread_end_loop(maln);
sr_free_aux(aux);
fclose(cache);
if(buf) free(buf);
thread_end_func(maln);

int main_aln(int argc, char **argv){
	ATOptions opt;
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	LibInserts *libs;
	FileReader *s1, *s2;
	Sequence *seq1, *seq2;
	FILE *log, *out;
	char *prefix, PREFIX[256], *date;
	uint32_t i, rid, rrid, n_vis;
	int ii, c, vis, is_gz;
	time_t tm;
	thread_preprocess(maln);
	if((c = parse_options(&opt, 0, argc, argv)) == -1) return usage_aln();
	if(c == argc || (argc - c) % 2 == 1) return usage_aln();
	prefix = opt.outf;
	sfs = init_sfv(4);
	for(ii=c;ii<argc;ii++){
		sf1 = next_ref_sfv(sfs);
		parse_seqfile(sf1, argv[ii]);
	}
	qsort(ref_sfv(sfs, 0), count_sfv(sfs), sizeof(SeqFile), cmp_seqfile);
	if(opt.min_lib == 0) opt.min_lib = ref_sfv(sfs, 0)->ins_size;
	for(i=0;2*i<count_sfv(sfs);i++){
		if(strcasecmp(ref_sfv(sfs, 2 * i)->prefix, ref_sfv(sfs, 2 * i + 1)->prefix)
			|| ref_sfv(sfs, 2 * i)->ins_size != ref_sfv(sfs, 2 * i + 1)->ins_size
			|| ref_sfv(sfs, 2 * i)->pair_idx >= ref_sfv(sfs, 2 * i + 1)->pair_idx){
			fprintf(stderr, " -- Wrong pair \"%s\" and \"%s\" in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	libs = guess_lib_inserts(sfs, opt.min_ins, opt.max_ins);
	sdb = sr_init_sdb(opt.kmer_size[0], opt.n_seed[0], opt.rd_len);
	sr_set_align_parameters(sdb, 1, opt.min_ol[0], opt.min_sm[0], opt.max_mm[0], 0);
	sr_set_filter_parameters(sdb, opt.low_cpx, opt.limit, 0);
	is_gz = 0;
	if(prefix == NULL){
		log = stderr;
		out = stdout;
	} else {
		sprintf(PREFIX, "%s.info", prefix);
		if((log = fopen(PREFIX, "w")) == NULL){
			fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", PREFIX, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			abort();
		}
		if(strlen(prefix) > 3 && strcmp(prefix + strlen(prefix) - 3, ".gz") == 0){
			is_gz = 1;
			sprintf(PREFIX, "gzip -c >%s.gz", prefix);
			if((out = popen(PREFIX, "w")) == NULL){
				fprintf(stderr, " -- Cannot invoke '%s' in %s -- %s:%d --\n", PREFIX, __FUNCTION__, __FILE__, __LINE__);
				fflush(stderr);
				abort();
			}
		} else {
			if((out = fopen(prefix, "w")) == NULL){
				fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", prefix, __FUNCTION__, __FILE__, __LINE__);
				fflush(stderr);
				abort();
			}
		}
	}
	rid = 0;
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	n_vis = 0;
	for(i=0;2*i<count_sfv(sfs);i++){
		set_u32list(libs->lib_offs, i, rid);
		seq1 = seq2 = NULL;
		sf1 = ref_sfv(sfs, 2 * i);
		sf2 = ref_sfv(sfs, 2 * i + 1);
		vis = (sf1->ins_size >= (int)opt.min_lib);
		if((s1 = fopen_filereader(sf1->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf1->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		if((s2 = fopen_filereader(sf2->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf2->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		rrid = rid;
		while(1){
			if(!((sf1->is_fq)? fread_fastq_adv(&seq1, s1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq1, s1, FASTA_FLAG_NO_NAME))) break;
			if(!((sf2->is_fq)? fread_fastq_adv(&seq2, s2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq2, s2, FASTA_FLAG_NO_NAME))) break;
			sr_push_sdb(sdb, seq1->seq.string, seq1->seq.size);
			sr_push_sdb(sdb, seq2->seq.string, seq2->seq.size);
			rid ++;
		}
		if(seq1) free_sequence(seq1);
		if(seq2) free_sequence(seq2);
		fclose_filereader(s1);
		fclose_filereader(s2);
		fprintf(log, "%s\t%d\t%d\t%d\t%s\n", sf1->prefix, sf1->ins_size, rrid, rid, vis? "yes" : "no");
		fflush(log);
		if(vis) n_vis = rid;
	}
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	sr_ready_sdb(sdb);
	thread_begin_init(maln, (int)opt.n_cpu);
	maln->sdb     = sdb;
	maln->rid_beg = 0;
	maln->rid_end = n_vis;
	maln->task    = 0;
	maln->out     = out;
	thread_end_init(maln);
	thread_begin_iter(maln);
	maln->task = 1; // indexing
	thread_wake(maln);
	thread_end_iter(maln);
	thread_waitfor_all_idle(maln);
	thread_begin_iter(maln);
	maln->task = 2; // alinging
	thread_wake(maln);
	thread_end_iter(maln);
	thread_waitfor_all_idle(maln);
	thread_begin_close(maln);
	thread_end_close(maln);
	sr_free_sdb(sdb);
	free_lib_inserts(libs);
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	if(log != stderr) fclose(log);
	if(is_gz) pclose(out);
	else if(out != stdout) fclose(out);
	return 0;
}


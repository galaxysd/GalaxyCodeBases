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
#include "sr_aln.h"
#include "local_assembly.h"
#include "file_reader.h"
#include "thread.h"
#include "anytag_aux.h"

int usage_all();


thread_beg_def(mall);
SR_SeqDB *sdb;
uint32_t rid_beg, rid_end;
uint32_t lmt;
LGraph *g;
int task;
FILE *cnsf, *samf, *msaf, *log;
thread_end_def(mall);

thread_beg_func(mall);
SR_SeqDB *sdb;
SR_AlnAux *aux1, *aux2;
SR_AlnHit *hit;
LGraph *g;
char seq1[256], seq2[256];
FILE *cnsf, *samf, *msaf, *log;
char *cns_buf, *sam_buf, *msa_buf;
size_t cns_size, sam_size, msa_size;
uint64_t rd_off1, rd_off2;
uint32_t i, j, rd_len1, rd_len2, t_idx, n_cpu, l, r, dir, n, m;
time_t t1, t2;
sdb = mall->sdb;
g = mall->g;
t_idx = mall->t_idx;
n_cpu = mall->n_cpu;
aux1 = sr_init_aux();
aux2 = sr_init_aux();
cns_buf = NULL;
sam_buf = NULL;
msa_buf = NULL;
cns_size = 0;
sam_size = 0;
msa_size = 0;
cnsf = open_memstream(&cns_buf, &cns_size);
samf = open_memstream(&sam_buf, &sam_size);
msaf = open_memstream(&msa_buf, &msa_size);
log  = mall->log;
thread_beg_loop(mall);
if(mall->task == 1){
	for(i=t_idx;i<sdb->n_idx;i+=n_cpu){
		fprintf(log, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx);
		fflush(log);
		sr_index_sdb(sdb, i);
		fprintf(log, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx);
		fflush(log);
	}
} else if(mall->task == 2){
	t1 = time(NULL);
	n = m = 0;
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
			fflush(log);
		}
		clear_sr_hitv(aux1->hits);
		sr_align_sdb(sdb, 2 * i, aux1);
		clear_sr_hitv(aux2->hits);
		sr_align_sdb(sdb, 2 * i + 1, aux2);
		l = count_sr_hitv(aux1->hits);
		r = count_sr_hitv(aux2->hits);
		if((l + r) == 0) continue;
		if(l + r > mall->lmt){
			if(l > mall->lmt / 2){
				if(r > mall->lmt / 2){
					l = mall->lmt / 2;
					r = mall->lmt / 2;
				} else {
					l = mall->lmt - r;
				}
			} else {
				r = mall->lmt - l;
			}
		}
		reset_lgraph(g);
		rd_off1 = sr_rdseq_offset(sdb, 2 * i);
		rd_off2 = sr_rdseq_offset(sdb, 2 * i + 1);
		rd_len1 = sr_rdseq_length(sdb, 2 * i);
		rd_len2 = sr_rdseq_length(sdb, 2 * i + 1);
		bits2seq(seq1, ref_u64list(sdb->rd_seqs, 0), rd_off1, rd_len1);
		bits2revseq(seq2, ref_u64list(sdb->rd_seqs, 0), rd_off2, rd_len2);
		push_lgraph(g, 2 * i + 0, seq1, rd_len1, 0);
		push_lgraph(g, 2 * i + 1, seq2, rd_len2, 1);
		for(j=0;j<l+r;j++){
			if(j < l){
				hit = ref_sr_hitv(aux1->hits, j);
				if(hit->rid2 == aux1->rid) continue;
			} else {
				hit = ref_sr_hitv(aux2->hits, j - l);
				if(hit->rid2 == aux2->rid) continue;
			}
			dir = hit->dir1 ^ hit->dir2;
			dir ^= hit->rid1 & 0x01;
			dir ^= hit->rid2 & 0x01;
			if(dir){
				rd_off1 = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 1);
				rd_off2 = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 0);
				rd_len1 = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 1);
				rd_len2 = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 0);
				bits2seq(seq1, ref_u64list(sdb->rd_seqs, 0), rd_off1, rd_len1);
				bits2revseq(seq2, ref_u64list(sdb->rd_seqs, 0), rd_off2, rd_len2);
				push_lgraph(g, 2 * hit->rid2 * 2 + 1, seq1, rd_len1, 0);
				push_lgraph(g, 2 * hit->rid2 * 2 + 0, seq2, rd_len2, 1);
			} else {
				rd_off1 = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 0);
				rd_off2 = sr_rdseq_offset(sdb, 2 * (hit->rid2 >> 1) + 1);
				rd_len1 = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 0);
				rd_len2 = sr_rdseq_length(sdb, 2 * (hit->rid2 >> 1) + 1);
				bits2seq(seq1, ref_u64list(sdb->rd_seqs, 0), rd_off1, rd_len1);
				bits2revseq(seq2, ref_u64list(sdb->rd_seqs, 0), rd_off2, rd_len2);
				push_lgraph(g, 2 * hit->rid2 * 2 + 0, seq1, rd_len1, 0);
				push_lgraph(g, 2 * hit->rid2 * 2 + 1, seq2, rd_len2, 1);
			}
		}
		align_lgraph(g);
		if(layout_lgraph(g)){
			consensus_lgraph(g);
			realign_lgraph(g);
			output_lgraph(g, cnsf, samf, msaf);
			n ++;
			if(n && (n % 1000) == 0){
				thread_beg_syn(mall);
				fflush(cnsf); fflush(samf); fflush(msaf);
				fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
				fwrite(sam_buf, ftell(samf), 1, mall->samf);
				fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
				fseek(cnsf, 0, SEEK_SET); fseek(samf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
				thread_end_syn(mall);
			}
		}
	}
	t2 = time(NULL);
	fprintf(log, "[THREAD%03u] %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
	fflush(log);
	if(n % 1000){
		thread_beg_syn(mall);
		fflush(cnsf); fflush(samf); fflush(msaf);
		fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
		fwrite(sam_buf, ftell(samf), 1, mall->samf);
		fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
		fseek(cnsf, 0, SEEK_SET); fseek(samf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
		thread_end_syn(mall);
	}
} else {
	microsleep(1);
}
thread_end_loop(mall);
sr_free_aux(aux1);
sr_free_aux(aux2);
fclose(cnsf); fclose(samf); fclose(msaf);
if(cns_buf) free(cns_buf);
if(sam_buf) free(sam_buf);
if(msa_buf) free(msa_buf);
thread_end_func(mall);

int main_all(int argc, char **argv){
	ATOptions opt;
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	FileReader *s1, *s2;
	Sequence *seq1, *seq2;
	FILE *log, *cnsf, *samf, *msaf;
	char *prefix, cmd[256], *date;
	uint32_t i, rid, rrid, n_vis;
	int ii, c, vis;
	time_t tm;
	thread_preprocess(mall);
	if((c = parse_options(&opt, 0, argc, argv)) == -1) return usage_all();
	if(c == argc || (argc - c) % 2 != 1) return usage_all();
	prefix = argv[c];
	sfs = init_sfv(4);
	for(ii=c+1;ii<argc;ii++){
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
	sdb = sr_init_sdb(opt.kmer_size[0], opt.n_seed[0], opt.rd_len);
	sr_set_align_parameters(sdb, 1, opt.min_ol[0], opt.min_sm[0], opt.max_mm[0], 0);
	sr_set_filter_parameters(sdb, opt.low_cpx, opt.limit, 0);
	sprintf(cmd, "%s.info", prefix);
	if((log = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
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
	rid = 0;
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	n_vis = 0;
	for(i=0;2*i<count_sfv(sfs);i++){
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
	thread_beg_init(mall, (int)opt.n_cpu);
	mall->sdb = sdb;
	mall->g = init_lgraph(opt.kmer_size[1], opt.rd_len, opt.min_ol[1], opt.min_sm[1], opt.max_mm[1], opt.gap_cutoff, opt.min_ins, opt.max_ins);
	mall->rid_beg = 0;
	mall->rid_end = n_vis;
	mall->log  = log;
	mall->cnsf = cnsf;
	mall->samf = samf;
	mall->msaf = msaf;
	mall->lmt  = opt.limit;
	mall->task    = 0;
	thread_end_init(mall);
	fprintf(log, "Indexing\n");
	fflush(log);
	thread_beg_iter(mall);
	mall->task = 1; // index
	thread_wake(mall);
	thread_end_iter(mall);
	thread_waitfor_all_idle(mall);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fprintf(log, "Local assembling\n");
	fflush(log);
	thread_beg_iter(mall);
	mall->task = 2; // alignand asm
	thread_wake(mall);
	thread_end_iter(mall);
	thread_waitfor_all_idle(mall);
	thread_beg_close(mall);
	free_lgraph(mall->g);
	thread_end_close(mall);
	sr_free_sdb(sdb);
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	fclose(log);
	fclose(cnsf);
	pclose(samf);
	pclose(msaf);
	return 0;
}


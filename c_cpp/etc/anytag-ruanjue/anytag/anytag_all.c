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

typedef struct {
	uint32_t rid;
	uint16_t offset:9, sign:1, n_mm:4, dir1:1, dir2:1;
} SimpleHit;

define_list(hitv, SimpleHit);

static inline void sr_hit2simp(SR_AlnHit *hit, SimpleHit *h, uint32_t ar_id){
	if(ar_id == hit->rid1){
		h->dir1 = hit->dir1; h->dir2 = hit->dir2; h->offset = hit->off; h->sign = 0; h->n_mm = hit->n_mm; h->rid = hit->rid2;
	} else {
		h->dir1 = !hit->dir2; h->dir2 = !hit->dir1; h->offset = hit->off; h->sign = 1; h->n_mm = hit->n_mm; h->rid = hit->rid1;
	}
}

static inline void write_sr_hit2simp(SR_AlnHit *hit, uint32_t ar_id, FILE *out){
	uint8_t byte;
	uint32_t rid;
	if(ar_id == hit->rid1){
		byte = (0x80U) | (hit->n_mm << 2) | ((hit->dir1 & 0x01) << 1) | (hit->dir2 & 0x01); fwrite(&byte, 1, 1, out);
		byte = hit->off; fwrite(&byte, 1, 1, out);
		rid = hit->rid2;
	} else {
		byte = (0x80U) | (hit->n_mm << 2) | ((!hit->dir2) << 1) | (!hit->dir1); fwrite(&byte, 1, 1, out);
		byte = hit->off | 0x80U; fwrite(&byte, 1, 1, out);
		rid = hit->rid1;
	}
	byte = (rid >>  0) & 0xFFU; fwrite(&byte, 1, 1, out);
	byte = (rid >>  8) & 0xFFU; fwrite(&byte, 1, 1, out);
	byte = (rid >> 16) & 0xFFU; fwrite(&byte, 1, 1, out);
	byte = (rid >> 24) & 0xFFU; fwrite(&byte, 1, 1, out);
}

static inline int read_simphit(SimpleHit *hit, FILE *in){
	uint8_t byte;
	if(fread(&byte, 1, 1, in) == 0) return 0;
	if(byte == 0) return 0;
	hit->dir1 = (byte >> 1) & 0x01;
	hit->dir2 = (byte >> 0) & 0x01;
	hit->n_mm = (byte >> 2) & 0x0F;
	hit->rid = 0;
	if(fread(&byte, 1, 1, in) == 0) return 0; hit->offset = byte & 0x7F; hit->sign = byte >> 7;
	if(fread(&byte, 1, 1, in) == 0) return 0; hit->rid |= ((uint32_t)byte) <<  0;
	if(fread(&byte, 1, 1, in) == 0) return 0; hit->rid |= ((uint32_t)byte) <<  8;
	if(fread(&byte, 1, 1, in) == 0) return 0; hit->rid |= ((uint32_t)byte) << 16;
	if(fread(&byte, 1, 1, in) == 0) return 0; hit->rid |= ((uint32_t)byte) << 24;
	return 1;
}

static inline uint32_t shuffle_hits(BitVec *flags, hitv *hits, uint32_t m){
	uint32_t i, j, n;
	double r;
	n = hits->size;
	if(m >= n){ ones_bitvec(flags); return n; }
	zeros_bitvec(flags);
	r = ((double)m) / n;
	i = 0; j = 0;
	while(j < m){
		if(get_bitvec(flags, i) == 0 && drand48() <= r){
			one_bitvec(flags, i); j ++;
		}
		if(i == n -1) i = 0;
		else i ++;
	}
	return j;
}

thread_beg_def(mall);
SR_SeqDB *sdb;
uint32_t rid_beg, rid_end;
ATOptions *opt;
LGraph *g;
uint32_t max_n_sr;
LibInserts *libs;
int task;
FILE *cnsf, *msaf, *alnf, *dotf, *log;
char *prefix;
BitVec *pair_flags;
thread_end_def(mall);

thread_beg_func(mall);
SR_SeqDB *sdb;
SR_AlnAux *aux1, *aux2;
LGraph *g;
LibInserts *libs;
hitv *hits1, *hits2;
SimpleHit H;
u32list *rid_used;
char seq1[1024], seq2[1024];
FILE *cnsf, *msaf, *alnf, *dotf, *log, *skef;
char *cns_buf, *msa_buf, *aln_buf, *dot_buf;
size_t cns_size, msa_size, aln_size, dot_size;
uint64_t rd_off1, rd_off2;
uint32_t i, j, rd_len1, rd_len2, t_idx, n_cpu, n_full, l, r, n, m;
int lib_id, min_ins, max_ins, off, lib_id0;
time_t t1, t2;
sdb = mall->sdb;
g = mall->g;
t_idx = mall->t_idx;
n_cpu = mall->n_cpu;
n_full = (mall->opt->n_full + mall->opt->n_cpu - 1) / mall->opt->n_cpu;
aux1 = sr_init_aux(1, MAX_RD_LEN, mall->opt->min_ol[0], mall->opt->max_mm[0], mall->opt->min_sm[0], mall->opt->allow_gap[0], 0, mall->opt->max_hits[0]);
aux2 = sr_init_aux(1, MAX_RD_LEN, mall->opt->min_ol[0], mall->opt->max_mm[0], mall->opt->min_sm[0], mall->opt->allow_gap[0], 0, mall->opt->max_hits[0]);
sr_fit_aux2sdb(aux1, sdb);
sr_fit_aux2sdb(aux2, sdb);
rid_used = init_u32list(1024);
hits1 = init_hitv(1024); hits2 = init_hitv(1024);
cns_buf = NULL; msa_buf = NULL; aln_buf = NULL; dot_buf = NULL;
cns_size = 0; msa_size = 0; aln_size = 0; dot_size = 0;
cnsf = open_memstream(&cns_buf, &cns_size);
msaf = open_memstream(&msa_buf, &msa_size);
alnf = open_memstream(&aln_buf, &aln_size);
dotf = open_memstream(&dot_buf, &dot_size);
log  = mall->log;
libs = mall->libs;
skef = NULL;
thread_beg_loop(mall);
if(mall->task == 1){
	sr_mask_low_complexity(sdb, n_cpu, t_idx);
} else if(mall->task == 2){
	for(i=t_idx;i<sdb->n_idx/2;i+=n_cpu){
		fprintf(log, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		sr_index_sdb(sdb, i);
		fprintf(log, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
	}
} else if(mall->task == 21){
	for(i=t_idx+sdb->n_idx/2;i<sdb->n_idx;i+=n_cpu){
		fprintf(log, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		sr_index_sdb(sdb, i);
		fprintf(log, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
	}
} else if(mall->task == 3){
	t1 = time(NULL);
	n = m = 0;
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL); fprintf(log, "[THREAD%03u] assembled %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(log);
		}
		if(get_bitvec(mall->pair_flags, i) == 1) continue;
		lib_id0 = lib_id =  get_read_lib(libs, 2 * i);
		min_ins = get_u32list(libs->min_ins, lib_id); max_ins = get_u32list(libs->max_ins, lib_id);
		g->min_ins = min_ins;
		g->max_ins = max_ins;
		sr_clear_aux(aux1); sr_align_sdb(sdb, 2 * i + 0, aux1);
		l = count_sr_hitv(aux1->hits);
		if(l > mall->max_n_sr) continue;
		sr_clear_aux(aux2); sr_align_sdb(sdb, 2 * i + 1, aux2);
		r = count_sr_hitv(aux2->hits);
		if(r > mall->max_n_sr) continue;
		clear_hitv(hits1);
		for(j=0;j<l;j++) sr_hit2simp(ref_sr_hitv(aux1->hits, j), next_ref_hitv(hits1), 2 * i);
		clear_hitv(hits2);
		for(j=0;j<r;j++) sr_hit2simp(ref_sr_hitv(aux2->hits, j), next_ref_hitv(hits2), 2 * i + 1);
		reset_lgraph(g);
		rd_off1 = sr_rdseq_offset(sdb, 2 * i); rd_len1 = sr_rdseq_length(sdb, 2 * i);
		rd_off2 = sr_rdseq_offset(sdb, 2 * i + 1); rd_len2 = sr_rdseq_length(sdb, 2 * i + 1);
		bits2seq(seq1, ref_u64list(sdb->rd_seqs, 0), rd_off1, rd_len1);
		bits2revseq(seq2, ref_u64list(sdb->rd_seqs, 0), rd_off2, rd_len2);
		push_lgraph(g, 2 * i + 0, seq1, rd_len1, 0, 0, 0, 0); push_lgraph(g, 2 * i + 1, seq2, rd_len2, 1, 1, 0, 0);
		if(m < n_full) fprintf(alnf, "A\t%u\t%u\t%u\t%s\t%s\n", i, l, r, seq1, seq2);
		for(j=0;j<l+r;j++){
			if(j < l){ H = get_hitv(hits1, j); }
			else { H = get_hitv(hits2, j - l); }
			lib_id = get_read_lib(libs, H.rid);
			if(lib_id < lib_id0) continue;
			min_ins = get_u32list(libs->min_ins, lib_id);
			max_ins = get_u32list(libs->max_ins, lib_id);
			rd_off1 = sr_rdseq_offset(sdb, H.rid); rd_len1 = sr_rdseq_length(sdb, H.rid);
			rd_off2 = sr_rdseq_offset(sdb, H.rid ^ 0x1U); rd_len2 = sr_rdseq_length(sdb, H.rid ^ 0x1U);
			off = (H.sign ^ H.dir1)? - ((int)H.offset) : (int)H.offset;
			if(j < l){
				bits2seq(seq1, sdb->rd_seqs->buffer, rd_off1, rd_len1);
				bits2revseq(seq2, sdb->rd_seqs->buffer, rd_off2, rd_len2);
				push_lgraph(g,        H.rid, seq1, rd_len1,  H.dir1, 0, off, off);
				push_lgraph(g, H.rid ^ 0x1U, seq2, rd_len2, !H.dir1, 0, min_ins + off - rd_len2, max_ins + off - rd_len2);
			} else {
				bits2revseq(seq1, sdb->rd_seqs->buffer, rd_off1, rd_len1);
				bits2seq(seq2, sdb->rd_seqs->buffer, rd_off2, rd_len2);
				push_lgraph(g,        H.rid, seq1, rd_len1, !H.dir1, 1, off, off);
				push_lgraph(g, H.rid ^ 0x1U, seq2, rd_len2,  H.dir1, 1, min_ins + off - rd_len2, max_ins + off - rd_len2);
			}
			if(m < n_full) fprintf(alnf, "%c\t%u\t%d\t%d\t%d\t%u\t%s\t%s\n", "LR"[j<l? 0 : 1], H.rid, off, min_ins, max_ins, H.n_mm, seq1, seq2);
		}
		ready_lgraph(g);
		if(layout_lgraph(g)){
			skeleton_lgraph(g);
			align_skeleton_lgraph(g, sdb);
			consensus_lgraph(g);
			for(j=0;j<g->sp_hits->size;j+=2){
				push_u32list(rid_used, ref_sr_hitv(g->sp_hits, j)->rid2 >> 1);
			}
			output_lgraph(g, cnsf, (m < n_full)? msaf : NULL);
			n ++;
		}
		if(m < n_full){ simplify_lgraph(g); print_dot_lgraph(g, dotf); }
		if((m & 0xFF) == 0){
			thread_beg_syn(mall);
			fflush(cnsf); fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf); fseek(cnsf, 0, SEEK_SET);
			fflush(msaf); fwrite(msa_buf, ftell(msaf), 1, mall->msaf); fseek(msaf, 0, SEEK_SET);
			fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET);
			fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET);
			for(j=0;j<rid_used->size;j++) one_bitvec(mall->pair_flags, get_u32list(rid_used, j));
			clear_u32list(rid_used);
			thread_end_syn(mall);
		}
	}
	t2 = time(NULL);
	fprintf(log, "[THREAD%03u] %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(log);
	if(1){
		thread_beg_syn(mall);
		fflush(cnsf); fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf); fseek(cnsf, 0, SEEK_SET);
		fflush(msaf); fwrite(msa_buf, ftell(msaf), 1, mall->msaf); fseek(msaf, 0, SEEK_SET);
		fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET);
		fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET);
		for(j=0;j<rid_used->size;j++) one_bitvec(mall->pair_flags, get_u32list(rid_used, j));
		clear_u32list(rid_used);
		thread_end_syn(mall);
	}
} else {
	microsleep(1);
}
thread_end_loop(mall);
sr_free_aux(aux1);
sr_free_aux(aux2);
free_hitv(hits1);
free_hitv(hits2);
free_u32list(rid_used);
fclose(cnsf); fclose(msaf);
if(alnf) fclose(alnf);
if(cns_buf) free(cns_buf);
if(msa_buf) free(msa_buf);
if(aln_buf) free(aln_buf);
thread_end_func(mall);

int debug_assembling(ATOptions opt, char *prefix){
	FileReader *fr;
	LGraph *g;
	FILE *cnsf, *msaf, *dotf;
	uint32_t id, n_lr, rd_len1, rd_len2, m, p;
	int n, off, max_ins, min_ins, p_dot;
	char *seq1, *seq2;
	clock_t t1, t2;
	char cmd[256];
	p_dot = 0;
	if((fr = fopen_filereader(opt.load_aln)) == NULL){
		fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", opt.load_aln, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	g = init_lgraph(&opt);
	sprintf(cmd, "%s.fasta", prefix);
	if((cnsf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "gzip -c -1 > %s.msa.gz", prefix);
	if((msaf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	if(p_dot){
		sprintf(cmd, "gzip -c -1 > %s.dot.gz", prefix);
		if((dotf = popen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			return 1;
		}
	}
	//print_options(&opt, stdout);
	m = p = 0;
	t1 = clock();
	while(1){
		if((n = fread_table(fr)) == -1) break;
		m ++;
		if((m % 100) == 0){
			t2 = clock();
			fprintf(stderr, "[%s] %u / %u clusters, %0.2f sec. \n", __FUNCTION__, p, m, 1.0 * (t2 - t1) / CLOCKS_PER_SEC); fflush(stderr);
		}
		reset_lgraph(g);
		id = atoll(get_col_str(fr, 1));
		n_lr = atoi(get_col_str(fr, 2)) + atoi(get_col_str(fr, 3));
		push_lgraph(g, 2 * id + 0, get_col_str(fr, 4), get_col_len(fr, 4), 0, 0, 0, 0);
		push_lgraph(g, 2 * id + 1, get_col_str(fr, 5), get_col_len(fr, 5), 1, 1, 0, 0);
		while(fread_table(fr) != -1){
			if(get_col_str(fr, 0)[0] == 'A'){ froll_back(fr); break; }
			id = atoll(get_col_str(fr, 1));
			off = atoi(get_col_str(fr, 2));
			min_ins = atoi(get_col_str(fr, 3));
			max_ins = atoi(get_col_str(fr, 4));
			seq1 = get_col_str(fr, 6);
			seq2 = get_col_str(fr, 7);
			rd_len1 = get_col_len(fr, 6);
			rd_len2 = get_col_len(fr, 7);
			push_lgraph(g,        id, seq1, rd_len1, 0, (get_col_str(fr, 0)[0] == 'R'), off, off);
			push_lgraph(g, id ^ 0x1U, seq2, rd_len2, 0, (get_col_str(fr, 0)[0] == 'R'), min_ins + off - rd_len2, max_ins + off - rd_len2);
		}
		ready_lgraph(g);
		//while(next_align_lgraph(g));
		if(layout_lgraph(g)){
			p ++;
			skeleton_lgraph(g);
			consensus_lgraph(g);
			output_lgraph(g, cnsf, msaf);
			fflush(cnsf); fflush(msaf);
		}
		if(p_dot){ simplify_lgraph(g); print_dot_lgraph(g, dotf); }
	}
	free_lgraph(g);
	fclose(cnsf);
	pclose(msaf);
	if(p_dot) pclose(dotf);
	fclose_filereader(fr);
	return 0;
}

int main_all(int argc, char **argv){
	ATOptions opt;
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	LibInserts *libs;
	FileReader *s1, *s2;
	BitVec *pair_flags;
	Sequence *seq1, *seq2;
	FILE *log, *cnsf, *msaf, *alnf, *dotf, *pair1, *pair2;
	char *prefix, cmd[512], *date, *outseq;
	uint64_t tot_len;
	uint32_t i, j, rid, rrid, n_vis, max_n_sr;
	int ii, c, vis;
	time_t tm;
	thread_preprocess(mall);
	if((c = parse_options(&opt, 0, argc, argv)) == -1) return usage_all();
	if(opt.genome_size == 0){
		fprintf(stdout, "Genome Size MUST be provided\n");
		return usage_all();
	}
	if(opt.load_aln){
		if(c + 1 < argc) return usage_all();
		return debug_assembling(opt, argv[c]);
	}
	if(c == argc || (argc - c) % 2 != 1) return usage_all();
	prefix = argv[c];
	sfs = init_sfv(4);
	for(ii=c+1;ii<argc;ii++){
		sf1 = next_ref_sfv(sfs);
		parse_seqfile(sf1, argv[ii], opt.var_size);
	}
	qsort(ref_sfv(sfs, 0), count_sfv(sfs), sizeof(SeqFile), cmp_seqfile);
	if(opt.min_lib == 0){
		if(sfs->size > 1){
			opt.min_lib = ref_sfv(sfs, sfs->size - 1)->ins_size;
		} else {
			opt.min_lib = ref_sfv(sfs, 0)->ins_size;
		}
	}
	for(i=0;2*i<count_sfv(sfs);i++){
		if(strcasecmp(ref_sfv(sfs, 2 * i)->prefix, ref_sfv(sfs, 2 * i + 1)->prefix)
			|| ref_sfv(sfs, 2 * i)->ins_size != ref_sfv(sfs, 2 * i + 1)->ins_size
			|| ref_sfv(sfs, 2 * i)->pair_idx >= ref_sfv(sfs, 2 * i + 1)->pair_idx){
			fprintf(stderr, " -- Wrong pair \"%s\" and \"%s\" in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	libs = guess_lib_inserts(sfs);
	sdb = sr_init_sdb(opt.kmer_size[0], opt.n_seed[0], opt.rd_len, opt.low_cpx);
	for(i=0;i<opt.n_seed[0];i++){
		if((opt.seed_skips[0] >> i) & 0x01) sr_mask_seed_sdb(sdb, i);
	}
	sprintf(cmd, "%s.info", prefix);
	if((log = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		abort();
	}
	sprintf(cmd, "%s.fasta", prefix);
	if((cnsf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.msa.gz", prefix);
	if((msaf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.aln.gz", prefix);
	if((alnf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.dot.gz", prefix);
	if((dotf = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
	}
	print_options(&opt, log);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date); fflush(log);
	for(i=0;i<libs->n_lib;i++){
		fprintf(log, "Lib%u: %u - %u bp\n", i, get_u32list(libs->min_ins, i), get_u32list(libs->max_ins, i));
	}
	fflush(log);
	n_vis = 0;
	tot_len = 0;
	rrid = rid = 0;
	for(i=0;2*i<count_sfv(sfs);i++){
		push_u32list(libs->lib_offs, rid * 2);
		seq1 = seq2 = NULL;
		sf1 = ref_sfv(sfs, 2 * i);
		sf2 = ref_sfv(sfs, 2 * i + 1);
		vis = (sf1->ins_size >= (int)opt.min_lib);
		if((s1 = fopen_filereader(sf1->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf1->name, __FUNCTION__, __FILE__, __LINE__); abort();
		}
		if((s2 = fopen_filereader(sf2->name)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", sf2->name, __FUNCTION__, __FILE__, __LINE__); abort();
		}
		rrid = rid;
		while(1){
			if(!((sf1->is_fq)? fread_fastq_adv(&seq1, s1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq1, s1, FASTA_FLAG_NO_NAME))) break;
			if(!((sf2->is_fq)? fread_fastq_adv(&seq2, s2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq2, s2, FASTA_FLAG_NO_NAME))) break;
			if(seq1->seq.size > (int)opt.rd_clip && seq2->seq.size > (int)opt.rd_clip){
				sr_push_sdb(sdb, seq1->seq.string + opt.rd_clip, seq1->seq.size - opt.rd_clip);
				sr_push_sdb(sdb, seq2->seq.string + opt.rd_clip, seq2->seq.size - opt.rd_clip);
				tot_len += seq1->seq.size - opt.rd_clip;
				tot_len += seq2->seq.size - opt.rd_clip;
			}
			rid ++;
		}
		if(seq1) free_sequence(seq1);
		if(seq2) free_sequence(seq2);
		fclose_filereader(s1);
		fclose_filereader(s2);
		push_u32list(libs->lib_cnts, (rid - rrid) * 2);
		fprintf(log, "%s\t%d\t%d\t%d\t%d\t%s\n", sf1->prefix, sf1->ins_size, rid - rrid, rrid, rid, vis? "yes" : "no"); fflush(log);
		if(vis) n_vis = rid;
	}
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	sr_ready_sdb(sdb);
	pair_flags = init_bitvec(rid);
	opt.inc_ar = 1; // always
	if(opt.inc_ar == 0){
		sr_set_n_ar(sdb, 2 * n_vis);
		max_n_sr = (2LLU * (rid - n_vis) * (2 * (sdb->max_rd_len - opt.min_ol[0] + 1)) / 2 / opt.genome_size) * 1.75;
	} else {
		max_n_sr = (2LLU * (rid) * (2 * (sdb->max_rd_len - opt.min_ol[0] + 1)) / 2 / opt.genome_size) * 1.75;
	}
	fprintf(log, "MAX_SR_NUM set to %d\n", max_n_sr); fflush(log);
	opt.limit = max_n_sr * 2;
	opt.max_hits[0] = max_n_sr * 1.2 + 1;
	thread_beg_init(mall, (int)opt.n_cpu);
	mall->opt = &opt;
	mall->max_n_sr = max_n_sr;
	mall->sdb = sdb;
	mall->pair_flags = pair_flags;
	mall->g = init_lgraph(&opt);
	mall->libs = libs;
	mall->rid_beg = 0;
	mall->rid_end = n_vis;
	mall->log  = log;
	mall->cnsf = cnsf;
	mall->msaf = msaf;
	mall->alnf = alnf;
	mall->dotf = dotf;
	mall->task    = 0;
	mall->prefix  = prefix;
	thread_end_init(mall);
	fprintf(log, "Filtering low complexity sequences\n"); fflush(log);
	thread_apply_all(mall, mall->task = 1);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fprintf(log, "Indexing\n"); fflush(log);
	thread_apply_all(mall, mall->task = 2);
	thread_apply_all(mall, mall->task = 21);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fprintf(log, "Local assembling\n"); fflush(log);
	thread_apply_all(mall, mall->task = 3);
	thread_beg_close(mall);
	free_lgraph(mall->g);
	thread_end_close(mall);
	rid = 0;
	rrid = 0;
	outseq = malloc(sdb->max_rd_len + 1);
	for(i=0;i<libs->n_lib;i++){
		sprintf(cmd, "%s.unused.ins%d.var%d.s1.fa", prefix, ref_sfv(sfs, i * 2 + 0)->ins_size, ref_sfv(sfs, i * 2 + 0)->var_size);
		if((pair1 = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		sprintf(cmd, "%s.unused.ins%d.var%d.s2.fa", prefix, ref_sfv(sfs, i * 2 + 1)->ins_size, ref_sfv(sfs, i * 2 + 1)->var_size);
		if((pair2 = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
		}
		rrid = get_u32list(libs->lib_cnts, i) / 2;
		for(j=0;j<rrid;j++){
			if(get_bitvec(pair_flags, j + rid)) continue;
			bits2seq(outseq, sdb->rd_seqs->buffer, sr_rdseq_offset(sdb, (j + rid) * 2 + 0), sr_rdseq_length(sdb, (j + rid) * 2 + 0));
			fprintf(pair1, ">%u/1\n%s\n", j, outseq);
			bits2seq(outseq, sdb->rd_seqs->buffer, sr_rdseq_offset(sdb, (j + rid) * 2 + 1), sr_rdseq_length(sdb, (j + rid) * 2 + 1));
			fprintf(pair2, ">%u/2\n%s\n", j, outseq);
		}
		rid += rrid;
		fclose(pair1);
		fclose(pair2);
	}
	free(outseq);
	sr_free_sdb(sdb);
	free_lib_inserts(libs);
	free_bitvec(pair_flags);
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date); fflush(log);
	fclose(log);
	fclose(cnsf);
	pclose(msaf);
	if(alnf) pclose(alnf);
	return 0;
}


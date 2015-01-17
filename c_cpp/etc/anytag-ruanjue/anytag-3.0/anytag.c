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
#include "ps_assembly.h"
#include "file_reader.h"
#include "thread.h"
#include "anytag_aux.h"
#include "counting_bloom_filter.h"

#define TASK_MASK	1
#define TASK_INDEX1	2
#define TASK_INDEX2	3
#define TASK_ALNASM	4
#define TASK_GAPCLOSE	5

int usage_all();

thread_beg_def(mall);
SR_SeqDB *sdb;
uint32_t rid_beg, rid_end;
ATOptions *opt;
Graph *g;
CBF *freqs;
uint32_t avg_n_sr;
LibInserts *libs;
int task;
FILE *cnsf, *msaf, *alnf, *dotf, *log;
char *prefix;
BitVec *pair_flags;
uint32_t once_id, once_mode;
String *scaff_seq, *scaff_head;
thread_end_def(mall);

thread_beg_func(mall);
SR_SeqDB *sdb;
SR_AlnAux *aux1, *aux2;
Graph *g;
LibInserts *libs;
SR_AlnHit H;
u32list *rid_used;
u32list *marks;
char seq1[1024], seq2[1024];
FILE *cnsf, *alnf, *msaf, *dotf, *log;
char *cns_buf, *msa_buf, *aln_buf, *dot_buf;
size_t cns_size, msa_size, aln_size, dot_size;
uint64_t rd_off1, rd_off2;
uint32_t i, j, k, b ,e, rd_len1, rd_len2, t_idx, n_cpu, n_full, l, r, n, m;
int lib_id, min_ins, max_ins, off, lib_id0;
time_t t1, t2;
sdb = mall->sdb;
g = mall->g;
t_idx = mall->t_idx;
n_cpu = mall->n_cpu;
n_full = (mall->opt->n_full + mall->opt->n_cpu - 1) / mall->opt->n_cpu;

marks = init_u32list(32);

aux1 = sr_init_aux();
sr_set_aux_strand(aux1, 1);
sr_set_aux_overlap(aux1, mall->opt->min_ol[0], sdb->rd_len);
sr_set_aux_min_similarity(aux1, mall->opt->min_sm[0]);
sr_set_aux_gap(aux1, 0);
sr_set_aux_max_hits(aux1, mall->avg_n_sr * ((mall->opt->_asm_flags[0] == 1)? 2 : 8));
sr_set_aux_cigar(aux1, 0);
sr_fit_aux2sdb(aux1, sdb);

aux2 = sr_init_aux();
sr_set_aux_strand(aux2, 1);
sr_set_aux_overlap(aux2, mall->opt->min_ol[0], sdb->rd_len);
sr_set_aux_min_similarity(aux2, mall->opt->min_sm[0]);
sr_set_aux_gap(aux2, 0);
sr_set_aux_max_hits(aux2, mall->avg_n_sr * ((mall->opt->_asm_flags[0] == 1)? 2 : 8));
sr_set_aux_cigar(aux2, 0);
sr_fit_aux2sdb(aux2, sdb);

rid_used = init_u32list(1024);
cns_buf = NULL; msa_buf = NULL; aln_buf = NULL; dot_buf = NULL;
cns_size = 0; msa_size = 0; aln_size = 0; dot_size = 0;
cnsf = open_memstream(&cns_buf, &cns_size);
msaf = open_memstream(&msa_buf, &msa_size);
alnf = open_memstream(&aln_buf, &aln_size);
dotf = open_memstream(&dot_buf, &dot_size);
log  = mall->log;
libs = mall->libs;

thread_beg_loop(mall);
if(mall->task == TASK_MASK){
	sr_mask_low_complexity(sdb, n_cpu, t_idx);
} else if(mall->task == TASK_INDEX1){
	for(i=t_idx;i<sdb->n_idx/2;i+=n_cpu){
		fprintf(log, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		fprintf(stderr, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(stderr);
		sr_index_sdb(sdb, i);
		fprintf(log, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		fprintf(stderr, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(stderr);
	}
} else if(mall->task == TASK_INDEX2){
	for(i=t_idx+sdb->n_idx/2;i<sdb->n_idx;i+=n_cpu){
		fprintf(log, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		fprintf(stderr, "[THREAD%03u] indexing %u/%u. \n", t_idx, i, sdb->n_idx); fflush(stderr);
		sr_index_sdb(sdb, i);
		fprintf(log, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(log);
		fprintf(stderr, "[THREAD%03u] indexed %u/%u. \n", t_idx, i, sdb->n_idx); fflush(stderr);
	}
} else if(mall->task == TASK_ALNASM){
	t1 = time(NULL);
	n = m = 0;
	if(mall->once_mode){
		i = mall->once_id;
		goto ONCE_ENTRY;
	}
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] assembled %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(log);
			fprintf(stderr, "[THREAD%03u] assembled %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(stderr);
		}
		if(get_bitvec(mall->pair_flags, i) == 1) continue;
ONCE_ENTRY:
		lib_id0 = lib_id =  get_read_lib(libs, 2 * i);
		min_ins = get_u32list(libs->min_ins, lib_id); max_ins = get_u32list(libs->max_ins, lib_id);
		set_insert_graph(g, min_ins, max_ins);
		sr_clear_aux(aux1);
		if(mall->opt->inc_ar || lib_id + 1U == libs->n_lib){
			sr_set_aux_hit_id_region(aux1, get_u32list(libs->lib_offs, lib_id), 0xFFFFFFFFU);
		} else {
			sr_set_aux_hit_id_region(aux1, get_u32list(libs->lib_offs, lib_id + 1), 0xFFFFFFFFU);
		}
		sr_aux_load(aux1, sdb, 2 * i + 0, sdb->rdseqs->bits, sr_rdseq_offset(sdb, 2 * i + 0), sdb->rd_len);
		sr_align_sdb(sdb, aux1);
		l = count_sr_hitv(aux1->hits);
		if(mall->opt->_asm_flags[0] == 1 && l >= 1.5 * mall->avg_n_sr){
			if(mall->once_mode) break;
			else continue;
		}
		sr_clear_aux(aux2);
		if(mall->opt->inc_ar || lib_id + 1U == libs->n_lib){
			sr_set_aux_hit_id_region(aux2, get_u32list(libs->lib_offs, lib_id), 0xFFFFFFFFU);
		} else {
			sr_set_aux_hit_id_region(aux2, get_u32list(libs->lib_offs, lib_id + 1), 0xFFFFFFFFU);
		}
		sr_aux_load(aux2, sdb, 2 * i + 1, sdb->rdseqs->bits, sr_rdseq_offset(sdb, 2 * i + 1), sdb->rd_len);
		sr_align_sdb(sdb, aux2);
		r = count_sr_hitv(aux2->hits);
		if(mall->opt->_asm_flags[0] == 1 && r >= 1.5 * mall->avg_n_sr){
			if(mall->once_mode) break;
			else continue;
		}
		if(mall->opt->_asm_flags[0] == 2 && (r >= 1.5 * mall->avg_n_sr && l >= 1.5 * mall->avg_n_sr)){
			if(mall->once_mode) break;
			else continue;
		}
		if(mall->opt->kmer_size[1] == 0){
			if(mall->once_mode) break;
			else continue; // ksize=0 to test the speed of clustering
		}
		reset_graph(g);
		rd_off1 = sr_rdseq_offset(sdb, 2 * i + 0); rd_len1 = sr_rdseq_length(sdb, 2 * i + 0);
		rd_off2 = sr_rdseq_offset(sdb, 2 * i + 1); rd_len2 = sr_rdseq_length(sdb, 2 * i + 1);
		bits2seq(seq1, sdb->rdseqs->bits, rd_off1, rd_len1);
		bits2revseq(seq2, sdb->rdseqs->bits, rd_off2, rd_len2);
		push_graph(g, 2 * i + 0, seq1, rd_len1, 0, 0); push_graph(g, 2 * i + 1, seq2, rd_len2, 1, 0);
		if(m < n_full) fprintf(alnf, "A\t%u\t%u\t%u\t%s\t%s\n", i, l, r, seq1, seq2);
		for(j=0;j<l+r;j++){
			if(j < l){ H = get_sr_hitv(aux1->hits, j); }
			else { H = get_sr_hitv(aux2->hits, j - l); }
			lib_id = get_read_lib(libs, H.rid);
			if(lib_id < lib_id0) continue;
			rd_off1 = sr_rdseq_offset(sdb, H.rid); rd_len1 = sr_rdseq_length(sdb, H.rid);
			rd_off2 = sr_rdseq_offset(sdb, H.rid ^ 0x1U); rd_len2 = sr_rdseq_length(sdb, H.rid ^ 0x1U);
			off = (H.dir1)? - ((int)H.off) : (int)H.off;
			if(j < l){
				bits2seq(seq1, sdb->rdseqs->bits, rd_off1, rd_len1);
				bits2revseq(seq2, sdb->rdseqs->bits, rd_off2, rd_len2);
				push_graph(g,        H.rid, seq1, rd_len1,  H.dir1, 1 | ((off < 0)? 4 : 0));
				push_graph(g, H.rid ^ 0x1U, seq2, rd_len2, !H.dir1, 0);
			} else {
				bits2revseq(seq1, sdb->rdseqs->bits, rd_off1, rd_len1);
				bits2seq(seq2, sdb->rdseqs->bits, rd_off2, rd_len2);
				push_graph(g,        H.rid, seq1, rd_len1, !H.dir1, 2 | ((off < 0)? 4 : 0));
				push_graph(g, H.rid ^ 0x1U, seq2, rd_len2,  H.dir1, 0);
			}
			if(m < n_full) fprintf(alnf, "%c\t%u\t%d\t%u\t%s\t%s\n", "LR"[j<l? 0 : 1], H.rid, off, H.n_mm, seq1, seq2);
		}
		ready_graph(g);
		index_graph(g);
		align_graph(g);
		simplify_graph(g);
		shave_graph(g);
		allpaths_graph(g);
		validate_paths_graph(g);
		consensus_graph(g);
		realign_graph(g, sdb, libs->lib_cnts, libs->lib_ins, mall->opt->re_cns);
		//if(mall->freqs) kmer_qc_graph(g, mall->freqs, mall->opt->qc_ksize, 2);
		output_contigs_graph(g, cnsf, (m < n_full)? msaf : NULL);
		if(g->cnss->size) n++;
		for(j=0;j<g->cnss->size;j++){
			for(k=0;k<g->cnss->buffer[j].path->size;k++){
				push_u32list(rid_used, ref_nodev(g->nodes, g->cnss->buffer[j].path->buffer[k].nid)->rid >> 1);
			}
		}
		if(m < n_full){
			if(g->cnss->size == 0){
				fprintf(dotf, "//Failed\n");
			} else if(g->cnss->size == 1){
				fprintf(dotf, "//Success\n");
			} else {
				fprintf(dotf, "//Multiple\n");
			}
			print_dot_graph(g, dotf);
		}
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
		if(mall->once_mode) goto ONCE_EXIT;
	}
	t2 = time(NULL);
	fprintf(log, "[THREAD%03u] finished %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(log);
	fprintf(stderr, "[THREAD%03u] finished %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1)); fflush(stderr);
ONCE_EXIT:
	if(1){
		thread_beg_syn(mall);
		fflush(cnsf); fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf); fflush(mall->cnsf); fseek(cnsf, 0, SEEK_SET);
		fflush(msaf); fwrite(msa_buf, ftell(msaf), 1, mall->msaf); fflush(mall->msaf); fseek(msaf, 0, SEEK_SET);
		fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fflush(mall->alnf); fseek(alnf, 0, SEEK_SET);
		fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fflush(mall->dotf); fseek(dotf, 0, SEEK_SET);
		for(j=0;j<rid_used->size;j++) one_bitvec(mall->pair_flags, get_u32list(rid_used, j));
		clear_u32list(rid_used);
		thread_end_syn(mall);
	}
} else if(mall->task == TASK_GAPCLOSE){
	clear_u32list(marks);
	j = 0;
	for(i=0;i<(unsigned)mall->scaff_seq->size;i++){
		if(mall->scaff_seq->string[i] == 'N'){
			if(j == 0){ j = 1; push_u32list(marks, i); }
		} else {
			if(j == 1){ j = 0; push_u32list(marks, i); }
		}
	}
	for(i=0;i+1<marks->size;i+=2){
		b = get_u32list(marks, i);
		e = get_u32list(marks, i + 1);
	}
} else {
	microsleep(1);
}
thread_end_loop(mall);
sr_free_aux(aux1);
sr_free_aux(aux2);
free_u32list(rid_used);
free_u32list(marks);
fclose(cnsf); fclose(msaf); fclose(alnf); fclose(dotf);
if(cns_buf) free(cns_buf);
if(msa_buf) free(msa_buf);
if(aln_buf) free(aln_buf);
if(dot_buf) free(dot_buf);
thread_end_func(mall);

void build_kmer_freq_table(SR_SeqDB *sdb, CBF *bf, uint32_t ksize){
	uint64_t i, k, r, kmask;
	uint64_t off;
	uint32_t j;
	kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32 - ksize) << 1);
	for(i=0;i<sdb->n_rd;i++){
		off = sr_rdseq_offset(sdb, i);
		k = 0;
		for(j=0;j<sdb->rd_len;j++){
			k = ((k << 2) | bits2bit(sdb->rdseqs->bits, off + j)) & kmask;
			if(j + 1 < ksize) continue;
			r = dna_rev_seq(k, ksize);
			if(k < r) r = k;
			put_cbf(bf, &r, sizeof(uint64_t));
		}
	}
}

int main_all(int argc, char **argv){
	ATOptions opt;
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	LibInserts *libs;
	FileReader *s1, *s2, *s3;
	CBF *freqs;
	BitVec *pair_flags;
	Sequence *seq1, *seq2;
	FILE *log, *cnsf, *msaf, *alnf, *dotf, *dmpf, *pair1, *pair2;
	char *prefix, cmd[512], *date, *outseq;
	uint64_t tot_len;
	uint32_t i, j, rid, rrid, n_vis, avg_n_sr, rdlen;
	float gcov;
	int ii, c, vis;
	time_t tm;
	thread_preprocess(mall);
	if((c = parse_options(&opt, 0, argc, argv)) == -1) return usage_all();
	if(opt.genome_size == 0){
		fprintf(stdout, " ** Genome Size MUST be provided ** \n");
		return usage_all();
	}
	if(c == argc || (argc - c) % 2 != 1){
		fprintf(stdout, " ** Please check output_prefix or input files ** \n");
		return usage_all();
	}
	prefix = argv[c];
	if(opt.debug > 1){
		prefix = alloca(strlen(argv[c]) + 6);
		sprintf(prefix, "%s.dump", argv[c]);
	} else prefix = argv[c];
	{
		sprintf(cmd, "%s.log", prefix);
		if((log = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			abort();
		}
		sprintf(cmd, "%s.fasta", prefix);
		if((cnsf = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
		}
		sprintf(cmd, "%s.aln", prefix);
		if((alnf = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
		}
		sprintf(cmd, "%s.msa", prefix);
		if((msaf = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
		}
		sprintf(cmd, "%s.dot", prefix);
		if((dotf = fopen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
		}
		fprintf(log, "Version: "VERSION" "UID"\n");
		fprintf(stderr, "Version: "VERSION" "UID"\n");
		print_options(&opt, log);
		print_options(&opt, stderr);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date); fflush(log);
		fprintf(stderr, "%s\n", date); fflush(stderr);
	}
	if(opt.debug > 1){
		fprintf(log, "Loading\n"); fflush(log);
		fprintf(stderr, "Loading\n"); fflush(stderr);
		sprintf(cmd, "%s", prefix);
		if((dmpf = fopen(cmd, "r")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
		}
		sdb = sr_load_sdb(dmpf);
		fread(&rid, sizeof(uint32_t), 1, dmpf);
		fread(&n_vis, sizeof(uint32_t), 1, dmpf);
		fread(&avg_n_sr, sizeof(uint32_t), 1, dmpf);
		libs = malloc(sizeof(LibInserts));
		fread(&libs->n_lib, sizeof(uint32_t), 1, dmpf);
		libs->lib_offs = load_u32list(dmpf);
		libs->lib_cnts = load_u32list(dmpf);
		libs->min_ins  = load_u32list(dmpf);
		libs->max_ins  = load_u32list(dmpf);
		libs->lib_ins  = load_u32list(dmpf);
		libs->lib_vars = load_u32list(dmpf);
		fclose(dmpf);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date);
		fprintf(stderr, "%s\n", date);
	} else {
		sfs = init_sfv(4);
		for(ii=c+1;ii<argc;ii++){
			sf1 = next_ref_sfv(sfs);
			parse_seqfile(sf1, argv[ii], 0.3f);
		}
		qsort(ref_sfv(sfs, 0), count_sfv(sfs), sizeof(SeqFile), cmp_seqfile);
		if(opt.min_ar == 0){
			if(sfs->size > 1){
				opt.min_ar = ref_sfv(sfs, sfs->size - 3)->ins_size;
			} else {
				opt.min_ar = ref_sfv(sfs, 0)->ins_size;
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
		for(i=0;i<libs->n_lib;i++){
			fprintf(log, "Lib%u: %u - %u bp\n", i, get_u32list(libs->min_ins, i), get_u32list(libs->max_ins, i));
			fprintf(stderr, "Lib%u: %u - %u bp\n", i, get_u32list(libs->min_ins, i), get_u32list(libs->max_ins, i));
		}
		sdb = sr_init_sdb(opt.kmer_size[0], opt.n_seed, 0, opt.low_cpx);
		for(i=0;i<opt.n_seed;i++){
			if((opt.seed_skips >> i) & 0x01) sr_mask_seed_sdb(sdb, i);
		}
		n_vis = 0;
		tot_len = 0;
		rrid = rid = 0;
		if(opt.rd_trim) rdlen = opt.rd_trim;
		else rdlen = 0;
		for(i=0;2*i<count_sfv(sfs);i++){
			push_u32list(libs->lib_offs, rid * 2);
			seq1 = seq2 = NULL;
			sf1 = ref_sfv(sfs, 2 * i);
			sf2 = ref_sfv(sfs, 2 * i + 1);
			vis = (sf1->ins_size >= (int)opt.min_ar);
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
				if(rdlen == 0){
					rdlen = (seq1->seq.size < seq2->seq.size)? seq1->seq.size : seq2->seq.size;
				}
				if(seq1->seq.size >= (int)(rdlen) && seq2->seq.size >= (int)(rdlen)){
					sr_push_sdb(sdb, seq1->seq.string, rdlen);
					sr_push_sdb(sdb, seq2->seq.string, rdlen);
					tot_len += rdlen + rdlen;
				} else {
					fprintf(stderr, " -- Reads SHOULD be the same length, or use -r to trim longer reads in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
				}
				rid ++;
			}
			if(seq1) free_sequence(seq1);
			if(seq2) free_sequence(seq2);
			fclose_filereader(s1);
			fclose_filereader(s2);
			push_u32list(libs->lib_cnts, (rid - rrid) * 2);
			fprintf(log, "%s\t%d\t%d\t%d\t%d\t%s\n", sf1->prefix, sf1->ins_size, rid - rrid, rrid, rid, vis? "yes" : "no"); fflush(log);
			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%s\n", sf1->prefix, sf1->ins_size, rid - rrid, rrid, rid, vis? "yes" : "no"); fflush(stderr);
			if(vis) n_vis = rid;
		}
		for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
		free_sfv(sfs);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date); fflush(log);
		fprintf(stderr, "%s\n", date); fflush(stderr);
		sr_ready_sdb(sdb);
	}

	gcov = (((uint64_t)sdb->n_rd) * sdb->rd_len) / opt.genome_size / (1024.0 * 1024.0);
	avg_n_sr = ((sdb->rd_len - opt.min_ol[0] + 1) * (gcov / sdb->rd_len)) * 8;
	fprintf(log, "Average coverage set to %0.1f\n", gcov); fflush(log);
	fprintf(stderr, "Average coverage  set to %0.1f\n", gcov); fflush(stderr);
	fprintf(log, "AVG_SR_NUM set to %d\n", avg_n_sr); fflush(log);
	fprintf(stderr, "AVG_SR_NUM set to %d\n", avg_n_sr); fflush(stderr);
	pair_flags = init_bitvec(rid);
	if(opt.qc_ksize){
		fprintf(log, "Begin to build kmer freq table ksize = %u\n", opt.qc_ksize); fflush(log);
		fprintf(stderr, "Begin to build kmer freq table ksize = %u\n", opt.qc_ksize); fflush(stderr);
		freqs = init_cbf(((uint64_t)sdb->n_rd) * (sdb->rd_len + 1 - opt.qc_ksize), 2, 3);
		build_kmer_freq_table(sdb, freqs, opt.qc_ksize);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date);
		fprintf(stderr, "%s\n", date);
	} else freqs = NULL;
	thread_beg_init(mall, (int)opt.n_cpu);
	mall->opt = &opt;
	mall->avg_n_sr = avg_n_sr;
	mall->sdb = sdb;
	mall->freqs = freqs;
	mall->pair_flags = pair_flags;
	mall->g = init_graph(opt.kmer_size[1], opt.min_ol[1], opt.min_sm[1]);
	mall->g->out_mode = (opt._asm_flags[0] == 2);
	mall->libs = libs;
	mall->rid_beg = (opt.debug == 2)? opt.debug_beg : 0;
	mall->rid_end = n_vis;
	mall->log  = log;
	mall->cnsf = cnsf;
	mall->alnf = alnf;
	mall->msaf = msaf;
	mall->dotf = dotf;
	mall->task    = 0;
	mall->once_mode = 0;
	mall->once_id   = 0;
	mall->prefix  = prefix;
	if(opt.scaff){ mall->scaff_seq = init_string(1024 * 100); mall->scaff_head = init_string(1024); }
	else mall->scaff_head = mall->scaff_seq = NULL;
	thread_end_init(mall);
	if(opt.debug > 1){
	} else {
		fprintf(log, "Filtering low complexity sequences\n"); fflush(log);
		fprintf(stderr, "Filtering low complexity sequences\n"); fflush(stderr);
		thread_apply_all(mall, mall->task = TASK_MASK);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date);
		fprintf(stderr, "%s\n", date);
		fprintf(log, "Indexing\n"); fflush(log);
		fprintf(stderr, "Indexing\n"); fflush(stderr);
		thread_apply_all(mall, mall->task = TASK_INDEX1);
		thread_apply_all(mall, mall->task = TASK_INDEX2);
		tm = time(NULL); date = asctime(localtime(&tm));
		fprintf(log, "%s\n", date);
		fprintf(stderr, "%s\n", date);
		if(opt.debug == 1){
			fprintf(log, "Dumping\n"); fflush(log);
			fprintf(stderr, "Dumping\n"); fflush(stderr);
			sprintf(cmd, "%s.dump", prefix);
			if((dmpf = fopen(cmd, "w")) == NULL){
				fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); return 1;
			}
			sr_dump_sdb(sdb, dmpf);
			fwrite(&rid, sizeof(uint32_t), 1, dmpf);
			fwrite(&n_vis, sizeof(uint32_t), 1, dmpf);
			fwrite(&avg_n_sr, sizeof(uint32_t), 1, dmpf);
			fwrite(&libs->n_lib, sizeof(uint32_t), 1, dmpf);
			dump_u32list(libs->lib_offs, dmpf);
			dump_u32list(libs->lib_cnts, dmpf);
			dump_u32list(libs->min_ins, dmpf);
			dump_u32list(libs->max_ins, dmpf);
			dump_u32list(libs->lib_ins, dmpf);
			dump_u32list(libs->lib_vars, dmpf);
			fclose(dmpf);
			tm = time(NULL); date = asctime(localtime(&tm));
			fprintf(log, "%s\n", date);
			fprintf(stderr, "%s\n", date);
		}
	}
	fprintf(log, "Local assembling\n"); fflush(log);
	fprintf(stderr, "Local assembling\n"); fflush(stderr);
	if(opt.debug == 3){
		s3 = stdin_filereader();
		while(fread_table(s3) != -1){
			thread_waitfor_one_idle(mall);
			mall->task = TASK_ALNASM;
			mall->once_mode = 1;
			mall->once_id = atoll(get_col_str(s3, 0));
			thread_wake(mall);
		}
		fclose_filereader(s3);
		thread_waitfor_all_idle(mall);
	} else if(opt.scaff){
		if((s3 = fopen_filereader(opt.scaff)) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", opt.scaff, __FUNCTION__, __FILE__, __LINE__); abort();
		}
		seq1 = NULL;
		while(1){
			if(fread_fasta(&seq1, s3)) break;
			thread_waitfor_one_idle(mall);
			mall->task = TASK_GAPCLOSE;
			clear_string(mall->scaff_head);
			clear_string(mall->scaff_seq);
			append_string(mall->scaff_head, seq1->name.string, seq1->name.size);
			if(seq1->comment.size){
				add_char_string(mall->scaff_head, ' ');
				append_string(mall->scaff_head, seq1->comment.string, seq1->comment.size);
			}
			append_string(mall->scaff_seq, seq1->seq.string, seq1->seq.size);
			thread_wake(mall);
		}
		fclose_filereader(s3);
		thread_waitfor_all_idle(mall);
	} else {
		thread_apply_all(mall, mall->task = TASK_ALNASM);
	}
	thread_beg_close(mall);
	free_graph(mall->g);
	if(mall->scaff_head) free_string(mall->scaff_head);
	if(mall->scaff_seq) free_string(mall->scaff_seq);
	thread_end_close(mall);
	if(freqs) free_cbf(freqs);
	rid = 0;
	rrid = 0;
	if(opt.debug != 3 && opt.scaff == NULL){
		outseq = malloc(sdb->rd_len + 1);
		for(i=0;i<libs->n_lib;i++){
			sprintf(cmd, "%s.unused.ins%d.var%d.s1.fa", prefix, get_u32list(libs->lib_ins, i), get_u32list(libs->lib_vars, i));
			if((pair1 = fopen(cmd, "w")) == NULL){
				fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			sprintf(cmd, "%s.unused.ins%d.var%d.s2.fa", prefix, get_u32list(libs->lib_ins, i), get_u32list(libs->lib_vars, i));
			if((pair2 = fopen(cmd, "w")) == NULL){
				fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__); fflush(stderr);
			}
			rrid = get_u32list(libs->lib_cnts, i) / 2;
			for(j=0;j<rrid;j++){
				if(get_bitvec(pair_flags, j + rid)) continue;
				bits2seq(outseq, sdb->rdseqs->bits, sr_rdseq_offset(sdb, (j + rid) * 2 + 0), sr_rdseq_length(sdb, (j + rid) * 2 + 0));
				fprintf(pair1, ">%u/1\n%s\n", j, outseq);
				bits2seq(outseq, sdb->rdseqs->bits, sr_rdseq_offset(sdb, (j + rid) * 2 + 1), sr_rdseq_length(sdb, (j + rid) * 2 + 1));
				fprintf(pair2, ">%u/2\n%s\n", j, outseq);
			}
			rid += rrid;
			fclose(pair1);
			fclose(pair2);
		}
		free(outseq);
	}
	sr_free_sdb(sdb);
	free_lib_inserts(libs);
	free_bitvec(pair_flags);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date); fflush(log);
	fprintf(stderr, "%s\n", date); fflush(stderr);
	fclose(log);
	fclose(cnsf);
	fclose(msaf);
	fclose(alnf);
	fclose(dotf);
	return 0;
}


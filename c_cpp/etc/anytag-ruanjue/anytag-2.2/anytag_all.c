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
	uint16_t offset:10, n_mm:4, dir1:1, dir2:1;
} SimpleHit;

define_list(hitv, SimpleHit);

thread_beg_def(mall);
SR_SeqDB *sdb;
uint32_t rid_beg, rid_end;
uint32_t lmt;
LGraph *g;
LibInserts *libs;
int task;
FILE *cnsf, *msaf, *alnf, *dotf, *log;
char *prefix;
uint16_t *n_hits;
thread_end_def(mall);

thread_beg_func(mall);
SR_SeqDB *sdb;
SR_AlnAux *aux1, *aux2;
SR_AlnHit *hit;
LGraph *g;
Supp SP;
LibInserts *libs;
hitv *hits1, *hits2;
SimpleHit H;
char seq1[1024], seq2[1024], chs[512];
FILE *cnsf, *msaf, *alnf, *dotf, *log, *skef;
FILE **aln_parts;
FILE **cns_parts;
char *cns_buf, *msa_buf, *aln_buf, *dot_buf;
size_t cns_size, msa_size, aln_size, dot_size;
uint64_t rd_off1, rd_off2;
uint32_t i, j, rd_len1, rd_len2, t_idx, n_cpu, l, r, n, m;
int lib_id, min_ins, max_ins, off;
uint16_t *n_hits;
uint8_t byte;
time_t t1, t2;
sdb = mall->sdb;
g = mall->g;
t_idx = mall->t_idx;
n_cpu = mall->n_cpu;
aux1 = sr_init_aux();
aux2 = sr_init_aux();
hits1 = init_hitv(1024);
hits2 = init_hitv(1024);
cns_buf = NULL;
msa_buf = NULL;
aln_buf = NULL;
dot_buf = NULL;
cns_size = 0;
msa_size = 0;
aln_size = 0;
dot_size = 0;
cnsf = open_memstream(&cns_buf, &cns_size);
msaf = open_memstream(&msa_buf, &msa_size);
if(mall->alnf) alnf = open_memstream(&aln_buf, &aln_size);
else alnf = NULL;
if(mall->dotf) dotf = open_memstream(&dot_buf, &dot_size);
else dotf = NULL;
log  = mall->log;
libs = mall->libs;
n_hits = mall->n_hits;
skef = NULL;
memset(&SP, 0, sizeof(Supp));
if(sdb->n_part < 2){
	aln_parts = NULL;
	cns_parts = NULL;
} else {
	aln_parts = malloc(sizeof(FILE*) * (sdb->n_part));
	cns_parts = malloc(sizeof(FILE*) * (sdb->n_part));
}
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
	sprintf(chs, "%s.aln.p%03u.t%03u", mall->prefix, sdb->part_idx, t_idx);
	if((aln_parts[sdb->part_idx] = fopen(chs, "w+")) == NULL){
		fprintf(stdout, " -- Cannot create %s in %s -- %s:%d --\n", chs, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
	}
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 100000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] aligned %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
			fflush(log);
		}
		clear_sr_hitv(aux1->hits);
		if(n_hits[2 * i + 0] < mall->lmt / 2) sr_align_sdb(sdb, 2 * i + 0, aux1);
		clear_sr_hitv(aux2->hits);
		if(n_hits[2 * i + 1] < mall->lmt / 2) sr_align_sdb(sdb, 2 * i + 1, aux2);
		l = count_sr_hitv(aux1->hits);
		if(l > mall->lmt / sdb->n_part) l = mall->lmt / sdb->n_part;
		r = count_sr_hitv(aux2->hits);
		if(r > mall->lmt / sdb->n_part) r = mall->lmt / sdb->n_part;
		if(l + r) n ++;
		for(j=0;j<l;j++){
			hit = ref_sr_hitv(aux1->hits, j);
			byte = (0x80U) | (hit->n_mm << 2) | ((hit->dir1 & 0x01) << 1) | (hit->dir2 & 0x01); fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = hit->off; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >>  0) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >>  8) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >> 16) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >> 24) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
		}
		n_hits[2 * i + 0] += l;
		byte = 0; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
		for(j=0;j<r;j++){
			hit = ref_sr_hitv(aux2->hits, j);
			byte = (0x80U) | (hit->n_mm << 2) | ((hit->dir1 & 0x01) << 1) | (hit->dir2 & 0x01); fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = hit->off; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >>  0) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >>  8) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >> 16) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
			byte = (hit->rid2 >> 24) & 0xFFU; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
		}
		n_hits[2 * i + 1] += r;
		byte = 0; fwrite(&byte, 1, 1, aln_parts[sdb->part_idx]);
	}
} else if(mall->task == 3){
	t1 = time(NULL);
	n = m = 0;
	if(sdb->n_part > 1){
		sprintf(chs, "%s.cns.p%03u.t%03u", mall->prefix, sdb->part_idx, t_idx);
		if((cns_parts[sdb->part_idx] = fopen(chs, "w+")) == NULL){
			fprintf(stdout, " -- Cannot create %s in %s -- %s:%d --\n", chs, __FUNCTION__, __FILE__, __LINE__); fflush(stdout); abort();
		}
		sprintf(chs, "%s.cns.ske.t%03u", mall->prefix, t_idx);
		if((skef = fopen(chs, "w+")) == NULL){
			fprintf(stdout, " -- Cannot create %s in %s -- %s:%d --\n", chs, __FUNCTION__, __FILE__, __LINE__); fflush(stdout); abort();
		}
	}
	for(i=0;i+1<sdb->n_part;i++){ fseek(aln_parts[i], 0, SEEK_SET); }
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] assembled %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
			fflush(log);
		}
		clear_sr_hitv(aux1->hits);
		if(n_hits[2 * i + 0] < mall->lmt / 2) sr_align_sdb(sdb, 2 * i + 0, aux1);
		l = count_sr_hitv(aux1->hits);
		if(l > mall->lmt / sdb->n_part) l = mall->lmt / sdb->n_part;
		clear_sr_hitv(aux2->hits);
		if(n_hits[2 * i + 1] < mall->lmt / 2) sr_align_sdb(sdb, 2 * i + 1, aux2);
		r = count_sr_hitv(aux2->hits);
		if(r > mall->lmt / sdb->n_part) r = mall->lmt / sdb->n_part;
		clear_hitv(hits1);
		for(j=0;j+1<sdb->n_part;j++){
			while(1){
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break;
				if(byte == 0) break;
				H.dir1 = (byte >> 1) & 0x01;
				H.dir2 = (byte >> 0) & 0x01;
				H.n_mm = (byte >> 2) & 0x0F;
				H.rid = 0;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.offset = byte;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) <<  0;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) <<  8;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) << 16;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) << 24;
				push_hitv(hits1, H);
			}
		}
		for(j=0;j<l;j++){
			hit = ref_sr_hitv(aux1->hits, j);
			H.dir1   = hit->dir1;
			H.dir2   = hit->dir2;
			H.offset = hit->off;
			H.rid    = hit->rid2;
			push_hitv(hits1, H);
		}
		if(count_hitv(hits1) > mall->lmt / 2) set_hitv_size(hits1, mall->lmt / 2);
		clear_hitv(hits2);
		for(j=0;j+1<sdb->n_part;j++){
			while(1){
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break;
				if(byte == 0) break;
				H.dir1 = (byte >> 1) & 0x01;
				H.dir2 = (byte >> 0) & 0x01;
				H.n_mm = (byte >> 2) & 0x0F;
				H.rid = 0;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.offset = byte;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) <<  0;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) <<  8;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) << 16;
				if(fread(&byte, 1, 1, aln_parts[j]) == 0) break; H.rid |= ((uint32_t)byte) << 24;
				push_hitv(hits2, H);
			}
		}
		for(j=0;j<r;j++){
			hit = ref_sr_hitv(aux2->hits, j);
			H.dir1   = hit->dir1;
			H.dir2   = hit->dir2;
			H.offset = hit->off;
			H.n_mm   = hit->n_mm;
			H.rid    = hit->rid2;
			push_hitv(hits2, H);
		}
		if(count_hitv(hits2) > mall->lmt / 2) set_hitv_size(hits2, mall->lmt / 2);
		if((i & 0xFF) == 0) renew_lgraph(g); 
		else reset_lgraph(g);
		rd_off1 = sr_rdseq_offset(sdb, 2 * i);
		rd_off2 = sr_rdseq_offset(sdb, 2 * i + 1);
		rd_len1 = sr_rdseq_length(sdb, 2 * i);
		rd_len2 = sr_rdseq_length(sdb, 2 * i + 1);
		bits2seq(seq1, ref_u64list(sdb->rd_seqs, 0), rd_off1, rd_len1);
		bits2revseq(seq2, ref_u64list(sdb->rd_seqs, 0), rd_off2, rd_len2);
		lib_id = get_read_lib(libs, 2 * i);
		min_ins = get_u32list(libs->min_ins, lib_id);
		max_ins = get_u32list(libs->max_ins, lib_id);
		push_lgraph(g, 2 * i + 0, seq1, rd_len1, 0, 0, 0, 0);
		push_lgraph(g, 2 * i + 1, seq2, rd_len2, 1, 1, 0, 0);
		l = count_hitv(hits1);
		r = count_hitv(hits2);
		if(alnf) fprintf(alnf, "A\t%u\t%u\t%u\t%s\t%s\n", i, l, r, seq1, seq2);
		for(j=0;j<l+r;j++){
			if(j < l){
				H = get_hitv(hits1, j);
			} else {
				H = get_hitv(hits2, j - l);
			}
			lib_id = get_read_lib(libs, H.rid);
			min_ins = get_u32list(libs->min_ins, lib_id);
			max_ins = get_u32list(libs->max_ins, lib_id);
			//dir = H.dir1 ^ H.dir2; // H.dir1 always equals H.dir2
			//dir ^= (j >= l); //hit->rid1 & 0x01;
			//dir ^= H.rid & 0x01;
			//at = (j >= l);
			rd_off1 = sr_rdseq_offset(sdb, H.rid);
			rd_len1 = sr_rdseq_length(sdb, H.rid);
			rd_off2 = sr_rdseq_offset(sdb, H.rid ^ 0x1U);
			rd_len2 = sr_rdseq_length(sdb, H.rid ^ 0x1U);
			off = H.dir1? - ((int)H.offset) : (int)H.offset;
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
			if(alnf) fprintf(alnf, "%c\t%u\t%d\t%d\t%d\t%u\t%s\t%s\n", "LR"[j<l? 0 : 1], H.rid, off, min_ins, max_ins, H.n_mm, seq1, seq2);
		}
		align_lgraph(g);
		if(dotf){ simplify_lgraph(g); print_dot_lgraph(g, dotf); }
		if(layout_lgraph(g)){
			skeleton_lgraph(g);
			align_skeleton_lgraph(g, sdb, g->aux);
			if(sdb->n_part > 1){
				dump_suppv(g->supps, cns_parts[sdb->part_idx]);
				fwrite(&SP, sizeof(Supp), 1, cns_parts[sdb->part_idx]);
				fwrite(&g->cns_id, sizeof(uint32_t), 1, skef);
				fwrite(&g->cns_len, sizeof(uint32_t), 1, skef);
				fwrite(&g->n_sr, sizeof(uint32_t), 1, skef);
				fwrite(&g->n_layout, sizeof(uint32_t), 1, skef);
				fwrite(g->skeleton->string, 1, g->skeleton->size, skef);
			} else {
				consensus_lgraph(g);
				output_lgraph(g, cnsf, msaf);
				n ++;
				if((n % 100) == 0){
					thread_beg_syn(mall);
					fflush(cnsf); fflush(msaf);
					fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
					fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
					fseek(cnsf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
					if(alnf){ fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET); }
					if(dotf){ fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET); }
					thread_end_syn(mall);
				}
			}
		}
	}
	t2 = time(NULL);
	fprintf(log, "[THREAD%03u] %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
	fflush(log);
	if(1){
		thread_beg_syn(mall);
		fflush(cnsf); fflush(msaf);
		fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
		fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
		fseek(cnsf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
		if(alnf){ fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET); }
		if(dotf){ fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET); }
		thread_end_syn(mall);
	}
	if(sdb->n_part > 1){
		fflush(skef);
		fflush(cns_parts[sdb->part_idx]);
		for(i=0;i+1<sdb->n_part;i++){
			fclose(aln_parts[i]);
			sprintf(chs, "%s.aln.p%03u.t%03u", mall->prefix, i, t_idx);
			unlink(chs);
		}
	}
} else if(mall->task == 4){
	t1 = time(NULL);
	sprintf(chs, "%s.cns.p%03u.t%03u", mall->prefix, sdb->part_idx, t_idx);
	if((cns_parts[sdb->part_idx] = fopen(chs, "w+")) == NULL){
		fprintf(stdout, " -- Cannot create %s in %s -- %s:%d --\n", chs, __FUNCTION__, __FILE__, __LINE__); fflush(stdout); abort();
	}
	m = 0;
	fseek(skef, 0, SEEK_SET);
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] aligning (long) %u, %u sec. \n", t_idx, m, (unsigned)(t2 - t1));
			fflush(log);
			renew_lgraph(g);
		} else reset_lgraph(g);
		fread(&g->cns_id, sizeof(uint32_t), 1, skef);
		fread(&g->cns_len, sizeof(uint32_t), 1, skef);
		fread(&g->n_sr, sizeof(uint32_t), 1, skef);
		fread(&g->n_layout, sizeof(uint32_t), 1, skef);
		encap_string(g->skeleton, g->cns_len);
		fread(g->skeleton->string, 1, g->cns_len, skef);
		g->skeleton->string[g->cns_len] = 0;
		g->skeleton->size = g->cns_len;
		align_skeleton_lgraph(g, sdb, g->aux);
		dump_suppv(g->supps, cns_parts[sdb->part_idx]);
		fwrite(&SP, sizeof(Supp), 1, cns_parts[sdb->part_idx]);
	}
} else if(mall->task == 5){
	t1 = time(NULL);
	fseek(skef, 0, SEEK_SET);
	n = m = 0;
	for(i=mall->rid_beg+t_idx;i<mall->rid_end;i+=n_cpu){
		m ++;
		if((m % 1000) == 0){
			t2 = time(NULL);
			fprintf(log, "[THREAD%03u] consensus %u / %u, %u sec. \n", t_idx, n, m, (unsigned)(t2 - t1));
			fflush(log);
			renew_lgraph(g);
		} else reset_lgraph(g);
		fread(&g->cns_id, sizeof(uint32_t), 1, skef);
		fread(&g->cns_len, sizeof(uint32_t), 1, skef);
		fread(&g->n_sr, sizeof(uint32_t), 1, skef);
		fread(&g->n_layout, sizeof(uint32_t), 1, skef);
		encap_string(g->skeleton, g->cns_len);
		fread(g->skeleton->string, 1, g->cns_len, skef);
		g->skeleton->string[g->cns_len] = 0;
		g->skeleton->size = g->cns_len;
		align_skeleton_lgraph(g, sdb, g->aux);
		for(j=0;j+1<sdb->n_part;j++){
			l = (j + sdb->n_part - 1) % sdb->n_part;
			while(fread(&SP, sizeof(Supp), 1, cns_parts[l]) == 1){ if(SP.len == 0) break; push_suppv(g->supps, SP); }
		}
		consensus_lgraph(g);
		output_lgraph(g, cnsf, msaf);
		n ++;
		if((n % 100) == 0){
			thread_beg_syn(mall);
			fflush(cnsf); fflush(msaf);
			fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
			fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
			fseek(cnsf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
			if(alnf){ fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET); }
			if(dotf){ fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET); }
			thread_end_syn(mall);
		}
	}
	if(1){
		thread_beg_syn(mall);
		fflush(cnsf); fflush(msaf);
		fwrite(cns_buf, ftell(cnsf), 1, mall->cnsf);
		fwrite(msa_buf, ftell(msaf), 1, mall->msaf);
		fseek(cnsf, 0, SEEK_SET); fseek(msaf, 0, SEEK_SET);
		if(alnf){ fflush(alnf); fwrite(aln_buf, ftell(alnf), 1, mall->alnf); fseek(alnf, 0, SEEK_SET); }
		if(dotf){ fflush(dotf); fwrite(dot_buf, ftell(dotf), 1, mall->dotf); fseek(dotf, 0, SEEK_SET); }
		thread_end_syn(mall);
	}
	sprintf(chs, "%s.cns.ske.t%03u", mall->prefix, t_idx);
	unlink(chs);
	for(j=0;j+1<sdb->n_part;j++){
		l = (j + sdb->n_part - 1) % sdb->n_part;
		fclose(cns_parts[l]);
		sprintf(chs, "%s.cns.p%03u.t%03u", mall->prefix, l, t_idx);
		unlink(chs);
	}
} else {
	microsleep(1);
}
thread_end_loop(mall);
sr_free_aux(aux1);
sr_free_aux(aux2);
free_hitv(hits1);
free_hitv(hits2);
fclose(cnsf); fclose(msaf);
if(alnf) fclose(alnf);
if(cns_buf) free(cns_buf);
if(msa_buf) free(msa_buf);
if(aln_buf) free(aln_buf);
if(sdb->n_part > 1){ free(aln_parts); free(cns_parts); }
thread_end_func(mall);

int debug_assembling(ATOptions opt, char *prefix){
	FileReader *fr;
	LGraph *g;
	FILE *cnsf, *msaf;
	uint32_t i, id, n_lr, rd_len1, rd_len2, m, p;
	int n, off, max_ins, min_ins;
	char *seq1, *seq2;
	time_t t1, t2;
	char cmd[256];
	if((fr = fopen_filereader(opt.load_aln)) == NULL){
		fprintf(stderr, " -- Cannot open '%s' in %s -- %s:%d --\n", opt.load_aln, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		abort();
	}
	g = init_lgraph(opt.kmer_size[1], opt.rd_len, opt.min_ol[1], opt.min_sm[1], opt.max_mm[1], opt.min_ins, opt.max_ins);
	sprintf(cmd, "%s.fasta", prefix);
	if((cnsf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "%s.msa", prefix);
	if((msaf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	m = p = 0;
	t1 = time(NULL);
	while(1){
		if((n = fread_table(fr)) == -1) break;
		m ++;
		if((m % 100) == 0){
			t2 = time(NULL);
			fprintf(stderr, "[%s] %u / %u clusters, %u sec. \n", __FUNCTION__, p, m, (unsigned)(t2 - t1)); fflush(stderr);
			renew_lgraph(g);
		} else reset_lgraph(g);
		id = atoll(get_col_str(fr, 1));
		n_lr = atoi(get_col_str(fr, 2)) + atoi(get_col_str(fr, 3));
		push_lgraph(g, 2 * id + 0, get_col_str(fr, 4), get_col_len(fr, 4), 0, 0, 0, 0);
		push_lgraph(g, 2 * id + 1, get_col_str(fr, 5), get_col_len(fr, 5), 1, 1, 0, 0);
		for(i=0;i<n_lr;i++){
			fread_table(fr);
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
		align_lgraph(g);
		if(layout_lgraph(g)){
			p ++;
			skeleton_lgraph(g);
			//align_skeleton_lgraph(g, sdb, g->aux);
			consensus_lgraph(g);
			output_lgraph(g, cnsf, msaf);
		}
	}
	free_lgraph(g);
	fclose(cnsf);
	fclose(msaf);
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
	Sequence *seq1, *seq2;
	FILE *log, *cnsf, *msaf, *alnf, *dotf;
	char *prefix, cmd[256], *date;
	uint32_t i, rid, rrid, n_vis, p;
	uint16_t *n_hits;
	int ii, c, vis;
	time_t tm;
	thread_preprocess(mall);
	if((c = parse_options(&opt, 0, argc, argv)) == -1) return usage_all();
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
	if(opt.min_lib == 0) opt.min_lib = ref_sfv(sfs, 0)->ins_size;
	for(i=0;2*i<count_sfv(sfs);i++){
		if(strcasecmp(ref_sfv(sfs, 2 * i)->prefix, ref_sfv(sfs, 2 * i + 1)->prefix)
			|| ref_sfv(sfs, 2 * i)->ins_size != ref_sfv(sfs, 2 * i + 1)->ins_size
			|| ref_sfv(sfs, 2 * i)->pair_idx >= ref_sfv(sfs, 2 * i + 1)->pair_idx){
			fprintf(stderr, " -- Wrong pair \"%s\" and \"%s\" in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	libs = guess_lib_inserts(sfs);
	sdb = sr_init_sdb(opt.kmer_size[0], opt.n_seed[0], opt.rd_len);
	sr_set_align_parameters(sdb, 1, opt.min_ol[0], opt.min_sm[0], opt.max_mm[0], 0);
	sr_set_filter_parameters(sdb, opt.low_cpx, opt.limit, 0);
	sr_set_n_part(sdb, opt.n_part);
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
	sprintf(cmd, "%s.msa", prefix);
	if((msaf = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	if(opt.dump_aln){
		sprintf(cmd, "gzip -1 -c >%s.aln.gz", prefix);
		if((alnf = popen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			return 1;
		}
		sprintf(cmd, "gzip -1 -c >%s.dot.gz", prefix);
		if((dotf = popen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
			fflush(stderr);
			return 1;
		}
	} else { alnf = NULL; dotf = NULL; }
	rid = 0;
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	for(i=0;i<libs->n_lib;i++){
		fprintf(log, "Lib%u: %u - %u bp\n", i, get_u32list(libs->min_ins, i), get_u32list(libs->max_ins, i));
	}
	fflush(log);
	n_vis = 0;
	for(i=0;2*i<count_sfv(sfs);i++){
		set_u32list(libs->lib_offs, i, rid * 2);
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
	if(opt.inc_ar == 0) sr_set_n_ar(sdb, 2 * n_vis);
	n_hits = malloc(sizeof(uint16_t) * 2 * n_vis);
	memset(n_hits, 0, sizeof(uint16_t) * 2 * n_vis);
	thread_beg_init(mall, (int)opt.n_cpu);
	mall->sdb = sdb;
	mall->g = init_lgraph(opt.kmer_size[1], opt.rd_len, opt.min_ol[1], opt.min_sm[1], opt.max_mm[1], opt.min_ins, opt.max_ins);
	mall->libs = libs;
	mall->rid_beg = 0;
	mall->rid_end = n_vis;
	mall->log  = log;
	mall->cnsf = cnsf;
	mall->msaf = msaf;
	mall->alnf = alnf;
	mall->dotf = dotf;
	mall->lmt  = opt.limit;
	mall->n_hits = n_hits;
	mall->task    = 0;
	mall->prefix  = prefix;
	thread_end_init(mall);
	if(opt.n_part > 1){
		for(p=0;p<opt.n_part;p++){
			sr_set_sr_part(sdb, p);
			fprintf(log, "Indexing [%u/%u]\n", p, opt.n_part);
			fflush(log);
			thread_beg_iter(mall);
			mall->task = 1; // index
			thread_wake(mall);
			thread_end_iter(mall);
			thread_waitfor_all_idle(mall);
			tm = time(NULL); date = asctime(localtime(&tm));
			fprintf(log, "%s\n", date);
			if(p + 1 == opt.n_part){
				fprintf(log, "Local assembling\n");
			} else {
				fprintf(log, "Writting temporary alignments\n");
			}
			fflush(log);
			thread_beg_iter(mall);
			mall->task = (p + 1 < opt.n_part)? 2 : 3; // align and asm
			thread_wake(mall);
			thread_end_iter(mall);
			thread_waitfor_all_idle(mall);
		}
		for(p=0;p+1<opt.n_part;p++){
			sr_set_sr_part(sdb, p);
			fprintf(log, "Indexing [%u/%u]\n", p, opt.n_part);
			fflush(log);
			thread_beg_iter(mall);
			mall->task = 1; // index
			thread_wake(mall);
			thread_end_iter(mall);
			thread_waitfor_all_idle(mall);
			tm = time(NULL); date = asctime(localtime(&tm));
			fprintf(log, "%s\n", date);
			if(p + 2 == opt.n_part){
				fprintf(log, "Consensus");
			} else {
				fprintf(log, "Aligning long sequence for consensus\n");
			}
			fflush(log);
			thread_beg_iter(mall);
			mall->task = (p + 2 < opt.n_part)? 4 : 5;
			thread_wake(mall);
			thread_end_iter(mall);
			thread_waitfor_all_idle(mall);
		}
	} else {
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
		mall->task = 3; // align and asm
		thread_wake(mall);
		thread_end_iter(mall);
		thread_waitfor_all_idle(mall);
	}
	thread_beg_close(mall);
	free_lgraph(mall->g);
	thread_end_close(mall);
	sr_free_sdb(sdb);
	free(n_hits);
	free_lib_inserts(libs);
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	tm = time(NULL); date = asctime(localtime(&tm));
	fprintf(log, "%s\n", date);
	fflush(log);
	fclose(log);
	fclose(cnsf);
	fclose(msaf);
	if(alnf) pclose(alnf);
	return 0;
}


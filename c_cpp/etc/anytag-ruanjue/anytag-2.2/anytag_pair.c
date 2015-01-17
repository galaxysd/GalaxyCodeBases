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
#include <stdint.h>
#include "anytag_aux.h"
#include "file_reader.h"
#include "sr_aln.h"
#include "hashset.h"
#include "list.h"
#include "string.h"
#include "thread.h"

int usage_pair();

define_list_core(u32v, uint32_t, uint16_t, 8);
define_list(u32vv, u32v*);

typedef struct {
	uint64_t id1:32, id2:32, dir1:1, dir2:1, closed:1, cnt:29;
} Link;

#define link_hashcode(link) u64hashcode((((uint64_t)(link).dir1) << 63) | (((uint64_t)(link).id1) << 32) | (((uint64_t)(link).dir2) << 31) | ((link).id2))
#define link_equals(lnk1, lnk2) (((lnk1).id1 == (lnk2).id1) && ((lnk1).id2 == (lnk2).id2) && ((lnk1).dir1 == (lnk2).dir1) && ((lnk1).dir2 == (lnk2).dir2))
define_hashset(lnkhash, Link, link_hashcode, link_equals);

define_list(linkv, Link);

typedef struct {
	uint64_t nameoff;
	uint32_t namelen:31;
	uint64_t seqoff;
	uint32_t seqlen;
	uint32_t closed:1;
} FIS;

define_list(fisv, FIS);

typedef struct {
	fisv  *fiss;
	u64list *seqs;
	uint64_t offset;
	String *names;
	lnkhash **hashs;
	linkv  *links;
	u32vv *hits;
	u8list *mms;
	FILE *out;
} MatePool;

int cmp_fis_link_func(const void *e1, const void *e2){
	return ((int)((Link*)e2)->cnt) - ((int)((Link*)e1)->cnt);
}

thread_beg_def(mpair);
MatePool *mp;
SR_SeqDB *sdb;
int task;
thread_end_def(mpair);

thread_beg_func(mpair);
uint32_t i, j, k, t_idx, n_cpu, dir;
uint32_t id1, id2, dir1, dir2, n, n_rd, mm;
int exists;
u8list *mms;
FIS *fis,*fis1, *fis2;
u32v *hits, *mate1, *mate2;
SR_SeqDB *sdb;
SR_AlnAux *aux;
SR_AlnHit *hit;
MatePool *mp;
Link LNK, *lnk, *lnk2;
lnkhash *hash;
sdb   = mpair->sdb;
n_rd  = sdb->n_rd;
mp    = mpair->mp;
t_idx = mpair->t_idx;
n_cpu = mpair->n_cpu;
mp->hashs[t_idx] = init_lnkhash(1023);
hash  = mp->hashs[t_idx];
aux   = sr_init_aux();
mms   = init_u8list(1024);
append_u8list(mms, mp->mms);
thread_beg_loop(mpair);
if(mpair->task == 1){ // indexing
	for(i=t_idx;i<sdb->n_idx;i+=n_cpu) sr_index_sdb(sdb, i);
} else if(mpair->task == 2){ // aligning
	n = 0;
	clear_sr_hitv(aux->hits);
	for(i=t_idx;i<count_fisv(mp->fiss);i+=n_cpu){
		fis = ref_fisv(mp->fiss, i);
		sr_align_long_sdb(sdb, aux, i, ref_u64list(mp->seqs, 0), fis->seqoff, fis->seqlen);
		if(count_sr_hitv(aux->hits) >= 16 * 1024){
			thread_beg_syn(mpair);
			for(j=0;j<count_sr_hitv(aux->hits);j++){
				hit = ref_sr_hitv(aux->hits, j);
				if(hit->n_mm > get_u8list(mms, hit->rid2)) continue;
				set_u8list(mms, hit->rid2, hit->n_mm);
				mm = hit->n_mm;
				mm <<= 28;
				dir = (hit->dir1 ^ hit->dir2)? 0x80000000U : 0x00000000U;
				hits = get_u32vv(mp->hits, hit->rid2);
				push_u32v(hits, hit->rid1 | mm | dir);
				n ++;
			}
			thread_end_syn(mpair);
			clear_sr_hitv(aux->hits);
		}
	}
	{
			thread_beg_syn(mpair);
			for(j=0;j<count_sr_hitv(aux->hits);j++){
				hit = ref_sr_hitv(aux->hits, j);
				if(hit->n_mm > get_u8list(mms, hit->rid2)) continue;
				set_u8list(mms, hit->rid2, hit->n_mm);
				mm = hit->n_mm;
				mm <<= 28;
				dir = (hit->dir1 ^ hit->dir2)? 0x80000000U : 0x00000000U;
				hits = get_u32vv(mp->hits, hit->rid2);
				push_u32v(hits, hit->rid1 | mm | dir);
				n ++;
			}
			thread_end_syn(mpair);
	}
	fprintf(stderr, "Thread[%u] %u hits\n", t_idx, n);
	fflush(stderr);
} else if(mpair->task == 3){ // One by one
	for(j=0;j<count_u8list(mms);j++){
		if(get_u8list(mms, j) < get_u8list(mp->mms, j)) set_u8list(mp->mms, j, get_u8list(mms, j));
	}
	free_u8list(mms);
} else if(mpair->task == 4){
	LNK.cnt = 1;
	LNK.closed = 0;
	if(t_idx == 0){
		n = 0;
		for(i=0;i<n_rd;i++){
			mate1 = get_u32vv(mp->hits, i);
			for(j=0;j<count_u32v(mate1);j++){
				id1 = get_u32v(mate1, j);
				mm = (id1 >> 28) & 0x7;
				id1 = id1 & 0x0FFFFFFFU;
				if(mm > get_u8list(mp->mms, id1)) continue;
				n ++;
			}
		}
		fprintf(stderr, "Total %u best hits\n", n);
		fflush(stderr);
	}
	n = 0;
	for(i=0;(2*i+1)<n_rd;i+=2){
		mate1 = get_u32vv(mp->hits, 2 * i + 0);
		mate2 = get_u32vv(mp->hits, 2 * i + 1);
		for(j=0;j<count_u32v(mate1);j++){
			id1 = get_u32v(mate1, j);
			mm = (id1 >> 28) & 0x7;
			dir1 = id1 >> 31;
			id1 = id1 & 0x0FFFFFFFU;
			if(mm > get_u8list(mp->mms, id1)) continue;
			for(k=0;k<count_u32v(mate2);k++){
				id2 = get_u32v(mate2, k);
				mm = (id2 >> 28) & 0x7;
				dir2 = id2 >> 31;
				id2 = id2 & 0x0FFFFFFFU;
				if(mm > get_u8list(mp->mms, id2)) continue;
				if(((id1 + id2)) % n_cpu != t_idx) continue;
				if(id1 == id2) continue;
				n ++;
				if(id1 < id2){
					LNK.id1  = id1;
					LNK.dir1 = dir1;
					LNK.id2  = id2;
					LNK.dir2 = dir2;
				} else {
					LNK.id1  = id2;
					LNK.dir1 = !dir2;
					LNK.id2  = id1;
					LNK.dir2 = !dir1;
				}
				lnk = prepare_lnkhash(hash, LNK, &exists);
				if(exists){
					lnk->cnt ++;
				} else {
					*lnk = LNK;
				}
			}
		}
	}
	fprintf(stderr, "Thread[%u] %u / %u links\n", t_idx, (unsigned)count_lnkhash(hash), n);
	fflush(stderr);
} else if(mpair->task == 5){
	if(1) continue;
} else if(mpair->task == 6){ // One by one
	while((lnk2 = ref_iter_lnkhash(hash))){
		push_linkv(mp->links, *lnk2);
	}
	free_lnkhash(hash);
	fprintf(stderr, "Thread[%u] %u links\n", t_idx, (unsigned)count_linkv(mp->links));
	fflush(stderr);
} else if(mpair->task == 7){
	if(t_idx != 0) continue;
	qsort(ref_linkv(mp->links, 0), count_linkv(mp->links), sizeof(Link), cmp_fis_link_func);
} else if(mpair->task == 8){
	if(t_idx != 0) continue;
	for(i=0;i<count_linkv(mp->links);i++){
		lnk = ref_linkv(mp->links, i);
		fis1 = ref_fisv(mp->fiss, lnk->id1);
		fis2 = ref_fisv(mp->fiss, lnk->id2);
		//if(fis1->closed || fis2->closed) continue;
		fprintf(mp->out, "%s\t%c\t%s\t%c\t%u\n", mp->names->string + fis1->nameoff, "+-"[lnk->dir1],
				mp->names->string + fis2->nameoff, "+-"[lnk->dir2], lnk->cnt);
		fis1->closed = 1;
		fis2->closed = 1;
	}
} else {
	microsleep(1);
}
thread_end_loop(mpair);
sr_free_aux(aux);
thread_end_func(mpair);

int main_pair(int argc, char **argv){
	SR_SeqDB *sdb;
	MatePool *mp;
	FIS *fis;
	Sequence *read1, *read2;
	FileReader *fr1, *fr2, *lr;
	FILE *out;
	float opt_s;
	uint32_t opt_x, opt_c, opt_t, opt_r, opt_n, opt_k, opt_m, opt_g;
	char *opt_i, *opt_o, *opt_1, *opt_2;
	int is_fq1, is_fq2, c, out_type;
	uint32_t i, n;
	char cmd[256];
	thread_preprocess(mpair);
	{
	opt_x = 0;
	opt_i = NULL;
	opt_o = NULL;
	opt_1 = NULL;
	opt_2 = NULL;
	opt_c = 10;
	opt_t = 1;
	opt_r = 0;
	opt_n = 3;
	opt_k = 15;
	opt_s = 0.95;
	opt_m = 6;
	opt_g = 0;
	is_fq1 = 1;
	is_fq2 = 1;
	}
	while((c = getopt(argc, argv, "hgx:i:o:1:2:c:t:r:n:k:s:m:")) != -1){
		switch(c){
			case 'h': return usage_pair();
			case 'g': opt_g = 1; break;
			case 'x': opt_x = atoi(optarg); break;
			case 'c': opt_c = atoi(optarg); break;
			case 't': opt_t = atoi(optarg); break;
			case 'r': opt_r = atoi(optarg); break;
			case 'n': opt_n = atoi(optarg); break;
			case 'k': opt_k = atoi(optarg); break;
			case 'm': opt_m = atoi(optarg); break;
			case 'i': opt_i = optarg; break;
			case 'o': opt_o = optarg; break;
			case '1': opt_1 = optarg; break;
			case '2': opt_2 = optarg; break;
			default: fprintf(stdout, "unknown option -%c\n", c); return usage_pair();
		}
	}
	if(opt_m > 7) opt_m = 7;
	if(opt_x == 0 || opt_i == NULL || opt_1 == NULL || opt_2 == NULL) return usage_pair();
	if((lr = fopen_filereader(opt_i)) == NULL){ fprintf(stdout, "Cannot open %s\n", opt_i); abort(); }
	if(strlen(opt_1) > 3 && strcasecmp(opt_1 + strlen(opt_1) - 3, ".fa") == 0) is_fq1 = 0;
	if(strlen(opt_1) > 6 && strcasecmp(opt_1 + strlen(opt_1) - 6, ".fa.gz") == 0) is_fq1 = 0;
	if(strlen(opt_2) > 3 && strcasecmp(opt_2 + strlen(opt_2) - 3, ".fa") == 0) is_fq2 = 0;
	if(strlen(opt_2) > 6 && strcasecmp(opt_2 + strlen(opt_2) - 6, ".fa.gz") == 0) is_fq2 = 0;
	if((fr1 = fopen_filereader(opt_1)) == NULL){ fprintf(stdout, "Cannot open %s\n", opt_1); abort(); }
	if((fr2 = fopen_filereader(opt_2)) == NULL){ fprintf(stdout, "Cannot open %s\n", opt_2); abort(); }
	mp = malloc(sizeof(MatePool));
	mp->fiss = init_fisv(1024);
	mp->hits = init_u32vv(1024);
	mp->seqs = init_u64list(1024);
	mp->offset = 64;
	push_u64list(mp->seqs, 0);
	push_u64list(mp->seqs, 0);
	mp->names = init_string(1024);
	mp->links = init_linkv(1024);
	mp->hashs = malloc(sizeof(lnkhash*) * opt_t);
	mp->mms   = init_u8list(1024);
	sdb = sr_init_sdb(opt_k, opt_n, opt_r);
	sr_set_align_parameters(sdb, 1 | 2, opt_r, opt_s, opt_m, opt_g);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Started\n");
	fflush(stderr);
	read1 = read2 = NULL;
	while(fread_fasta(&read1, lr)){
		fis = next_ref_fisv(mp->fiss);
		fis->seqlen = read1->seq.size;
		fis->seqoff = mp->offset;
		fis->closed = 0;
		mp->offset += fis->seqlen;
		encap_u64list(mp->seqs, (mp->offset + 31 + 64) / 32);
		seq2bits(ref_u64list(mp->seqs, 0), fis->seqoff, read1->seq.string, fis->seqlen);
		fis->nameoff = mp->names->size;
		fis->namelen = read1->name.size;
		append_string(mp->names, read1->name.string, fis->namelen + 1);
	}
	fclose_filereader(lr);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Loaded %u FIS\n", (unsigned)count_fisv(mp->fiss));
	fflush(stderr);
	if(count_fisv(mp->fiss) > 0x0FFFFFFFU){
		fprintf(stderr, "Too many FIS %u > %u\n", (unsigned)count_fisv(mp->fiss), 0x0FFFFFFFU);
		fflush(stderr);
		abort();
	}
	while(1){
		if(!(is_fq1? fread_fastq_adv(&read1, fr1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&read1, fr1, FASTA_FLAG_NO_NAME))) break;
		if(!(is_fq2? fread_fastq_adv(&read2, fr2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&read2, fr2, FASTA_FLAG_NO_NAME))) break;
		sr_push_sdb(sdb, read1->seq.string, read1->seq.size);
		sr_push_sdb(sdb, read2->seq.string, read2->seq.size);
		push_u32vv(mp->hits, init_u32v(4));
		push_u32vv(mp->hits, init_u32v(4));
		push_u8list(mp->mms, opt_m + 1);
		push_u8list(mp->mms, opt_m + 1);
	}
	fclose_filereader(fr1);
	fclose_filereader(fr2);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Loaded %u mate pairs\n", (unsigned)count_u32vv(mp->hits) / 2);
	fflush(stderr);
	if(opt_o == NULL){
		out = stdout;
		out_type = 0;
	} else if(strlen(opt_o) > 3 && strcasecmp(opt_o + strlen(opt_o) - 3, ".gz") == 0){
		out_type = 1;
		sprintf(cmd, "gzip -1 -c >%s", opt_o);
		out = popen(cmd, "w");
	} else {
		out_type = 2;
		out = fopen(opt_o, "w");
	}
	mp->out = out;
	sr_ready_sdb(sdb);
	{
	thread_beg_init(mpair, (int)opt_t);
	mpair->mp = mp;
	mpair->sdb = sdb;
	mpair->task = 0;
	thread_end_init(mpair);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 1;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Indexed\n");
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 2;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Aligned\n");
	fflush(stderr);
	sr_free_sdb(sdb);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 3;
	thread_wake(mpair);
	thread_waitfor_idle(mpair);
	thread_end_iter(mpair);
	n = 0;
	for(i=0;i<count_u8list(mp->mms);i++){ if(get_u8list(mp->mms, i) <= opt_m) n ++; }
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Minimum mismatches found for %u matched reads\n", n);
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 4;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	for(i=0;i<count_u32vv(mp->hits);i++) free_u32v(get_u32vv(mp->hits, i));
	free_u32vv(mp->hits);
	mp->hits = NULL;
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Bundled\n");
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 5;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Merged\n");
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 6;
	thread_wake(mpair);
	thread_waitfor_idle(mpair);
	thread_end_iter(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Dumped hash->vector\n");
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 7;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Sorted\n");
	fflush(stderr);
	}
	{
	thread_beg_iter(mpair);
	mpair->task = 8;
	thread_wake(mpair);
	thread_end_iter(mpair);
	thread_waitfor_all_idle(mpair);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Paired\n");
	fflush(stderr);
	}
	thread_beg_close(mpair);
	thread_end_close(mpair);
	if(out_type == 1) pclose(out);
	else if(out_type == 2) fclose(out);
	free_u64list(mp->seqs);
	free_string(mp->names);
	free_fisv(mp->fiss);
	free_linkv(mp->links);
	free(mp);
	fprintf(stderr, "%s\n", date());
	fprintf(stderr, "Finished\n");
	fflush(stderr);
	return 0;
}


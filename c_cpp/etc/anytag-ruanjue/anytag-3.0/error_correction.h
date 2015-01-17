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
 
/*
 * Example:
 * sdb = init sdb, load sequences, and index.
 * ReadCorrector corr = init_corrector(sdb, 4, 0.25);
 * int rid = 7;
 * hits = short read align.
 * if(do_corrector(corr, 7, hits)){
 *   if(corr->corrs->size){
 *     fprintf("Check and Corrected: ");
 *     for(i=0;i<corr->corrs->size;i++){
 *       cr = ref_corrv(corr->corrs, i);
 *       printf(", %d(%c -> %c)", cr->pos, "ACGT"[cr->base[0]], "ACGT"[cr->base[1]]);
 *       printf("\n");
 *     }
 *   } else {
 *     printf("Checked but not corrected\n");
 *   }
 * } else {
 *   printf("No enough reads aligned to do correction\n");
 * }
 */
#ifndef __SHORT_READS_CORRECTION_RJ_H
#define __SHORT_READS_CORRECTION_RJ_H

#include "sr_aln.h"

typedef struct {
	uint32_t error_cnt, major_cnt;
	uint16_t pos;
	uint8_t  base[2];
} corr_t;

define_list(corrv, corr_t);

typedef struct {
	SR_SeqDB *sdb;
	uint32_t max_cnt;
	float    max_freq;
	uint8_t  rd_bases[256];
	uint32_t rd_len, msa_bases[256][4];
	corrv    *corrs;
} ReadCorrector;

static inline ReadCorrector* init_corrector(SR_SeqDB *sdb, uint32_t max_cnt, float max_freq){
	ReadCorrector *corr;
	corr = malloc(sizeof(ReadCorrector));
	corr->sdb = sdb;
	corr->max_cnt = max_cnt;
	corr->max_freq = max_freq;
	corr->corrs = init_corrv(12);
	return corr;
}

static inline void free_corrector(ReadCorrector *corr){
	free_corrv(corr->corrs);
	free(corr);
}

static inline int do_corrector(ReadCorrector *corr, uint32_t rid, sr_hitv *hits){
	SR_AlnHit *hit;
	corr_t *cr;
	uint64_t off;
	uint32_t i, j, len, b,e, dir, max_cnt, err_cnt;
	uint8_t c;
	int s;
	clear_corrv(corr->corrs);
	corr->rd_len = sr_rdseq_length(corr->sdb, rid);
	off = sr_rdseq_offset(corr->sdb, rid);
	for(i=0;i<corr->rd_len;i++){
		corr->rd_bases[i] = bits2bit(corr->sdb->rd_seqs->buffer, off + i);
		corr->msa_bases[i][0] = 0;
		corr->msa_bases[i][1] = 0;
		corr->msa_bases[i][2] = 0;
		corr->msa_bases[i][3] = 0;
	}
	if(hits->size * corr->max_freq <= 1) return 0;
	for(i=0;i<hits->size;i++){
		hit = ref_sr_hitv(hits, i);
		if(hit->is_gap) continue;
		if(hit->rid1 == rid){
			off = sr_rdseq_offset(corr->sdb, hit->rid2);
			len = sr_rdseq_length(corr->sdb, hit->rid2);
			dir = hit->dir1 ^ hit->dir2;
			if(hit->dir1){
				s = - hit->off;
			} else {
				s = hit->off;
			}
		} else {
			off = sr_rdseq_offset(corr->sdb, hit->rid1);
			len = sr_rdseq_length(corr->sdb, hit->rid1);
			dir = hit->dir1 ^ hit->dir2;
			if(hit->dir2){
				s = hit->off;
			} else {
				s = - hit->off;
			}
		}
		if(s < 0){
			b = - s;
		} else {
			b = 0;
		}
		if(s + (int)len < (int)corr->rd_len){
			e = len;
		} else {
			e = corr->rd_len - s;
		}
		if(dir){
			for(j=b;j<e;j++){
				c = bits2revbit(corr->sdb->rd_seqs->buffer, off + len - j - 1);
				corr->msa_bases[j + s][c] ++;
			}
		} else {
			for(j=b;j<e;j++){
				c = bits2bit(corr->sdb->rd_seqs->buffer, off + j);
				corr->msa_bases[j + s][c] ++;
			}
		}
	}
	for(i=0;i<corr->rd_len;i++){
		c = 0;
		for(j=1;j<4;j++) if(corr->msa_bases[i][j] > corr->msa_bases[i][c]) c = j;
		if(c == corr->rd_bases[i]) continue;
		max_cnt = corr->msa_bases[i][c];
		err_cnt = corr->msa_bases[i][corr->rd_bases[i]] + 1;
		if(err_cnt > corr->max_cnt) continue;
		if(max_cnt * corr->max_freq < err_cnt) continue;
		cr = next_ref_corrv(corr->corrs);
		cr->pos = i;
		cr->base[0] = corr->rd_bases[i];
		cr->base[1] = c;
		cr->error_cnt = err_cnt;
		cr->major_cnt = max_cnt;
		corr->rd_bases[i] = c;
	}
	return 1;
}

#endif

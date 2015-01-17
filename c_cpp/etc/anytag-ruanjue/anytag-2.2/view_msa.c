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

#include <stdio.h>
#include "list.h"
#include "dna.h"
#include "local_assembly.h"

int main(){
	u64list *bits;
	Supp *sp;
	uint32_t cns_id, cns_len, n_cns, min_bcov, min_pcov, n_sr, n_layout;
	uint32_t seq_offset, off, i, j, n, m;
	uint64_t head;
	char chs[2000], rdseq[500];
	bits = init_u64list(1024);
	while(fread(&head, sizeof(uint64_t), 1, stdin) == 1){
		n = head >> 32;
		cns_id = head & 0xFFFFFFFFU;
		clear_u64list(bits);
		push_u64list(bits, head);
		fread(bits->buffer + 1, sizeof(uint64_t), n - 1, stdin);
		cns_len = bits->buffer[1] >> 32;
		n_cns = bits->buffer[1] & 0xFFFFFFFFU;
		min_bcov = bits->buffer[2] >> 32;
		min_pcov = bits->buffer[2] & 0xFFFFFFFFU;
		n_sr = bits->buffer[3] >> 32;
		n_layout = bits->buffer[3] & 0xFFFFFFFFU;
		seq_offset = 4 + n_cns;
		fprintf(stdout, ">ps%u len=%u n_rd=%u,%u,%u min_base_cov=%u min_pair_cov=%u\n", cns_id, cns_len, n_cns, n_layout, n_sr, min_bcov, min_pcov);
		bits2seq(chs, bits->buffer + seq_offset, 0, cns_len);
		off = cns_len;
		fprintf(stdout, "ps%08u\t%s\n", cns_id, chs);
		for(i=0;i<n_cns;i++){
			sp = (Supp*)(bits->buffer + 4 + i);
			bits2seq(rdseq, bits->buffer + seq_offset, off, sp->len);
			m = 0;
			for(j=0;j<sp->len;j++){
				if(rdseq[j] != chs[sp->off + j]){
					rdseq[j] += 'a' - 'A';
					m ++;
				}
			}
			off += sp->len;
			fprintf(stdout, "%010u\t", sp->rid);
			for(j=0;j<sp->off;j++) fputc('-', stdout);
			fprintf(stdout, "%s", rdseq);
			for(j=sp->off+sp->len;j<cns_len;j++) fputc('-', stdout);
			fprintf(stdout, "\t%c\t%u\t%u\n", "+-"[sp->dir], sp->off, m);
		}
	}
	free_u64list(bits);
	return 0;
}

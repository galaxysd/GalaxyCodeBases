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
 
#include "all_path.h"

int cmp_trace(const void *e1, const void *e2, void *ref){
	Trace *t1, *t2;
	t1 = (Trace*)e1;
	t2 = (Trace*)e2;
	if(t1->h_off == t2->h_off) return 0;
	if(t1->h_off > t2->h_off) return 1;
	else return -1;
	ref = ref;
}

MyPath* init_mypath(uint32_t n_cpu, uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm,  int allow_gap, uint32_t min_ins, uint32_t max_ins, uint32_t max_n_rd){
	MyPath *mp;
	mp = malloc(sizeof(MyPath));
	mp->sdb    = sr_init_sdb(NULL, n_cpu, kmer_size, rd_len);
	sr_set_align_parameters(mp->sdb, 1, min_ol, min_sm, max_mm, allow_gap);
	sr_set_filter_parameters(mp->sdb, 0, 1024);
	mp->sdb->gap_init = 12;
	mp->nodes  = init_nodev(12);
	mp->path   = init_u32v(12);
	mp->cigars = init_cigarvv(12);
	mp->cns    = init_string(1024);
	mp->comps  = init_compv(1024);
	mp->matrix = init_strv(12);
	mp->begs   = init_u32list(1024);
	mp->ends   = init_u32list(1024);
	mp->fwds   = init_u32list(12);
	mp->revs   = init_u32list(12);
	mp->traces = init_tracev(1024);
	mp->heap   = init_heap(cmp_trace, mp);
	mp->pair_score = 0;
	mp->cns_len  = 0;
	mp->p1     = 0;
	mp->p2     = 1;
	mp->min_ins  = min_ins;
	mp->max_ins  = max_ins;
	mp->max_n_rd = max_n_rd;
	mp->sn     = 0;
	mp->exec_simp  = 0;
	mp->output_aln = 0;
	mp->output_dot = 0;
	mp->output_cns = 0;
	return mp;
}

void free_mypath(MyPath *mp){
	Node *n;
	uint32_t i;
	for(i=0;i<count_cigarvv(mp->cigars);i++){
		cigarv_free(ref_cigarvv(mp->cigars, i));
	}
	free_cigarvv(mp->cigars);
	for(i=0;i<count_nodev(mp->nodes);i++){
		n = ref_nodev(mp->nodes, i);
		free_u32v(n->links);
	}
	sr_free_sdb(mp->sdb);
	free_compv(mp->comps);
	free_nodev(mp->nodes);
	free_u32v(mp->path);
	free_string(mp->cns);
	free_u32list(mp->fwds);
	free_u32list(mp->revs);
	free_tracev(mp->traces);
	free_heap(mp->heap);
	for(i=0;i<count_strv(mp->matrix);i++) free_string(get_strv(mp->matrix, i));
	free_strv(mp->matrix);
	free_u32list(mp->begs);
	free_u32list(mp->ends);
	free(mp);
}

void reset_mypath(MyPath *mp){
	Node *n;
	uint32_t i;
	sr_reset_sdb(mp->sdb);
	for(i=0;i<count_nodev(mp->nodes);i++){
		n = ref_nodev(mp->nodes, i);
		free_u32v(n->links);
	}
	clear_nodev(mp->nodes);
	mp->pair_score = 0;
	mp->cns_len = 0;
	mp->p1 = 0;
	mp->p2 = 1;
	mp->sn ++;
}

void push_mypath(MyPath *mp, uint32_t rid, char *seq, uint8_t seqlen, int dir){
	Node *n;
	sr_push_sdb(mp->sdb, seq, seqlen, 1);
	n = next_ref_nodev(mp->nodes);
	n->rid   = rid;
	n->len   = seqlen;
	n->used  = 0;
	n->inasm = 0;
	n->dir   = dir;
	n->off   = 0;
	n->idx   = 0;
	n->prev_idx   = 0;
	n->links = init_u32v(2);
}

int cmp_sr_alnhit(const void *e1, const void *e2){
	SR_AlnHit *h1, *h2;
	h1 = (SR_AlnHit*)e1;
	h2 = (SR_AlnHit*)e2;
	if(h1->n_ol == h2->n_ol){
		if(h1->n_mm == h2->n_mm) return 0;
		else if(h1->n_mm < h2->n_mm) return 1;
		else return -1;
	} else if(h1->n_ol < h2->n_ol) return 1;
	else return -1;
}

void _output_dot_graph(MyPath *mp, FILE *out){
	Node *n;
	SR_AlnHit *hit;
	uint32_t i;
	char cstr[64];
	if(count_nodev(mp->nodes) == 0) return;
	fprintf(out, "digraph G%d {\n", ref_nodev(mp->nodes, mp->p1)->rid >> 1);
	for(i=0;i<count_nodev(mp->nodes);i++){
		n = ref_nodev(mp->nodes, i);
		fprintf(out, " N%03u\n", i);
		//fprintf(out, " N%03u\n", n->rid);
	}
	for(i=0;i<count_sr_hitv(mp->sdb->hits);i++){
		hit = ref_sr_hitv(mp->sdb->hits, i);
		if(hit->n_cigar == 0) continue;
		cigars2string(hit->cigars, hit->n_cigar, cstr);
		fprintf(out, " N%03d -> N%03d [label=\"%s\"]\n", hit->rid1, hit->rid2, cstr);
		//fprintf(out, " N%03d -> N%03d [label=\"%s\"]\n", ref_nodev(mp->nodes, hit->rid1)->rid, ref_nodev(mp->nodes, hit->rid2)->rid, cstr);
	}
	fprintf(out, "}\n");
	fflush(out);
}

int _n_overlapped(MyPath *mp, uint32_t nid1, uint32_t nid2){
	Node *n;
	SR_AlnHit *hit;
	uint32_t i;
	if(nid1 == nid2) return 0;
	n = ref_nodev(mp->nodes, nid1);
	for(i=0;i<count_u32v(n->links);i++){
		hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, i));
		if(hit->n_cigar == 0) continue;
		if(hit->rid2 == nid2) return (n->len - hit->off);
	}
	return 0;
}

void simplify_mypath_core(MyPath *mp, u32list *ids, uint32_t dir){
	SR_AlnHit *hit;
	uint32_t i, id1, id2, idx, ov;
	int j;
	for(i=count_u32list(ids)-2;i>0;i-=2){
		id1 = get_u32list(ids, i);
		idx = get_u32list(ids, i + 1);
		for(j=i-2;j>=0;j-=2){
			id2 = get_u32list(ids, j);
			ov = dir? _n_overlapped(mp, id1, id2) : _n_overlapped(mp, id2, id1);
			if(ov == 0) continue;
			hit = ref_sr_hitv(mp->sdb->hits, idx);
			hit->n_cigar = 0;
			break;
		}
	}
}

void simplify_mypath(MyPath *mp){
	Node *n;
	SR_AlnHit *hit;
	uint32_t i, j, dir, id;
	for(i=0;i<count_nodev(mp->nodes);i++){
		n = ref_nodev(mp->nodes, i);
		clear_u32list(mp->fwds);
		for(j=0;j<count_u32v(n->links);j++){
			hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, j));
			if(hit->n_cigar == 0) continue;
			if(hit->rid1 == i){
				dir = hit->dir1;
				id = hit->rid2;
			} else {
				dir = !hit->dir2;
				id = hit->rid1;
			}
			if(id == i) continue;
			if(dir){
			} else {
				push_u32list(mp->fwds, id);
				push_u32list(mp->fwds, get_u32v(n->links, j));
			}
		}
		if(count_u32list(mp->fwds)) simplify_mypath_core(mp, mp->fwds, 0);
	}
}

void aln_mypath(MyPath *mp){
	SR_AlnHit *hit;
	Node *n1;
	Node *n2;
	uint32_t i, t;
	if(count_nodev(mp->nodes) >= mp->max_n_rd){
		return;
	} else {
		sr_ready_sdb(mp->sdb);
		sr_aln_sdb(mp->sdb);
	}
	qsort(ref_sr_hitv(mp->sdb->hits, 0), count_sr_hitv(mp->sdb->hits), sizeof(SR_AlnHit), cmp_sr_alnhit);
	for(i=0;i<count_sr_hitv(mp->sdb->hits);i++){
		hit = ref_sr_hitv(mp->sdb->hits, i);
		if(mp->output_aln) sr_output_hit(mp->sdb, hit, stderr);
		if(hit->dir1){
			hit->dir1 = hit->dir2 = 0;
			t = hit->rid1; hit->rid1 = hit->rid2; hit->rid2 = t;
			n1 = ref_nodev(mp->nodes, hit->rid1);
			n2 = ref_nodev(mp->nodes, hit->rid2);
			if(n1->len + hit->off < n2->len) continue;
			hit->off = n1->len - n2->len + hit->off;
		} else {
			n1 = ref_nodev(mp->nodes, hit->rid1);
			n2 = ref_nodev(mp->nodes, hit->rid2);
		}
		push_u32v(n1->links, i);
		push_u32v(n2->links, i);
	}
	if(mp->exec_simp) simplify_mypath(mp);
	if(mp->output_dot) _output_dot_graph(mp, stderr);
	//fprintf(stderr, "\r%8u", ref_nodev(mp->nodes, mp->p1)->rid);
	//fflush(stderr);
}

void push_trace_mp(MyPath *mp, uint32_t p_off, uint32_t h_off, uint32_t nid1, uint32_t nid2, uint32_t idx){
	Trace *t;
	t = ref_next_tracev(mp->traces);
	t->p_off = p_off;
	t->h_off = h_off;
	t->nid1  = nid1;
	t->nid2  = nid2;
	t->idx   = idx;
	push_heap(mp->heap, t);
}

Trace* pop_trace_mp(MyPath *mp){
	return (Trace*)pop_heap(mp->heap);
}

int solve_mypath(MyPath *mp){
	Node *n, *nn;
	SR_AlnHit *hit;
	Trace *t;
	uint32_t nid, i;
	if(count_nodev(mp->nodes) == 0) return 0;
	if(count_nodev(mp->nodes) >= mp->max_n_rd) return 0;
	nid = mp->p1;
	n = ref_nodev(mp->nodes, nid);
	n->used = 1;
	n->off  = 0;
	n->prev = 0xFFFFFFFFU;
	clear_heap(mp->heap);
	clear_tracev(mp->traces);
	encap_tracev(mp->traces, count_sr_hitv(mp->sdb->hits));
	for(i=0;i<count_u32v(n->links);i++){
		hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, i));
		if(hit->n_cigar == 0) continue;
		if(hit->rid2 == nid) continue;
		push_trace_mp(mp, n->off, hit->off, nid,  hit->rid2, i);
	}
	while((t = pop_trace_mp(mp))){
		n  = ref_nodev(mp->nodes, t->nid2);
		if(n->used) continue;
		n->prev_idx = t->idx;
		n->used     = 1;
		n->prev     = t->nid1;
		n->off      = t->p_off + t->h_off;
		if(t->nid2 == mp->p2){
			if(t->p_off + t->h_off + n->len < mp->min_ins) return 0;
			else return 1;
		}
		for(i=0;i<count_u32v(n->links);i++){
			hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, i));
			if(hit->n_cigar == 0) continue;
			if(hit->rid2 == t->nid2) continue;
			nn = ref_nodev(mp->nodes, hit->rid2);
			if(nn->used) continue;
			if(n->off + (uint32_t)hit->off + nn->len > mp->max_ins) continue;
			push_trace_mp(mp, n->off, hit->off, t->nid2, hit->rid2, i);
		}
	}
	return 0;
}

uint32_t cal_pair_score(MyPath *mp){
	Node *n;
	uint32_t i, n_f, score, nid, len1, len2;
	score = 0;
	n_f = 0;
	len1 = ref_nodev(mp->nodes, 0)->len;
	len2 = ref_nodev(mp->nodes, 0)->len;
	for(i=1;i+1<count_u32v(mp->path);i++){
		nid = get_u32v(mp->path, i);
		n = ref_nodev(mp->nodes, nid);
		if(n->off < len1) continue;
		if(n->off + len2 >= mp->cns_len) continue;
		if((nid & 0x01)) score += n_f;
		else n_f ++;
	}
	return score;
}

void do_consensus_mypath(MyPath *mp){
	SR_AlnHit *hit;
	Node *n, *nn;
	cigarv *cs, *sc;
	Compose *cp;
	AlnCigar *ac;
	String *str;
	uint32_t i, j, c, off, np, nm, nid, tlen, cns_len;
	char seq[256];
	clear_u32v(mp->path);
	push_u32v(mp->path, mp->p2);
	n = ref_nodev(mp->nodes, mp->p2);
	n->inasm = 1;
	while(n->prev != 0xFFFFFFFFU){
		push_u32v(mp->path, n->prev);
		nn = ref_nodev(mp->nodes, n->prev);
		nn->idx = n->prev_idx;
		n = nn;
		n->inasm = 1;
	}
	reverse_u32v(mp->path);
	off = 0;
	hit = NULL;
	for(i=0;i<count_u32v(mp->path);i++){
		nid = get_u32v(mp->path, i);
		n = ref_nodev(mp->nodes, nid);
		if(count_cigarvv(mp->cigars) <= i){
			cs = next_ref_cigarvv(mp->cigars);
			cigarv_init(cs, 4);
		} else {
			cs = ref_cigarvv(mp->cigars, i);
			clear_cigarv(cs);
		}
		//char cstr[256];
		if(hit == NULL){
			encap_cigarv(cs, 3);
			cs->size = append_cigars(ref_cigarv(cs, 0), count_cigarv(cs), ALN_CIGAR_TYPE_MAT, sr_rdseq_length(mp->sdb, nid));
			hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, n->idx));
			n->off = 0;
		} else {
			sc = ref_cigarvv(mp->cigars, i - 1);
			nn = ref_nodev(mp->nodes, get_u32v(mp->path, i - 1));
			encap_cigarv(cs, hit->n_cigar + count_cigarv(sc));
			cs->size = cat_cigars(ref_cigarv(cs, 0), 0, hit->cigars, hit->n_cigar);
			//cigars2string(ref_cigarv(cs, 0), cs->size, cstr);
			//fprintf(stdout, "%08d %s\n", get_u32v(mp->path, i), cstr);
			flip_cigars(ref_cigarv(cs, 0), count_cigarv(cs));
			np = compile_cigars(mp->pat, hit->cigars, hit->n_cigar, ref_cigarv(sc, 0), count_cigarv(sc), 0);
			nm = apply_cigars(mp->mat, ref_cigarv(cs, 0), count_cigarv(cs), mp->pat, np);
			cs->size = cat_cigars(ref_cigarv(cs, 0), 0, mp->mat, nm);
			//cigars2string(ref_cigarv(cs, 0), cs->size, cstr);
			//fprintf(stdout, "%08d %s\n", get_u32v(mp->path, i), cstr);
			np = compile_cigars(mp->pat, ref_cigarv(sc, 0), count_cigarv(sc), hit->cigars, hit->n_cigar, 0);
			for(j=0;j<i;j++){
				sc = ref_cigarvv(mp->cigars, j);
				nm = apply_cigars(mp->mat, ref_cigarv(sc, 0), count_cigarv(sc), mp->pat, np);
				if(nm > count_cigarv(sc)) encap_cigarv(sc, nm - count_cigarv(sc));
				sc->size = cat_cigars(ref_cigarv(sc, 0), 0, mp->mat, nm);
				//char aln[256];
				//bits2seq(seq, ref_u64list(mp->sdb->rd_seqs, 0), sr_rdseq_offset(mp->sdb, get_u32v(mp->path, j)), sr_rdseq_length(mp->sdb, get_u32v(mp->path, j)));
				//cigars_seq2aln(aln, ref_cigarv(sc, 0), count_cigarv(sc), 0, seq);
				//fprintf(stdout, "%08d %s\n", get_u32v(mp->path, j), aln);
			}
			//fprintf(stdout, "\n");
			sc = ref_cigarvv(mp->cigars, i - 1);
			nn->off = select_cigars_seqlen(ref_cigarv(sc, 0), count_cigarv(sc), 1, 0) - 1;
			n->off  = select_cigars_seqlen(ref_cigarv(cs, 0), count_cigarv(cs), 1, 0) - 1;
			//fprintf(stdout, "%08d %d\n", get_u32v(mp->path, i-1), nn->off);
			//fprintf(stdout, "%08d %d\n", get_u32v(mp->path, i), n->off);
			//fprintf(stdout, "\n");
			if(i + 1 < count_u32v(mp->path)) hit = ref_sr_hitv(mp->sdb->hits, get_u32v(n->links, n->idx));
			else hit = NULL;
		}
	}
	cns_len = 0;
	for(i=0;i<count_u32v(mp->path);i++){
		cs = ref_cigarvv(mp->cigars, i);
		cigars_lengths(ref_cigarv(cs, 0), count_cigarv(cs), (int*)&tlen, NULL, NULL);
		if(tlen > cns_len) cns_len = tlen;
	}
	for(i=0;i<count_u32v(mp->path);i++){
		cs = ref_cigarvv(mp->cigars, i);
		cigars_lengths(ref_cigarv(cs, 0), count_cigarv(cs), (int*)&tlen, NULL, NULL);
		encap_cigarv(cs, 1);
		cs->size = append_cigars(ref_cigarv(cs, 0), count_cigarv(cs), ALN_CIGAR_TYPE_DEL, cns_len - tlen);
	}
	clear_u32list(mp->begs);
	clear_u32list(mp->ends);
	for(i=0;i<count_u32v(mp->path);i++){
		cs = ref_cigarvv(mp->cigars, i);
		for(j=0;j<count_cigarv(cs);j++){
			ac = ref_cigarv(cs, j);
			if(ac->type == ALN_CIGAR_TYPE_INS) ac->type = ALN_CIGAR_TYPE_MAT;
		}
		push_u32list(mp->begs, select_cigars_seqlen(ref_cigarv(cs, 0), count_cigarv(cs), 0, 0));
		push_u32list(mp->ends, select_cigars_seqlen(ref_cigarv(cs, 0), count_cigarv(cs), sr_rdseq_length(mp->sdb, get_u32v(mp->path, i)) - 1, 0));
		if(i >= count_strv(mp->matrix)){
			str = init_string(1024);
			push_strv(mp->matrix, str);
		} else {
			str = get_strv(mp->matrix, i);
			clear_string(str);
		}
		encap_string(str, cns_len);
		bits2seq(seq, ref_u64list(mp->sdb->rd_seqs, 0), sr_rdseq_offset(mp->sdb, get_u32v(mp->path, i)), sr_rdseq_length(mp->sdb, get_u32v(mp->path, i)));
		str->size = cigars_seq2aln(str->string, ref_cigarv(cs, 0), count_cigarv(cs), 0, seq);
		if(mp->output_cns) fprintf(stderr, "%08d %s\n", get_u32v(mp->path, i), str->string);
	}
	clear_string(mp->cns);
	clear_compv(mp->comps);
	for(j=0;j<cns_len;j++){
		cp = next_ref_compv(mp->comps);
		cp->acgtn[0] = cp->acgtn[1] = cp->acgtn[2] = cp->acgtn[3] = cp->acgtn[4] = 0;
		for(i=0;i<count_u32v(mp->path);i++){
			if(j < get_u32list(mp->begs, i)) continue;
			if(j > get_u32list(mp->ends, i)) continue;
			str = get_strv(mp->matrix, i);
			cp->acgtn[base_bit_table[(int)str->string[j]]] ++;
		}
		c = 0;
		for(i=1;i<4;i++){
			if(cp->acgtn[i] > cp->acgtn[c]){ c = i; }
		}
		add_char_string(mp->cns, bit_base_table[c]);
	}
	mp->cns_len = cns_len;
	mp->pair_score = cal_pair_score(mp);
}

void output_consensus_mypath(MyPath *mp, FILE *lr_seqs, FILE *lr_stcs){
	cigarv *cs;
	Node *n;
	Compose *cp;
	uint32_t i;
	fprintf(lr_seqs, ">longread%u len=%u n_rd=%d pair_score=%d\n", ref_nodev(mp->nodes, mp->p1)->rid >> 1, mp->cns_len, (unsigned)count_u32v(mp->path), mp->pair_score);
	for(i=0;i<mp->cns_len;i++){
		fputc(mp->cns->string[i], lr_seqs);
		if((i + 1) % 100 == 0) fputc('\n', lr_seqs);
	}
	if((i) % 100) fputc('\n', lr_seqs);
	fprintf(lr_stcs, "B %u %u %u %u\n", ref_nodev(mp->nodes, mp->p1)->rid >> 1, mp->cns_len, (unsigned)count_u32v(mp->path), mp->pair_score);
	fprintf(lr_stcs, "F pos\t-\tA\tC\tG\tT\n");
	for(i=0;i<count_compv(mp->comps);i++){
		cp = ref_compv(mp->comps, i);
		fprintf(lr_stcs, "F %u\t%u\t%u\t%u\t%u\t%u\n", i, cp->acgtn[4], cp->acgtn[0], cp->acgtn[1], cp->acgtn[2], cp->acgtn[3]);
	}
	for(i=0;i<count_u32v(mp->path);i++){
		cs = ref_cigarvv(mp->cigars, i);
		n = ref_nodev(mp->nodes, get_u32v(mp->path, i));
		cigars2string(ref_cigarv(cs, 0), count_cigarv(cs), mp->aux_str);
		fprintf(lr_stcs, "I %u\t%c\t%s\t", n->rid, "+-"[n->dir], mp->aux_str);
		bits2seq(mp->aux_str, ref_u64list(mp->sdb->rd_seqs, 0), sr_rdseq_offset(mp->sdb, get_u32v(mp->path, i)), sr_rdseq_length(mp->sdb, get_u32v(mp->path, i)));
		fprintf(lr_stcs, "%s\n", mp->aux_str);
	}
	for(i=0;i<count_nodev(mp->nodes);i++){
		n = ref_nodev(mp->nodes, i);
		if(n->inasm) continue;
		bits2seq(mp->aux_str, ref_u64list(mp->sdb->rd_seqs, 0), sr_rdseq_offset(mp->sdb, i), sr_rdseq_length(mp->sdb, i));
		fprintf(lr_stcs, "O %u\t%c\t-\t%s\n", n->rid, "+-"[n->dir], mp->aux_str);
	}
}

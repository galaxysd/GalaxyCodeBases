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
 
#include "local_assembly.h"

int cmp_trace(const void *e1, const void *e2, void *ref){
	Trace *t1, *t2;
	t1 = (Trace*)e1;
	t2 = (Trace*)e2;
	if(t1->n_ol == t2->n_ol) return 0;
	if(t1->n_ol <  t2->n_ol) return 1;
	else return -1;
	ref = ref; // useless
}

LGraph* init_lgraph(uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm, uint32_t gap_cutoff, uint32_t min_ins, uint32_t max_ins){
	LGraph *g;
	g = malloc(sizeof(LGraph));
	g->sdb         = sr_init_sdb(kmer_size, 4, rd_len);
	sr_set_align_parameters(g->sdb, 1, min_ol, min_sm, max_mm, 0);
	sr_set_filter_parameters(g->sdb, 0, 1024, 1);
	g->aux         = sr_init_aux();
	g->gap_cutoff  = gap_cutoff;
	g->nodes       = init_nodev(1024);
	g->m_node      = 0;
	g->trackers[0] = malloc(sizeof(Tracker));
	g->trackers[0]->trace_dir = 0;
	g->trackers[0]->traces    = init_tracev(1024);
	g->trackers[0]->m_trace   = 0;
	g->trackers[0]->heap      = init_heap(cmp_trace, g);
	g->trackers[1] = malloc(sizeof(Tracker));
	g->trackers[1]->trace_dir = 1;
	g->trackers[1]->traces    = init_tracev(1024);
	g->trackers[1]->m_trace   = 0;
	g->trackers[1]->heap      = init_heap(cmp_trace, g);
	g->p1          = 0;
	g->p2          = 1;
	g->lays        = init_layv(1024);
	g->cms         = init_cigarm(1024);
	g->srs         = init_strv(1024);
	g->cns_seq     = init_string(1024);
	g->ords        = init_u32list(1024);
	g->aux_ords    = init_u32list(1024);
	g->cns_index   = init_khash(1023);
	g->cns_id      = 0;
	g->cns_len     = 0;
	g->min_ins     = min_ins;
	g->max_ins     = max_ins;
	return g;
}

void free_lgraph(LGraph *g){
	uint32_t i;
	sr_free_sdb(g->sdb);
	sr_free_aux(g->aux);
	for(i=0;i<count_nodev(g->nodes);i++){
		free_edgev(ref_nodev(g->nodes, i)->edges[0]);
		free_edgev(ref_nodev(g->nodes, i)->edges[1]);
	}
	free_nodev(g->nodes);
	free_layv(g->lays);
	for(i=0;i<count_cigarm(g->cms);i++) free_cigarv(get_cigarm(g->cms, i));
	free_cigarm(g->cms);
	for(i=0;i<count_strv(g->srs);i++) free_string(get_strv(g->srs, i));
	free_strv(g->srs);
	free_string(g->cns_seq);
	free_u32list(g->ords);
	free_u32list(g->aux_ords);
	for(i=0;i<g->trackers[0]->m_trace;i++) free(get_tracev(g->trackers[0]->traces, i));
	free_tracev(g->trackers[0]->traces);
	for(i=0;i<g->trackers[1]->m_trace;i++) free(get_tracev(g->trackers[1]->traces, i));
	free_tracev(g->trackers[1]->traces);
	free_heap(g->trackers[0]->heap);
	free_heap(g->trackers[1]->heap);
	free_khash(g->cns_index);
	free(g->trackers[0]);
	free(g->trackers[1]);
	free(g);
}

void reset_lgraph(LGraph *g){
	sr_reset_sdb(g->sdb);
	clear_nodev(g->nodes);
	clear_sr_hitv(g->aux->hits);
	clear_layv(g->lays);
	clear_string(g->cns_seq);
	g->cns_len = 0;
	g->cns_id = 0;
}

void push_lgraph(LGraph *g, uint32_t rid, char *seq, uint8_t seqlen, uint32_t seqdir, uint32_t anchor_type, int min_ins, int max_ins){
	Node *node;
	node = next_ref_nodev(g->nodes);
	if(count_nodev(g->nodes) > g->m_node){
		node->edges[0] = init_edgev(64);
		node->edges[1] = init_edgev(64);
		g->m_node = count_nodev(g->nodes);
	} else {
		clear_edgev(node->edges[0]);
		clear_edgev(node->edges[1]);
	}
	node->rid     = rid;
	node->seqdir  = seqdir;
	node->len     = seqlen;
	node->anchor_type = anchor_type;
	node->min_ins = min_ins;
	node->max_ins = max_ins;
	node->bt      = 0;
	node->bt_idx  = 0;
	node->bt_off  = 0;
	node->bt_dir  = 0;
	node->visit   = 0;
	node->in_asm  = 0;
	node->asm_off = 0;
	sr_push_sdb(g->sdb, seq, seqlen);
}

void align_lgraph(LGraph *g){
	SR_AlnHit *hit;
	Node *n1, *n2;
	Edge *e;
	uint32_t i, offb, offe, rid1, rid2, n_ol;
	if(count_nodev(g->nodes) < 2) return;
	g->cns_id = ref_nodev(g->nodes, 0)->rid >> 1;
	sr_ready_sdb(g->sdb);
	//if(count_nodev(g->nodes) >= g->gap_cutoff){
		//g->sdb->allow_gap = 0;
	//} else {
		//g->sdb->allow_gap = 1;
	//}
	for(i=0;i<g->sdb->n_idx;i++) sr_index_sdb(g->sdb, i);
	clear_sr_hitv(g->aux->hits);
	for(i=0;i<g->sdb->n_rd;i++) sr_align_sdb(g->sdb, i, g->aux);
	for(i=0;i<count_sr_hitv(g->aux->hits);i++){
		hit = ref_sr_hitv(g->aux->hits, i);
		if(hit->dir1 ^ hit->dir2) continue; // nerver happen
		//sr_output_hit(g->sdb, hit, stdout);
		n_ol = hit->n_ol;
		if(hit->dir1){
			rid1 = hit->rid2;
			rid2 = hit->rid1;
			offe = hit->off;
			if(hit->off + ref_nodev(g->nodes, hit->rid2)->len < ref_nodev(g->nodes, hit->rid1)->len) continue;
			offb = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
		} else {
			rid1 = hit->rid1;
			rid2 = hit->rid2;
			offb = hit->off;
			if(hit->off + ref_nodev(g->nodes, hit->rid2)->len < ref_nodev(g->nodes, hit->rid1)->len) continue;
			offe = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
		}
		n1 = ref_nodev(g->nodes, rid1);
		n2 = ref_nodev(g->nodes, rid2);
		if(offb == offe && offb == 0){
			e = next_ref_edgev(n1->edges[0]); e->nid = rid2; e->n_ol = n_ol; e->off = 0; e->hit_idx = i;
			e = next_ref_edgev(n1->edges[1]); e->nid = rid2; e->n_ol = n_ol; e->off = 0; e->hit_idx = i;
			e = next_ref_edgev(n2->edges[0]); e->nid = rid1; e->n_ol = n_ol; e->off = 0; e->hit_idx = i;
			e = next_ref_edgev(n2->edges[1]); e->nid = rid1; e->n_ol = n_ol; e->off = 0; e->hit_idx = i;
		} else {
			e = next_ref_edgev(n1->edges[0]); e->nid = rid2; e->n_ol = n_ol; e->off = offb; e->hit_idx = i;
			e = next_ref_edgev(n2->edges[1]); e->nid = rid1; e->n_ol = n_ol; e->off = offe; e->hit_idx = i;
		}
	}
}

void push_tracker(Tracker *tracker, uint32_t src_nid, uint32_t hit_idx, uint32_t nid, uint32_t n_ol, uint32_t off){
	Trace *t;
	if(count_tracev(tracker->traces) >= tracker->m_trace){
		t = malloc(sizeof(Trace));
		push_tracev(tracker->traces, t);
		tracker->m_trace = count_tracev(tracker->traces);
	} else {
		t = get_tracev(tracker->traces, tracker->traces->size ++);
	}
	t->src_nid = src_nid;
	t->hit_idx = hit_idx;
	t->nid     = nid;
	t->n_ol    = n_ol;
	t->off     = off;
	push_heap(tracker->heap, t);
}

Trace* pop_tracker(Tracker *tracker){ return (Trace*)pop_heap(tracker->heap); }

int is_illegal_path(LGraph *g, uint32_t nid, int path_off, uint32_t path_dir){
	int a, b;
	Node *n;
	n = ref_nodev(g->nodes, nid);
	a = n->min_ins;
	b = n->max_ins;
	if((n->anchor_type & 0x01) == path_dir){
		if(path_off < a || path_off > b) return 0;
	} else {
		if(((int)g->max_ins) - (path_off + n->len) < a) return 0;
		if(((int)g->min_ins) - (path_off + n->len) > b) return 0;
	}
	return 1;
}

int layout_lgraph(LGraph *g){
	Node *n;
	Edge *e;
	Trace *t1, *t2;
	Layout *l;
	uint32_t i, bp1, bp2, hid;
	if(count_nodev(g->nodes) < 2) return 0;
	clear_heap(g->trackers[0]->heap);
	clear_tracev(g->trackers[0]->traces);
	clear_heap(g->trackers[1]->heap);
	clear_tracev(g->trackers[1]->traces);
	n = ref_nodev(g->nodes, g->p1);
	n->visit  = 0;
	n->bt     = g->p1;
	n->bt_off = 0;
	n->bt_dir = 0;
	push_tracker(g->trackers[0], g->p1, 0, g->p1, n->len, 0);
	n = ref_nodev(g->nodes, g->p2);
	n->visit  = 0;
	n->bt     = g->p2;
	n->bt_off = 0;
	n->bt_dir = 1;
	push_tracker(g->trackers[1], g->p2, 0, g->p2, n->len, 0);
	bp1 = bp2 = 0;
	hid = 0;
	while((t1 = pop_tracker(g->trackers[0])) && (t2 = pop_tracker(g->trackers[1]))){
		n = ref_nodev(g->nodes, t1->nid);
		if(n->visit == 0){
			n->visit = 1;
			n->bt     = t1->src_nid;
			n->bt_idx = t1->hit_idx;
			n->bt_off = t1->off;
			n->bt_dir = 0;
			for(i=0;i<count_edgev(n->edges[0]);i++){
				e = ref_edgev(n->edges[0], i);
				if(is_illegal_path(g, e->nid, e->off + n->bt_off, 0) == 0) continue;
				push_tracker(g->trackers[0], t1->nid, e->hit_idx, e->nid, e->n_ol, e->off + n->bt_off);
			}
		} else if(n->bt_dir == 1){
			if(n->bt_off + t1->off + n->len >= g->min_ins && n->bt_off + t1->off + n->len <= g->max_ins){
				bp1 = t1->src_nid;
				bp2 = t1->nid;
				hid = t1->hit_idx;
				break;
			}
		}
		n = ref_nodev(g->nodes, t2->nid);
		if(n->visit == 0){
			n->visit = 1;
			n->bt     = t2->src_nid;
			n->bt_idx = t2->hit_idx;
			n->bt_off = t2->off;
			n->bt_dir = 1;
			for(i=0;i<count_edgev(n->edges[1]);i++){
				e = ref_edgev(n->edges[1], i);
				if(is_illegal_path(g, e->nid, e->off + n->bt_off, 1) == 0) continue;
				push_tracker(g->trackers[1], t2->nid, e->hit_idx, e->nid, e->n_ol, e->off + n->bt_off);
			}
		} else if(n->bt_dir == 0){
			if(n->bt_off + t2->off + n->len >= g->min_ins && n->bt_off + t2->off + n->len <= g->max_ins){
				bp1 = t2->nid;
				bp2 = t2->src_nid;
				hid = t2->hit_idx;
				break;
			}
		}
	}
	if(bp1 == bp2) return 0;
	clear_layv(g->lays);
	while(1){
		n = ref_nodev(g->nodes, bp1);
		l = next_ref_layv(g->lays);
		l->nid = bp1;
		l->hit_idx = n->bt_idx;
		if(n->bt == bp1) break;
		bp1 = n->bt;
	}
	reverse_layv(g->lays);
	while(1){
		n = ref_nodev(g->nodes, bp2);
		l = next_ref_layv(g->lays);
		l->nid = bp2;
		l->hit_idx = hid;
		if(n->bt == bp2) break;
		bp2 = n->bt;
		hid = n->bt_idx;
	}
	return 1;
}

void consensus_lgraph(LGraph *g){
	uint32_t i, j, dir, np, nm, tlen, acgtn[5], c;
	Node *n;
	Layout *lay;
	cigarv *cs, *sc;
	AlnCigar *ac, pat[64], mat[64], ori[64];
	SR_AlnHit *hit;
	String *str;
	char seq[256];
	for(i=count_cigarm(g->cms);i<count_nodev(g->nodes);i++) push_cigarm(g->cms, init_cigarv(64));
	for(i=count_strv(g->srs);i<count_nodev(g->nodes);i++) push_strv(g->srs, init_string(1024));
	clear_u32list(g->ords);
	for(i=0;i<count_layv(g->lays);i++){
		lay = ref_layv(g->lays, i);
		push_u32list(g->ords, lay->nid);
		n   = ref_nodev(g->nodes, lay->nid);
		n->in_asm = 1;
		cs  = get_cigarm(g->cms, lay->nid);
		clear_cigarv(cs);
		if(i == 0){
			encap_cigarv(cs, 1);
			cs->size = append_cigars(ref_cigarv(cs, 0), 0, ALN_CIGAR_TYPE_MAT, n->len);
			continue;
		}
		hit = ref_sr_hitv(g->aux->hits, lay->hit_idx);
		dir = hit->dir1 ^ (hit->rid1 == lay->nid);
		if(dir){
			cat_cigars(ori, 0, hit->cigars, hit->n_cigar);
			reverse_cigars(ori, hit->n_cigar);
			flip_cigars(ori, hit->n_cigar);
		} else {
			cat_cigars(ori, 0, hit->cigars, hit->n_cigar);
		}
		//sr_output_hit(g->sdb, hit, stdout);
		sc = get_cigarm(g->cms, ref_layv(g->lays, i - 1)->nid);
		//cigars2string(ref_cigarv(sc, 0), count_cigarv(sc), seq); fprintf(stdout, "%s\n", seq);
		encap_cigarv(cs, hit->n_cigar + count_cigarv(sc));
		cs->size = cat_cigars(ref_cigarv(cs, 0), 0, ori, hit->n_cigar);
		flip_cigars(ref_cigarv(cs, 0), count_cigarv(cs));
		np = compile_cigars(pat, ori, hit->n_cigar, ref_cigarv(sc, 0), count_cigarv(sc), 0);
		nm = apply_cigars(mat, ref_cigarv(cs, 0), count_cigarv(cs), pat, np);
		cs->size = cat_cigars(ref_cigarv(cs, 0), 0, mat, nm);
		//cigars2string(ref_cigarv(cs, 0), count_cigarv(cs), seq); fprintf(stdout, "%010u\t%s\n", ref_nodev(g->nodes, lay->nid)->rid, seq);
		np = compile_cigars(pat, ref_cigarv(sc, 0), count_cigarv(sc), ori, hit->n_cigar, 0);
		for(j=0;j<i;j++){
			sc = get_cigarm(g->cms, ref_layv(g->lays, j)->nid);
			nm = apply_cigars(mat, ref_cigarv(sc, 0), count_cigarv(sc), pat, np);
			encap_cigarv(sc, nm);
			sc->size = cat_cigars(ref_cigarv(sc, 0), 0, mat, nm);
		}
	}
	g->cns_len = 0;
	for(i=count_layv(g->lays)-1;i&&i+10>count_layv(g->lays);i--){
		cs = get_cigarm(g->cms, ref_layv(g->lays, i)->nid);
		cigars_lengths(ref_cigarv(cs, 0), count_cigarv(cs), (int*)&tlen, NULL, NULL);
		if(tlen > g->cns_len) g->cns_len = tlen;
	}
	for(i=0;i<count_layv(g->lays);i++){
		cs = get_cigarm(g->cms, ref_layv(g->lays, i)->nid);
		cs->size = refine_cigars(ref_cigarv(cs, 0), count_cigarv(cs));
		for(j=0;j<count_cigarv(cs);j++){
			ac = ref_cigarv(cs, j);
			if(ac->type == ALN_CIGAR_TYPE_INS) ac->type = ALN_CIGAR_TYPE_MAT;
		}
		cigars_lengths(ref_cigarv(cs, 0), count_cigarv(cs), (int*)&tlen, NULL, NULL);
		if(tlen < g->cns_len){
			encap_cigarv(cs, 1);
			cs->size = append_cigars(ref_cigarv(cs, 0), count_cigarv(cs), ALN_CIGAR_TYPE_DEL, g->cns_len - tlen);
		}
		ref_nodev(g->nodes, ref_layv(g->lays, i)->nid)->bt_off = select_cigars_seqlen(ref_cigarv(cs, 0), count_cigarv(cs), 0, 0);
		str = get_strv(g->srs, ref_layv(g->lays, i)->nid);
		clear_string(str);
		encap_string(str, g->cns_len);
		bits2seq(seq, ref_u64list(g->sdb->rd_seqs, 0), 
				sr_rdseq_offset(g->sdb, ref_layv(g->lays, i)->nid),
				sr_rdseq_length(g->sdb, ref_layv(g->lays, i)->nid));
		str->size = cigars_seq2aln(str->string, ref_cigarv(cs, 0), count_cigarv(cs), 0, seq);
		//fprintf(stdout, "%010u\t%s\n", ref_layv(g->lays, i)->nid, str->string);
	}
	clear_string(g->cns_seq);
	for(i=0;i<g->cns_len;i++){
		acgtn[0] = 0;
		acgtn[1] = 0;
		acgtn[2] = 0;
		acgtn[3] = 0;
		acgtn[4] = 0;
		for(j=0;j<count_layv(g->lays);j++){
			str = get_strv(g->srs, ref_layv(g->lays, j)->nid);
			acgtn[base_bit_table[(int)str->string[i]]] ++;
		}
		acgtn[4] = 0;
		c = 4;
		for(j=0;j<4;j++){
			if(acgtn[j] > acgtn[c]) c = j;
		}
		add_char_string(g->cns_seq, bit_base_table[c]);
	}
}

int aln_rd2cns_dp(LGraph *g, char *rdseq, uint32_t rdlen, cigarv *cs){
	AlnAln *aln;
	path_t *p;
	AlnCigar *ac;
	uint32_t i, mm;
	int ret;
	aln = aln_stdaln_aux(rdseq, g->cns_seq->string, &g->sdb->ap, 1, 0, rdlen, g->cns_seq->size);
	clear_cigarv(cs);
	p = aln->path + aln->path_len - 1;
	ac = NULL;
	mm = 0;
	while(p >= aln->path){
		if(p->ctype == FROM_M && rdseq[p->i] != g->cns_seq->string[p->j]) mm ++;
		if(ac && p->ctype == ac->type) ac->len ++;
		else {
			ac = next_ref_cigarv(cs);
			ac->type = p->ctype;
			ac->len = 1;
		}
		p --;
	}
	for(i=0;i<count_cigarv(cs);i++){
		ac = ref_cigarv(cs, i);
		switch(ac->type){
			case FROM_M: ac->type = ALN_CIGAR_TYPE_MAT; break;
			case FROM_I: ac->type = ALN_CIGAR_TYPE_DEL; break;
			case FROM_D: ac->type = ALN_CIGAR_TYPE_CLIP1; break;
		}
	}
	ret = 1;
	// Locate inside cns sequence
	if(count_cigarv(cs)){
		ac = ref_cigarv(cs, 0);
		if(ac->type == ALN_CIGAR_TYPE_CLIP1) ret = 0;
		ac = ref_cigarv(cs, count_cigarv(cs) - 1);
		if(ac->type == ALN_CIGAR_TYPE_CLIP1) ret = 0;
	} else ret = 0;
	// Mismatch
	if(mm > rdlen * (1 - g->sdb->min_similarity)) ret = 0;
	aln_free_AlnAln(aln);
	return ret;
}

int aln_rd2cns(LGraph *g, char *rdseq, int rdlen, cigarv *cs, i32hash *hash){
	AlnCigar *ac;
	int i, mm, min, x, y, pos, p;
	min = rdlen * (1 - g->sdb->min_similarity) + 1;
	p = 10000;
	pos = 0;
	while(iter_i32hash(hash, &pos)){
		if(pos < 0) x = -pos;
		else x = 0;
		if(pos + rdlen > (int)g->cns_len) y = g->cns_len - pos;
		else y = rdlen;
		mm = 0;
		for(i=x;i<y;i++){
			if(g->cns_seq->string[i+pos] != rdseq[i]) mm ++;
		}
		if(mm < min){ p = pos; min = mm; }
	}
	if(min < (int)(rdlen * (1 - g->sdb->min_similarity) + 1)){
		clear_cigarv(cs);
		if(p < 0) x = -p;
		else x = 0;
		if(p + rdlen > (int)g->cns_len) y = g->cns_len - p;
		else y = rdlen;
		if(x){
			ac = next_ref_cigarv(cs);
			ac->type = ALN_CIGAR_TYPE_CLIP1;
			ac->len  = x;
		} else if(p > 0){
			ac = next_ref_cigarv(cs);
			ac->type = ALN_CIGAR_TYPE_DEL;
			ac->len  = p;
		}
		ac = next_ref_cigarv(cs);
		ac->type = ALN_CIGAR_TYPE_MAT;
		ac->len  = y - x;
		if(y < rdlen){
			ac = next_ref_cigarv(cs);
			ac->type = ALN_CIGAR_TYPE_CLIP1;
			ac->len  = rdlen - y;
		} else if(p + rdlen < (int)g->cns_len){
			ac = next_ref_cigarv(cs);
			ac->type = ALN_CIGAR_TYPE_DEL;
			ac->len  = ((int)g->cns_len) - (p + rdlen);
		}
		return 1;
	} else if(count_i32hash(hash) >= 5){
		return aln_rd2cns_dp(g, rdseq, rdlen, cs);
	} else return 0;
}

int get_cigarv_offset(uint32_t nid, LGraph *g){
	AlnCigar *ac;
	ac = ref_cigarv(get_cigarm(g->cms, nid), 0);
	switch(ac->type){
		case ALN_CIGAR_TYPE_CLIP1: return - ((int)ac->len);
		case ALN_CIGAR_TYPE_DEL  : return ac->len;
		default: return 0;
	}
}

int cmp_cigarv_off_func(uint32_t nid1, uint32_t nid2, void *obj){
	return get_cigarv_offset(nid1, (LGraph*)obj) - get_cigarv_offset(nid2, (LGraph*)obj);
}

define_quick_sort(qsort_cg_ords, uint32_t, cmp_cigarv_off_func);

void realign_lgraph(LGraph *g){
	cigarv *cs;
	String *str;
	Kmer *kmer, KMER;
	i32hash *hash;
	char rdseq[256];
	uint32_t i, j, k, x, ksize, kmask, nid, nid1, nid2, rdlen, acgtn[5], c, mm;
	int flag, offset;
	ksize = 7;
	kmask = (1LLU << (ksize << 1)) - 1;
	k = 0;
	KMER.kmer = 0;
	KMER.pos  = 0;
	clear_khash(g->cns_index);
	for(i=0;i+ksize<=g->cns_len;i++){
		k = (k << 2) | (base_bit_table[(int)g->cns_seq->string[i]]);
		if(i + 1 < ksize) continue;
		k = k & kmask;
		KMER.kmer = k;
		KMER.pos  = i + 1 - ksize;
		put_khash(g->cns_index, KMER);
	}
	clear_u32list(g->aux_ords);
	hash = init_i32hash(1023);
	for(nid=0;nid<count_nodev(g->nodes);nid+=2){
		flag = 0;
		for(j=0;j<2;j++){
			if(ref_nodev(g->nodes, nid + j)->in_asm){ flag |= (1 << j);  continue; }
			cs = get_cigarm(g->cms, nid + j);
			clear_cigarv(cs);
			rdlen = sr_rdseq_length(g->sdb, nid + j);
			bits2seq(rdseq, ref_u64list(g->sdb->rd_seqs, 0), sr_rdseq_offset(g->sdb, nid + j), rdlen);
			clear_i32hash(hash);
			k = 0;
			for(i=0;i+ksize<=rdlen;i++){
				k = (k << 2) | (base_bit_table[(int)rdseq[i]]);
				if(i + 1 < ksize) continue;
				k = k & kmask;
				KMER.kmer = k;
				kmer = get_khash(g->cns_index, KMER);
				if(kmer == NULL) continue;
				//if(kmer->pos < i + 1 - ksize) continue;
				//if(kmer->pos + (rdlen - (i + 1 - ksize)) > g->cns_len) continue;
				put_i32hash(hash, ((int)kmer->pos) - (i + 1 - ksize));
			}
			if(count_i32hash(hash) == 0) continue;
			if(aln_rd2cns(g, rdseq, rdlen, cs, hash)) flag |= (1 << j);
		}
		if(flag == 3){
			push_u32list(g->aux_ords, nid + 0);
			push_u32list(g->aux_ords, nid + 1);
		}
	}
	for(i=0;i<count_u32list(g->aux_ords);i++){
		nid = get_u32list(g->aux_ords, i);
		ref_nodev(g->nodes, nid)->asm_off = get_cigarv_offset(nid, g);
	}
	qsort_cg_ords(ref_u32list(g->aux_ords, 0), count_u32list(g->aux_ords), g);
	x = 0;
	clear_i32hash(hash);
	for(i=1;i<=count_u32list(g->aux_ords);i++){
		if(i == count_u32list(g->aux_ords) || ref_nodev(g->nodes, get_u32list(g->aux_ords, i))->asm_off != ref_nodev(g->nodes, get_u32list(g->aux_ords, x))->asm_off){
			if(i - x > 1){
				for(j=x;j<i;j++){
					nid1 = get_u32list(g->aux_ords, j);
					for(k=j+1;k<i;k++){
						nid2 = get_u32list(g->aux_ords, k);
						if(ref_nodev(g->nodes, nid1 ^ 0x1U)->asm_off != ref_nodev(g->nodes, nid2 ^ 0x1U)->asm_off) continue;
						if(strcmp(get_strv(g->srs, nid1)->string, get_strv(g->srs, nid2)->string) != 0) continue;
						if(strcmp(get_strv(g->srs, nid1 ^ 0x1U)->string, get_strv(g->srs, nid2 ^ 0x1U)->string) != 0) continue;
						put_i32hash(hash, nid2);
						put_i32hash(hash, nid2 ^ 0x1U);
					}
				}
			}
			x = i;
		}
	}
	clear_u32list(g->ords);
	for(i=0;i<count_u32list(g->aux_ords);i++){
		nid = get_u32list(g->aux_ords, i);
		if(exists_i32hash(hash, nid)) continue;
		push_u32list(g->ords, nid);
	}
	free_i32hash(hash);
	for(i=0;i<count_u32list(g->ords);i++){
		nid = get_u32list(g->ords, i);
		cs  = get_cigarm(g->cms, nid);
		if(ref_nodev(g->nodes, nid)->in_asm) continue;
		ref_nodev(g->nodes, nid)->in_asm = 1;
		str = get_strv(g->srs, nid);
		clear_string(str);
		encap_string(str, g->cns_len);
		bits2seq(rdseq, ref_u64list(g->sdb->rd_seqs, 0), sr_rdseq_offset(g->sdb, nid), sr_rdseq_length(g->sdb, nid));
		str->size = cigars_seq2aln(str->string, ref_cigarv(cs, 0), count_cigarv(cs), 0, rdseq);
	}
	clear_string(g->cns_seq);
	for(i=0;i<g->cns_len;i++){
		acgtn[0] = 0;
		acgtn[1] = 0;
		acgtn[2] = 0;
		acgtn[3] = 0;
		acgtn[4] = 0;
		for(j=0;j<count_u32list(g->ords);j++){
			str = get_strv(g->srs, get_u32list(g->ords, j));
			acgtn[base_bit_table[(int)str->string[i]]] ++;
		}
		acgtn[4] = 0;
		c = 4;
		for(j=0;j<4;j++){
			if(acgtn[j] > acgtn[c]) c = j;
		}
		add_char_string(g->cns_seq, bit_base_table[c]);
	}
	for(i=0;i<count_u32list(g->ords);i++){
		str = get_strv(g->srs, get_u32list(g->ords, i));
		mm = 0;
		for(j=0;j<g->cns_len;j++){
			if(str->string[j] != '-' && str->string[j] != g->cns_seq->string[j]){
				str->string[j] += 'a' - 'A'; mm ++;
			}
		}
		ref_nodev(g->nodes, get_u32list(g->ords, i))->asm_mm = mm;
		cs = get_cigarm(g->cms, get_u32list(g->ords, i));
		offset = 0;
		for(j=0;j<count_cigarv(cs);j++){
			if(ref_cigarv(cs, j)->type == ALN_CIGAR_TYPE_CLIP1){
				offset = ((int)ref_cigarv(cs, j)->len) - ((int)ref_nodev(g->nodes, get_u32list(g->ords, i))->len);
				break;
			} else if(ref_cigarv(cs, j)->type == ALN_CIGAR_TYPE_DEL){
				offset += ref_cigarv(cs, j)->len;
			} else break;
		}
		ref_nodev(g->nodes, get_u32list(g->ords, i))->asm_off = offset;
	}
}

void output_lgraph(LGraph *g, FILE *cnsf, FILE *samf, FILE *msaf){
	uint32_t i, nid, b, e, f;
	Node *n, *nn;
	cigarv *cs;
	char ac[128], seq[256];
	if(cnsf){
		fprintf(cnsf, ">fis%u len=%u n_rd=%u pair_score=%u\n", g->cns_id, g->cns_len, (unsigned)count_u32list(g->ords), 0);
		for(i=0;i<g->cns_len;i++){
			fputc(g->cns_seq->string[i], cnsf);
			if(((i + 1) % 100) == 0) fputc('\n', cnsf);
		}
		if(i % 100) fputc('\n', cnsf);
	}
	if(samf){
		for(i=0;i<count_u32list(g->ords);i++){
			nid = get_u32list(g->ords, i);
			n = ref_nodev(g->nodes, nid);
			if(nid & 0x01){
				f = 0x80;
				nn = ref_nodev(g->nodes, nid - 1);
			} else {
				f = 0x40;
				nn = ref_nodev(g->nodes, nid + 1);
			}
			if(n->seqdir) f |= 0x10;
			if(nn->in_asm){
				if(nn->seqdir) f |= 0x20;
			} else {
				f |= 0x08;
			}
			cs = get_cigarm(g->cms, nid);
			b = 0;
			while(ref_cigarv(cs, b)->type == ALN_CIGAR_TYPE_DEL) b ++;
			e = count_cigarv(cs) - 1;
			while(ref_cigarv(cs, e)->type == ALN_CIGAR_TYPE_DEL) e --;
			cigars2string(ref_cigarv(cs, b), e - b + 1, ac);
			bits2seq(seq, ref_u64list(g->sdb->rd_seqs, 0), sr_rdseq_offset(g->sdb, nid), sr_rdseq_length(g->sdb, nid));
			fprintf(samf, "%u\t%u\tfis%u\t%u\t%u\t%s\t%c\t%u\t%d\t%s\t*\tXT:A:U\tNM:i;%u\n",
					n->rid,
					f,
					g->cns_id,
					n->bt_off,
					255,
					ac,
					nn->in_asm? '=' : '*',
					nn->in_asm? nn->bt_off : 0,
					0,
					seq,
					n->asm_mm
				   );
		}
	}
	if(msaf){
		fprintf(msaf, ">fis%u len=%u n_rd=%u pair_score=%u\n", g->cns_id, g->cns_len, (unsigned)count_u32list(g->ords), 0);
		for(i=0;i<count_u32list(g->ords);i++){
			cs = get_cigarm(g->cms, get_u32list(g->ords, i));
			cigars2string(ref_cigarv(cs, 0), count_cigarv(cs), ac);
			n = ref_nodev(g->nodes, get_u32list(g->ords, i));
			fprintf(msaf, "%010u %s\t%s\t%d\t%u:%d:%d\n", ref_nodev(g->nodes, get_u32list(g->ords, i))->rid, get_strv(g->srs, get_u32list(g->ords, i))->string, ac, n->asm_off, n->anchor_type, n->min_ins, n->max_ins);
		}
	}
}


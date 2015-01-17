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
 

static inline rt_next_node(ReadTracker *rt){
	Node *n;
	if(rt->n_node < count_nlist(rt->nodes)){
		n = ref_nlist(rt->nodes, rt->n_node);
		clear_linkv(n->links);
	} else {
		n = next_ref_nlist(rt->nodes);
		n->links = init_linkv(4);
	}
	rt->n_node ++;
	n->k = NULL;
	n->kmer = 0;
	n->dir = 0;
	n->closed = 0;
	n->x = 0;
	return n;
}

uint32_t clone_traces(ReadTracker *rt, uint32_t dst, uint32_t src){
	//TODO clone
	return 0;
}

uint32_t link_traces(ReadTracker *rt, uint32_t tid){
	//TODO link
	return 0;
}

uint32_t graph_search_readpath_core(ReadTracker *rt, uint32_t seed_tid, uint32_t tdir, uint32_t roff){
	Kmer *k, K;
	Nod NOD, *nod;
	Node *n;
	Link *link;
	Trace *t, *tt;
	uint32_t i, tid, ttid;
	clear_u32list(rt->ends[tdir]);
	clear_u32list(rt->stack);
	push(rt->stack, seed_tid);
	while(pop_u32list(rt->stack, &tid)){
		t = ref_tlist(rt->traces, tid);
		n = ref_nlist(rt->nodes, t->nid);
		if(tdir? (n->x == 0) : (n->x == rt->rlen - 1)){
			push_u32list(rt->ends[tdir], tid);
			link_traces(rt, tid);
			continue;
		}
		if(n->closed){
			for(i=0;i<count_linkv(n->links);i++){
				link = ref_linkv(n->links, i);
				ttid = clone_traces(rt, tid, link->tid);
				if(ref_tlist(rt->traces, ttid)->score <= rt->max_score){
					push_u32list(rt->stack, ttids[1]);
				} else {
					push_u32list(rt->ends[tdir], tid);
				}
			}
			continue;
		}
	}
}

uint32_t graph_search_readpath(ReadTracker *rt, uint32_t rid, uint32_t rlen, char *seqs){
	Kmer *k, K;
	Node *n;
	Nod  NOD, *nod;
	Trace *t, *tt;
	uint64_t kk, v;
	uint32_t i, j, tdir, x;
	if(rlen < rt->kmer_size) return 0;
	rt->rid  = rid;
	rt->rlen = rlen;
	rt->seqs = seqs;
	memset(&K, 0, sizeof(Kmer));
	memset(&NOD, 0, sizeof(Nod));
	for(rt->rdir=0;rt->rdir<2;rt->rdir++){
		clear_u8list(rt->seqs);
		for(i=0;i<rt->rlen;i++) push_u8list(rt->seqs, rt->rdir? ((~base_bit_table[(int)seqs[rlen - i - 1]]) & 0x03) : (base_bit_table[(int)seqs[i]]));
		rt->n_node = 0;
		clear_nhash(rt->hash);
		rt->kmer = 0;
		for(i=0;i+1<rt->kmer-size;i++) rt->kmer = (rt->kmer << 2) | get_u8list(rt->seqs, i);
		for(rt->roff=0;rt->roff+rt->kmer_size<=rt->rlen;rt->roff++){
			NOD.x = rt->roff;
			rt->kmer = ((rt->kmer << 2) | get_u8list(rt->seqs, rt->roff + rt->kmer_size)) & rt->kmer_mask;
			for(tdir=0;tdir<2;tdir++){
				NOD.dir = tdir;
				NOD.kmer = rt->kmer;
				nod = get_nhash(rt->hash, NOD);
				if(nod == NULL){
					K.kmer = rt->kmer;
					k = get_khash(rt->g, K);
					if(k == NULL) continue;
					NOD.nid = rt->n_node;
					nod = put_nhash(rt->hash, NOD);
					n = rt_next_node(rt);
					n->kmer = rt->kmer;
					n->k    = k;
					n->x    = rt->roff;
					n->dir  = rt->rdir;
					n->closed = 0;
				}
				n = ref_nlist(rt->nodes, nod->nid);
			}
		}
	}
}

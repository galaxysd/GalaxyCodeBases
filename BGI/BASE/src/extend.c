/*

   extend.c        extend in one direction of BASE.

#    Copyright (C) 2015, The University of Hong Kong.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 3
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  03110-1301, USA.

    Date   : 28th Jan 2016
    Author : Binghang Liu, bhliu@cs.hku.hk
*/

#include "extend.h"

/*
    Get reads id and keep them in reads.
*/

int set_read_ids_plus(ULL l, ULL r, int* tCount, int level, int node, int size, int cons, int threadId)
{
    ULL num = r-l+1;
    int c = *tCount;
    ULL i;
    if(num > MAXDEPTH)
    {
        fprintf(stderr, "Index PLus Error: left: %d right %d\n", l, r);
        return 0;
    }
    for(i=0; i< num; ++i)
    {
        ULL id = idx2BWT[size]->readIDtable[l + i - idx2BWT[size]->bwt->cumulativeFreq[4]];
        g_reads[threadId][c+i].id = id;
        g_reads[threadId][c+i].strand = 1;
        g_reads[threadId][c+i].size = size;
        g_reads[threadId][c+i].level = level;
        g_reads[threadId][c+i].node = node;
        g_reads[threadId][c+i].cons = cons;
    }
    *tCount = c + num;
    return 1;
}

int set_read_ids_minus(ULL l, ULL r, int* tCount, int level, int node, int size, int cons, int threadId)
{
    ExactRange er;
    er.saL = l; er.saR = r; er.strand = 0;
    if(r-l+1 > MAXDEPTH)
    {
        fprintf(stderr, "Index Minus Error: left: %d right %d\n", l, r);
        return 0;
    }
    ReadInf ri[r-l+2];
    int add_num = extractReadInf(idx2BWT[size], er, ri, r-l+2);
    int i, b = (*tCount);

    for( i=0; i< add_num; i++)
    {
        g_reads[threadId][b+i].id = ri[i].read_id;
        g_reads[threadId][b+i].level = level;
        g_reads[threadId][b+i].strand = 0;
        g_reads[threadId][b+i].node = node;
        g_reads[threadId][b+i].size = size;
        g_reads[threadId][b+i].cons = cons;
    }
    (*tCount) = b + add_num;
    return 1;
}

/*
    Add new bases, including ACGT$.
*/
int add_base(Base* b, int base, Base *s)
{
    ULL l, r, ll, rr, ll_rev, rr_rev, depth=0;
    int i;

    for(i=0; i< bwtCount; i++)
    {
        s[i].plus = 0; s[i].minus = 0;
        s[i].saL = 0; s[i].saR = 0; s[i].saLL = 0; s[i].saRR = 0; s[i].saLLrev = 0; s[i].saRRrev = 0;
        if(b[i].plus > 0 && b[i].saR >= b[i].saL && b[i].saR - b[i].saL < MAXDEPTH)
        {
            l = b[i].saL; r = b[i].saR;
            BWTSARangeBackwardFor(idx2BWT[i], base ,&l, &r);
            if(l <= r && r-l <= b[i].saR - b[i].saL)
            {
                s[i].saL = l;
                s[i].saR = r;
                s[i].plus = r - l + 1;
            }
        }
        if(b[i].minus > 0 && b[i].saRR >= b[i].saLL && b[i].saRR - b[i].saLL < MAXDEPTH)
        {
            ll = b[i].saLL; rr = b[i].saRR; ll_rev = b[i].saLLrev; rr_rev = b[i].saRRrev;
            BWTSARangeBackwardRev(idx2BWT[i], 3-base ,&ll, &rr, &ll_rev, &rr_rev);
            if(ll <= rr && rr-ll <= b[i].saRR - b[i].saLL)
            {
                s[i].saLL = ll;
                s[i].saRR = rr;
                s[i].saLLrev = ll_rev;
                s[i].saRRrev = rr_rev;
                s[i].minus = rr -ll + 1;
            }
        }
        depth += s[i].plus + s[i].minus;
    }

    return depth;
}

int add_end(Base* b, int* read_num, int level, int node, int cons, int threadId)
{
    ULL l, r, ll, rr, ll_rev, rr_rev, total=0;
    int i, flag;

    for(i=0; i< bwtCount; i++)
    {
        if(b[i].plus > 0 && b[i].saR >= b[i].saL && b[i].saR - b[i].saL < MAXDEPTH)
        {
            l = b[i].saL; r = b[i].saR;
            BWTSARangeBackward(idx2BWT[i], 4 ,&l, &r);
            if(l <= r && r-l <= b[i].saR - b[i].saL)
            {
                flag = set_read_ids_plus(l, r, read_num, level, node, i, cons, threadId);
                if(flag == 1)
                    total += r - l + 1;
            }
        }
        if(b[i].minus > 0 && b[i].saRR >= b[i].saLL && b[i].saRR - b[i].saLL < MAXDEPTH)
        {
            ll = b[i].saLL; rr = b[i].saRR; ll_rev = b[i].saLLrev; rr_rev = b[i].saRRrev;
            BWTSARangeBackwardRev(idx2BWT[i], 4 ,&ll, &rr, &ll_rev, &rr_rev);
            if(ll <= rr && rr-ll <= b[i].saRR - b[i].saLL)
            {
                flag = set_read_ids_minus(ll, rr, read_num, level, node, i, cons, threadId);
                if(flag == 1)
                    total += rr - ll + 1;
            }
        }
    }
    return total;
}

int get_real_depth(Base* b)
{
    int i;
    ULL depth = 0;
    for(i=0; i< bwtCount; i++)
        depth += b[i].plus + b[i].minus;
    return depth; 
}

int get_strand_depth(Base* b, int strand)
{
    int i;
    ULL depth = 0;
    for(i=0; i< bwtCount; i++)
    {
        if(strand == 1)
            depth += b[i].plus;
        else
            depth += b[i].minus;
    }
    return depth;
}
void set_blank_base(Base* to)
{
    int i;
    for(i=0; i< bwtCount; i++)
    {
        to[i].plus = 0;
        to[i].minus = 0;
    }
}

void print_inf(FILE* fp, Node* nodes, Lev* levs, ReadId* reads, int* s, int cur_len);
void print_base(FILE* fp, Base* b)
{
    int i;
    for(i=0; i< bwtCount; i++)
    {
        fprintf(fp, "i: %d l %lld r %lld ll %lld rr %lld plus %d minus %d\n", i, b->saL, b->saR, b->saLL, b->saRR, b->plus, b->minus); 
    }
}
/*
    construct the backward extension tree.
*/
int backward_search(FILE* logs, int* ss, int threadId)
{
    int num[3]={ss[0], ss[1], ss[2]};
    int i,j, k;
    ULL l, r, ll, rr, ll_rev, rr_rev;

    //define left to speed up.
    unsigned long long left, flag;

    g_levs[threadId][num[1]].nnum=0;

    int max_i=0, max_j=0, max_depth=0, count, cons, node_max_i = -1, node_max_j = 0, node_max_depth =0;
    //check ACGT.
    for(j=num[0]-1; j>=0; j--)
    {
        if(g_nodes[threadId][j].is_end == 1)
            continue;
        if(g_nodes[threadId][j].level != num[1]-1)
            continue;

        left = get_real_depth(g_nodes[threadId][j].b);
        if(left == 0)
        {
            fprintf(stderr, "j=%d, b[0].l=%d, b[0].r=%d, b[0].ll=%d, b[0].rr=%d\n", j, g_nodes[threadId][j].b[0].saL, g_nodes[threadId][j].b[0].saR, g_nodes[threadId][j].b[0].saLL, g_nodes[threadId][j].b[0].saRR);
            exit(7);
        }
    	cons = 0;
    	if(j == 0 && num[1] == 1)
            cons = 1;
        if(j == g_levs[threadId][num[1]-1].node && g_levs[threadId][ num[1] - g_nodes[threadId][j].len ].node == j)
            cons = 1;

        flag = add_end(g_nodes[threadId][j].b, &(num[2]), num[1], j, cons, threadId);
        left -= flag;
        g_levs[threadId][num[1]].c[4] += flag;
        
        if(left <= 0)
        {
            g_nodes[threadId][j].is_end = 1;
            if(max_depth == 0 || max_depth < flag)
            {
                max_j = j;
                max_i = 4;
                max_depth = flag;
            }
            continue;
        }

        count=0;
        int now_max_depth=0, now_max_i=0, now_noend_depth = left;
        for(i=0; i<4; i++)
        {
            if(left == 0)
                break;
            flag = add_base(g_nodes[threadId][j].b, i, sbase[threadId][i]);
            if(flag == 0)
                continue;

            left -= flag;
            g_levs[threadId][num[1]].c[i] += flag;
            count++;
            if(now_max_depth < flag)
            {
                now_max_i = i;
                now_max_depth = flag;
            }
        }
    	if(count == 0)
        {
            fprintf(stderr, "Error: total %d, left %d, c[0] %d c[1] %d c[2] %d c[3] %d\n", get_real_depth(g_nodes[threadId][j].b), left, g_levs[threadId][num[1]].c[0], g_levs[threadId][num[1]].c[1], g_levs[threadId][num[1]].c[2], g_levs[threadId][num[1]].c[3]);
            print_base(stderr, g_nodes[threadId][j].b);
            exit(2);
        }
        if(count == 1)
        {
            g_nodes[threadId][j].seq[g_nodes[threadId][j].len] = dnaChar[now_max_i];
            g_nodes[threadId][j].len ++;
            g_nodes[threadId][j].seq[g_nodes[threadId][j].len] = '\0';
            g_nodes[threadId][j].level = num[1];
            set_base(g_nodes[threadId][j].b, sbase[threadId][now_max_i]);
            g_levs[threadId][num[1]].nnum++;

            if(max_depth < now_max_depth)
            {
                max_i = now_max_i; max_j = j;
                max_depth = now_max_depth;
            }else if(max_depth == now_max_depth)
            {   if(max_i != now_max_i)
                {
                    if(g_levs[threadId][ss[1]-1].node == j)
                    {
                        max_i = now_max_i; max_j = j;
                    }
                }
            }else if(max_i == now_max_i)
            {
                if(g_levs[threadId][ss[1]-1].node == j)
                {
                    max_j = j;
                }
            }
            if(j == g_levs[threadId][ss[1]-1].node)
            {
                node_max_i = now_max_i; node_max_j = j;
                node_max_depth = now_max_depth;
            }
            continue;
        }

        for(k=0; k < i; k++)
        {
            int real_depth = get_real_depth(sbase[threadId][k]);
            if(real_depth > 0)
            {
                set_base(g_nodes[threadId][num[0]].b, sbase[threadId][k]);
                g_nodes[threadId][num[0]].seq[0] = dnaChar[k];
                g_nodes[threadId][num[0]].seq[1] = '\0';
                g_nodes[threadId][num[0]].len = 1;
                g_nodes[threadId][num[0]].level = num[1];
                g_nodes[threadId][num[0]].parent = j;
                g_nodes[threadId][num[0]].depth = real_depth;
                g_nodes[threadId][num[0]].is_end = 0;
                g_nodes[threadId][num[0]].conf = 0;

                if(max_depth < real_depth)
                {
                    max_i = k; max_j = num[0];
                    max_depth = g_nodes[threadId][num[0]].depth;
                }else if(max_depth == real_depth)
                {    if(max_i != k)
                    {
                        if(g_levs[threadId][ss[1]-1].node == j && g_levs[threadId][ss[1]-1].node != g_nodes[threadId][max_j].parent)
                        {
                            max_i = k; max_j = num[0];
                        }
                    }
                }else{
                    if(max_i == k && g_levs[threadId][ss[1]-1].node == j && g_levs[threadId][ss[1]-1].node != g_nodes[threadId][max_j].parent)
                    {
                        max_j = num[0];
                    }
                }

                num[0]++;
                g_levs[threadId][num[1]].nnum++;
            }
        }
        
        if(now_noend_depth * 0.75 > now_max_depth)
        {
            g_nodes[threadId][j].conf = 1;
        }
    }

    //improve the treatment of heterozygous.
    if(ss[1] > 1 && num[0] == ss[0] && max_j != g_levs[threadId][ss[1]-1].node) //no change of node number.
    {
        if(node_max_i == max_i)
        {
            max_j = g_levs[threadId][ss[1]-1].node;
            max_depth = node_max_depth;
        }
    }

    if(g_levs[threadId][num[1]].nnum == 0)
    {
        if(g_levs[threadId][num[1]].c[4] > 0)
        {
            g_levs[threadId][num[1]].node = max_j;
            g_levs[threadId][num[1]].base = max_i;
            g_levs[threadId][num[1]].depth = 0;
            
            ss[0] = num[0];
            ss[2] = num[2];
            ss[1] = ss[1] + 1;
        }
        return 0;
    }   
    g_levs[threadId][num[1]].base = max_i;
    g_levs[threadId][num[1]].depth = max_depth;
    g_levs[threadId][num[1]].node = max_j;
    
    ss[0] = num[0];
    ss[2] = num[2];
    ss[1] = ss[1] + 1;

    return 1;
}

/*
    call consensus of the extended region.
*/

int pe_check_pure(ULL id, int size, int pos, int strand, int *real_size, int threadId)
{
    *real_size = 0;
    LL pid = get_pair_id(id);
if(DEBUG)
{
    fprintf(stdout, "Check %lld,%lld,%d,%d ", id, pid, strand, pos);
}
    if(strand == 0)
    {
if(DEBUG)  fprintf(stdout, "Minus\n");
        return 0;
    }
    if(is_used_bwt_sparse(pid, size) == 0)
    {
if(DEBUG)  fprintf(stdout, "NotUsed\n");
        return 0;
    }
    int i = tmpCount[threadId]-1;
    int lowerbound = pos - (InsertSize[size] - Read_length) - 3*SD[size];
    for(; i>=0; i--)
    {
        if(tmpRead[threadId][i].size != size)
            continue;
        if(tmpRead[threadId][i].pos < lowerbound)
        {
            if(DEBUG)  fprintf(stdout, "OutLowerBound\n");
            return 0;
        }
        if(tmpRead[threadId][i].id == pid)
        {
            *real_size = pos - tmpRead[threadId][i].pos + 1 + Read_length;
if(DEBUG)   fprintf(stdout, "Find %d ", *real_size);
            if(tmpRead[threadId][i].cons != 1 || tmpRead[threadId][i].unique != 1)
            {
if(DEBUG)       fprintf(stdout, "CONS%dUnique%d\n", tmpRead[threadId][i].cons, tmpRead[threadId][i].unique);
                return 0;
            }
            *real_size = pos - tmpRead[threadId][i].pos + 1 + Read_length;
if(DEBUG)   fprintf(stdout, "\n");
            return 1;
        }
    }
if(DEBUG)  fprintf(stdout, "NotFound\n");
    return 0;
}

int base2int(char b)
{
    if(b > 'T')
        b = b - ('T' - 't');

    switch (b)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        return 0;
    }
    return 0;
}

int get_subopt_brother(Node* nodes, int* bs, int start, int bi, int level, int num)
{
    int i, sub_i=-1;
    for(i = start-1; i>= 0; i--)
    {
        if(nodes[i].level < level)
            break;
        if(bs[i] == bi || i == bi)
        {
            if(sub_i == -1)
                sub_i = i;
            else
            {
                if(nodes[sub_i].depth < nodes[i].depth)
                    sub_i = i;
            }
        }
    }
    for(i=start +1; i< num; i++)
    {
        if(nodes[i].level - nodes[i].len > level)
            break;
        if(bs[i] == bi || i == bi)
        {
            if(sub_i == -1)
                sub_i = i;
            else{
                if(nodes[sub_i].depth < nodes[i].depth)
                    sub_i = i;
            }
        }
    }
    return sub_i;
}

void clean_inf(FILE* logs, Node* nodes, Lev* levs, ReadId* reads, int* s, int* bs, int bi)
{
    //only update the levs, s[1], reads.node
    int i, j; 
    for(i=0; i< s[2]; i++)
    {
        if(bs[reads[i].node] != bi)
        {
            reads[i].node = 0;
            reads[i].level = 0;
        }
    }

    //update s[1]
    for(i=1; i< s[1]; i++)
    {
        if(bs[levs[i].node] != bi)
        {
            j = get_subopt_brother(nodes, bs, levs[i].node, bi, i, s[0]);
            
            if(j == -1)
            {
                fprintf(stderr, "Cut directly %d->%d\n", s[1], i);
                s[1] = i;
                break;
            }
            int p = nodes[j].len - (nodes[j].level - i + 1) ;
            int newbase = base2int(nodes[j].seq[p]);
            int oldnode = levs[i].node;
            int oldbase = levs[i].base;
            levs[i].node = j;

            if(newbase == levs[i].base)
            {
                levs[i].depth = levs[i].c[levs[i].base] - levs[i].depth;
                levs[i].c[levs[i].base] = levs[i].depth;
                if(levs[i].depth > nodes[j].depth)
                {
                    levs[i].depth = nodes[j].depth; //would be overestimated.
                }
            }else{
                int k ;
                for(k=0; k<4; k++)
                {
                    if(k != newbase)
                        levs[i].c[k]=0;
                }
                levs[i].base = newbase;
                levs[i].depth = levs[i].c[levs[i].base];
                if(levs[i].depth > nodes[j].depth)
                    levs[i].depth = nodes[j].depth; //would be overestimated.
            }
            s[1]=i+1;
            break; //solve one error per time.

            if(levs[i].depth < Error_depth)
            {
                s[1] = i+1;
                break;
            }
        }else{
	    int k ;
	    for(k=0; k<4; k++)
	    {
	        if(k != levs[i].base)
	            levs[i].c[k]=0;
            }
            s[1]=i+1;
            break;
	}
    }

}

int get_ancestor(Node* nodes, int i)
{
   int nowi=i;
   while(nodes[nowi].parent != 0)
       nowi = nodes[nowi].parent;

   return nowi;
}

int fill_branches(Node* nodes, int n, int* bs)
{
    int branch_count=1;
    int i;
    for(i=1; i< n; i++)
    {
        if(nodes[i].parent != 0)
        {
            i--;
            break;
        }
        branch_count++;
        bs[i] = i;
    }

    for(; i< n; i++)
    {
        bs[i] = get_ancestor(nodes, i);
    }

    return branch_count;
}

void clean_path(FILE* logs, Node* nodes, Lev* levs, ReadId* reads, int* s, int ancestor )
{
    int i;
    int bs[s[0]];

    for(i=0; i< s[0]; i++)
        bs[i]=0;

    int branch_count = fill_branches(nodes, s[0], bs);

    clean_inf(logs, nodes, levs, reads, s, bs, ancestor);
}

int solve_conf(FILE* logs, Node* nodes, Lev* levs, ReadId* reads, int* s, int seed_len, int ctg_len, int threadId, int diff)
{
    int i, total_pe_count=0, real_size;
    int bs[s[0]];

    for(i=0; i< s[0]; i++)
        bs[i]=0;

    int branch_count = fill_branches(nodes, s[0], bs);
    int count[branch_count];
    int pecount[branch_count];

    for(i=0; i< branch_count; i++)
    {
        count[i]=0; pecount[i]=0;
    }

    for(i=0; i< s[2]; i++)
    {
        if(reads[i].node == 0)
            continue;
        count[bs[reads[i].node]]++;
        if(pe_check_pure(reads[i].id, reads[i].size, reads[i].level + ctg_len, reads[i].strand, &real_size, threadId) == 1)
        {
if(DEBUG)   fprintf(logs, "Got: %d %d %d %d\n", reads[i].id, reads[i].node, bs[reads[i].node], real_size);
            pecount[bs[reads[i].node]]++;
            total_pe_count++;
        }
    }

    int sup_i=-1, sup_c=0, c=0, total_c=0, total_pe_c = 0;
    for(i=1; i< branch_count; i++)
    {
        total_pe_c += pecount[i];
        total_c += count[i];
        if(count[i] >= Error_depth && pecount[i] > 0)
        {
            c++;
if(DEBUG)   fprintf(logs, "confs: %d:%d:%d\n", i, count[i], pecount[i]);
            //total_c += count[i];
            //total_pe_c += pecount[i];
            if(sup_c < pecount[i])
            {
                sup_c = pecount[i];
                sup_i = i;
            }
        }
    }

    if(c == 0)
        return 0;

if(DEBUG)    fprintf(logs, "CONF:%d,%d,%d,%d,%d,%d,%d;\n", c, branch_count, sup_i, sup_c, count[sup_i], total_c, total_pe_c);
    if( (c == 1 && (sup_c > Error_depth || (sup_c == Error_depth && total_pe_c == sup_c))&& sup_c >= count[sup_i]/10)//&& (sup_c >= 2 * Consensus_depth || (sup_c < 2 * Consensus_depth && sup_c == total_pe_count )))// && total_pe_count == sup_c)
        || (c > 1 && sup_c > Error_depth && sup_c >= count[sup_i]/5 && sup_c > double(total_pe_c) * 0.9 && total_pe_c - sup_c < Error_depth)
      )
    {
        clean_inf(logs, nodes, levs, reads, s, bs, sup_i);
        return 1;
    }

    //try to solve heterozygous bubble.
    if(Solve_Hete == 1 && diff < 5 && c == 2 && sup_c > 2 && seed_len <= depth_half_length && nodes[0].depth <= Upper_bound * double(Expect_depth*(Read_length-seed_len+1))/double(Read_length))
    {
        clean_inf(logs, nodes, levs, reads, s, bs, sup_i);
        return 1;
    }
    return 0;
}

void get_check_seq(char* seq, char* ctg_seq, Lev* levs, int lnum, int ctg_len)
{
    int i, j=0, k=0;

    for(i=ctg_len-1; i>=0 &&k < Read_length; k++, i--)
    {
        seq[k] = ctg_seq[i];
    }
    seq[k]='\0';

if(DEBUG) fprintf(stdout, "check seq:\n%s\n", seq);

}

int is_full_cons(ULL read_id, int size, int read_strand, int level, char* check_seq, int check_seq_len, int seed_len, int lnum)
{
    int start, end, startp;
    
    if(read_strand == 0)
    {
        start = level + seed_len -1;
        end = Read_length-1;
        startp = seed_len;
        return minus_extract_and_compare(idx2BWT[size], read_id, start, end, check_seq, check_seq_len, startp);
    }else{
        start = Read_length - (level + seed_len);
        end = 0;
        startp = seed_len + start;
        if(startp == Read_length)
            return 1;
        return plus_extract_and_compare(idx2BWT[size], read_id, start, end, check_seq, check_seq_len, startp);
    }
}

int is_full_consensus(ULL read_id, int size, int read_strand, int read_level, char* check_seq, int check_seq_len, int seed_len, int lnum)
{
    char seq[Read_length + 1];
    int qual[Read_length + 1];
//is_full_cons(read_id, size, read_strand, read_level, check_seq, check_seq_len, seed_len, lnum);
    get_id_seq(seq, qual, read_id, size);

    int i, j;
    if(read_strand == 0)
    {
        int rqual[Read_length + 1];
        char rcseq[Read_length + 1];
        reverse_com(seq, rcseq, Read_length);

        for(i=0; i< Read_length; i++)
            rqual[i] = qual[Read_length-i-1];

        i = read_level+seed_len-1; j = seed_len;
if(DEBUG)    fprintf(stdout, "Minu Read: %s, pos %d,%c; Cons %d,%c\n", rcseq, i, rcseq[i], j, check_seq[j]);

        for(; i< Read_length; i++, j++)
        {
            if(rcseq[i] != check_seq[j] && rqual[i] == 1)
            {
if(DEBUG)   fprintf(stdout, "i: %d,%c->%c,%d\n", i, check_seq[j], rcseq[i], rqual[i]);
                return 0;
            }
        }
    }else{

        i = read_level+seed_len-1; j = seed_len;
if(DEBUG)    fprintf(stdout, "Plus Read: %s, pos %d,%c; Cons %d,%c\n", seq, i, seq[i], j, check_seq[j]);

        for(; i< Read_length; i++, j++)
        {
            if(seq[i] != check_seq[j] && qual[i] == 1)
	    {
if(DEBUG)   fprintf(stdout, "i: %d,%c->%c,%d\n", i, check_seq[j], seq[i], qual[i]);
                return 0;
	    }
        }
    }
if(DEBUG)    fflush(stdout);
//    is_full_cons(read_id, size, read_strand, read_level, check_seq, check_seq_len, seed_len, lnum);
    return 1;
}

int mask_last_reads(FILE* logs, int cut_level, int* s, char* ctg_seq, int keep_used, int* keep_used_extension, int *has_used_read, int ctg_len, int iter, int unique, int threadId, int has_repeat_break, int seed_len)
{
    int i, j=tmpCount[threadId], used, cons, cut = cut_level, tid=-1, check_read=0, use_read=0;
    char check_seq[2*Read_length + 2];

    get_check_seq(check_seq, ctg_seq, g_levs[threadId], s[1], ctg_len);
    int check_seq_len = strlen(check_seq), cons_check;

    for(i=0; i< s[2]; i++)
    {
        if(i > 0 && g_reads[threadId][i].node == 0 && g_reads[threadId][i].level == 0)
            continue;
        if(g_reads[threadId][i].level > cut_level)
            break;

        if(j >= MaxReadCount)
        {
            fprintf(stderr, "Too many reads %d\n", j);
	        return -3;
        }
        tmpRead[threadId][j].id     = g_reads[threadId][i].id;
        tmpRead[threadId][j].strand = g_reads[threadId][i].strand;
        tmpRead[threadId][j].size   = g_reads[threadId][i].size;
        tmpRead[threadId][j].pos    = g_reads[threadId][i].level - 1 + ctg_len; //start from 1 for level.
        tmpRead[threadId][j].iter   = iter;
        tmpRead[threadId][j].unique = unique;
        tmpRead[threadId][j].cons   = g_reads[threadId][i].cons;
        tmpRead[threadId][j].keep = 0;

        cons_check = double(Read_length)/double(g_reads[threadId][i].level + seed_len) > Upper_bound ? 1 : 0;
if(DEBUG)        fprintf(logs, "read %d cons %d pos %d \n", g_reads[threadId][i].id, g_reads[threadId][i].cons, tmpRead[threadId][j].pos);
        if(g_reads[threadId][i].cons == 1 && tmpRead[threadId][j].pos > Read_length)
        {
            check_read++;
            //if(has_repeat_break == 1 && 
            if(cons_check == 1 && is_full_cons(g_reads[threadId][i].id, g_reads[threadId][i].size, g_reads[threadId][i].strand, g_reads[threadId][i].level, check_seq, check_seq_len, seed_len, s[1]) == 0)
            {
                tmpRead[threadId][j].cons = 0;
                g_reads[threadId][i].cons = 0;
                continue;
            }
            use_read++;

            if(is_used_bwt_sparse(g_reads[threadId][i].id, g_reads[threadId][i].size) == 0)
            {
                endInf[threadId].ur ++;
                tmpRead[threadId][j].keep = 1;
if(DEBUG)       fprintf(logs, "read %d cons %d pos %d \n", g_reads[threadId][i].id, g_reads[threadId][i].cons, tmpRead[threadId][j].pos);
                check_and_set_bwt_sparse_thread_and_used(g_reads[threadId][i].id, g_reads[threadId][i].size, threadId);
                //    fprintf(stderr, "UR %d %lld %d %d\n", threadId, reads[i].id, is_used_bwt_sparse(reads[i].id, reads[i].size), get_threadId_bwt_sparse(reads[i].id, reads[i].size));
            }else{
                if(cons_check == 0 && is_full_cons(g_reads[threadId][i].id, g_reads[threadId][i].size, g_reads[threadId][i].strand, g_reads[threadId][i].level, check_seq, check_seq_len, seed_len, s[1]) == 0)
		{
                    g_reads[threadId][i].cons = 0;
                    tmpRead[threadId][j].cons = 0;
                    use_read--;
                    continue;
		}
                tid = get_threadId_bwt_sparse(g_reads[threadId][i].id, g_reads[threadId][i].size); 
                if(tid != -1 && tid != threadId)
                {
                    tmpCount[threadId] = j;
                    //fprintf(stderr, "CONFMeet:%d %d %d %d %d %d\n", reads[i].id, reads[i].size, threadId, j, tid, is_used_bwt_sparse(reads[i].id, reads[i].size));
                    free_thread_tmpRead(threadId, 1);
                    return -100;
                }
if(DEBUG)       fprintf(logs, " (%d,%d,%d)",g_reads[threadId][i].id, g_reads[threadId][i].level, tid);
                if(tid == threadId)
                    *keep_used_extension = 1;

                cut = g_reads[threadId][i].level - 1;
                if(cut <= 0)
                {
                    tmpCount[threadId] = j;
                    return cut;
                }
                *has_used_read = 1;
                break;
            }
        }
        j++;
    }
    tmpCount[threadId] = j;

    if(use_read > 0)
    {
        if(four_iter[threadId] > 0)
        {
             four_iter[threadId] = 0; four_iter_checked_read[threadId]=0; four_iter_used_read[threadId] = 0;
        }
    }else
    {
        if(check_read > 5)
        {
if(DEBUG)   fprintf(stdout, "Maybe Same seed problem. check read %d and use read %d\n", check_read, use_read);
            return -2; 
        }

        four_iter[threadId] ++;
        four_iter_checked_read[threadId] += check_read;
        if(four_iter[threadId] >= 4 && four_iter_checked_read[threadId] >= 10 && four_iter_used_read[threadId] == 0)
        {
if(DEBUG)   fprintf(stdout, "Maybe Same seed problem. iter check number %d check read %d and use read %d\n", four_iter[threadId], check_read, use_read);
            return -2;
        }
        if(four_iter[threadId] >= 10)
        {
if(DEBUG)   fprintf(stdout, "Maybe Same seed problem. iter check number %d check read %d\n", four_iter[threadId], four_iter_checked_read[threadId], use_read);
            return -2;
        }
    }
    

    //fit for the cut level. 
    if(i < s[2])
    {
        s[2] = i;
    }

    for(j=s[0]-1; j>=0; j--)
    {
        if(g_nodes[threadId][j].level <= cut)
            continue;
        //fprintf(stderr, "Cut: %d: %d, level%d, len %d\n", j, cut, nodes[j].level, nodes[j].len);

        int left = int(g_nodes[threadId][j].level) - cut ;
        if(left <=0)
            continue;
        if(left >= g_nodes[threadId][j].len)
            left = g_nodes[threadId][j].len;
        g_nodes[threadId][j].len -= left;
        g_nodes[threadId][j].seq[ g_nodes[threadId][j].len ]='\0';
	//fprintf(stderr, "[%d;%d]",j,left);
    }
    s[1] = cut + 1;
    
   return cut; 
}

int is_quick_decrease(int d1, int d2)
{
    if(d2 < 5)
        return d2*2-1 <= d1;
    return d2*1.5 <= d1;
}

int has_branches(Lev &lev)
{
    return lev.c[0] + lev.c[1] + lev.c[2] + lev.c[3] > lev.depth;
}


int get_sub_node(Node* nodes, int level, int p_node, char base, int num)
{
    int i, max_i=-1, max=0, high_num=0;

    for(i=p_node+1; i< num; i++)
    {
        if(nodes[i].level < level-1)
            continue;
        if(nodes[i].level - nodes[i].len == level-1 && nodes[i].seq[0] != base) 
        {
            if(nodes[i].depth > Error_depth)
                high_num++;
            if(nodes[i].depth > max)
            {
                max = nodes[i].depth;
                max_i = i;
            }
        }
    }

    if(high_num > 1)
        return 0;

    return max_i; 
}

void get_sub_seq(Node* nodes, Lev* levels, int p, int num, char* subseq)
{
    int type[num], max_depth[num], max_i[num], max_same[num], type_count[num], nodes_count[num];
    int i, j, k;
    for(i=0; i< num; i++)
    {
        type[i]=0;
        type_count[i]=0;
        nodes_count[i]=0;
        max_depth[i]=0;
        max_i[i]=0;
        max_same[i]=0;
    }

    type[p]=1;
    for(i=p+1; i< num; i++)
    {
        if(type[nodes[i].parent] == 1)
        {
           type[i]=1;
           type_count[nodes[i].parent]++;
           nodes_count[nodes[i].parent] += nodes[i].depth;
           if( 
                  ( 
                     (nodes[i].depth > max_depth[nodes[i].parent])
                     || (nodes[i].depth == max_depth[nodes[i].parent] && max_same[nodes[i].parent] == 0 && nodes[i].seq[0] == dnaChar[levels[nodes[i].level +1 - nodes[i].len].base])
                     || (nodes[i].depth == max_depth[nodes[i].parent] && max_same[nodes[i].parent] == 0 && nodes[i].len > nodes[ max_i[nodes[i].parent] ].len)
                  )
             )
           {
               max_depth[nodes[i].parent] = nodes[i].depth;
               max_i[nodes[i].parent] = i;
               if(nodes[i].seq[0] == dnaChar[levels[nodes[i].level +1 - nodes[i].len].base])
                   max_same[nodes[i].parent] = 1;
           }
        }
    }

    j=0;
    for(i=p; i< num; i++)
    {
        if(i > p && type[i] == 1 && (type[nodes[i].parent] == 0 || max_i[nodes[i].parent] != i))
            type[i] = 0;

        if(type[i] == 0 )
            continue;
        if(i>p && type_count[i] > 0 && type_count[i] == nodes_count[i])
        {
            //fprintf(stdout, "i %d: type count %d, node count %d\n", i, type_count[i], nodes_count[i]);
            break;
        }
        //fprintf(stdout, ">%d,%d", i, nodes[i].len);
        for(k=0; k < nodes[i].len; k++, j++)
            subseq[j] = nodes[i].seq[k];
        if(nodes[i].conf == 1)
            break;
    }
    //fprintf(stderr, "\n");
    subseq[j]='\0';

}

int get_first_branch(Node* nodes, int num, int now_branch_node, int current_node, char* s1)
{
    int flag[num], type[num];
    int i;
    flag[current_node]=1;
    for(i=0; i< num; i++)
    {
        flag[i]=0; type[i]=0;
    }
    int node=current_node;
    while(!(node == 0 && nodes[node].parent ==0) )
    {
        node = nodes[node].parent;
        flag[node] = 1;
    }

    node = now_branch_node;
    while(flag[node] == 0)
    {
        type[node]=1;
        if(flag[nodes[node].parent] == 1)
            break;
        node = nodes[node].parent;
    }
    int j=0, k=0;
    for(i=node; i< num; i++)
    {
        if(type[i] == 0 )
            continue;
        for(k=0; k < nodes[i].len; k++, j++)
        {
            s1[j] = nodes[i].seq[k];
        }
    }
    s1[j]='\0';
    return node;
}

int get_branch_node(FILE* logs, Node* nodes, Lev* levs, int p, int b, int num)
{
    int p_node = nodes[levs[p].node].parent;
    

    if(p_node != levs[p-1].node)
        return 0;
    
    int i;
    char base = dnaChar[b];
    for(i=p_node+1; i< num; i++)
    {
        if(nodes[i].level < p-1 || nodes[i].parent != p_node)
            continue;
        if(nodes[i].level - nodes[i].len + 1 == p && nodes[i].seq[0] == base) 
        {
            return i;
        }
    }
    return 0;
}


int get_sub_node_id2(FILE* logs, Node* nodes, Lev* levs, int p, int b, int* s, char* s1)
{
    int p_node = nodes[levs[p].node].parent;
    int i, max_i=-1, max=0, high_num=0;
    char base = dnaChar[b];//dnaChar[levs[p].base];

    for(i=1; i< s[0]; i++)
    {
        //fprintf(stderr, "%d %d %d %d %d\n", i, p, nodes[i].level, nodes[i].parent, nodes[nodes[i].parent].level);
        if(nodes[i].level < p || i == levs[p].node || nodes[nodes[i].parent].level > p - 1)
            continue;
        //fprintf(stderr, "%d %d %d %c \n", i, p, nodes[nodes[i].parent].level, nodes[i].seq[p - nodes[nodes[i].parent].level -1]);
        if(nodes[i].seq[p - nodes[nodes[i].parent].level -1] == base)
        {
            if(nodes[i].depth > Error_depth)
                high_num++;
            if(nodes[i].depth > max)
            {
                max = nodes[i].depth;
                max_i = i;
            }
        }
    }
    int n_node = max_i; //now branch node.
    if(n_node < 0)
    {
       fflush(logs);
       fprintf(stderr, "Subnode: %d, %d, %c\n", p_node, p, base);
       exit(1);
    }

    int c_node = levs[p].node; //current node.
    int fb_node = get_first_branch(nodes, s[0], n_node, c_node, s1); //first branch node.

    //fprintf(logs, "nobranch %d, nonode %d, firstbranch %d\n", n_node, c_node, fb_node);
    return fb_node;
}

int get_sub_node_id(Node* nodes, Lev* levs, int p, int* s)
{
    int p_node = levs[p-1].node;
    if(nodes[p_node].level != p-1)
        return 0;

    int node = get_sub_node(nodes, p, p_node, dnaChar[levs[p].base], s[0]);
    if(node < 0)
        return 0;
    return node;
}

int get_pair_seq(FILE* logs, Node* nodes, Lev* levs, int p, int* s, char* s1, char* s2, int sub_node)
{

    if(strlen(s1) == 0)
        get_sub_seq(nodes, levs, sub_node, s[0], s1);

    int i, j, dif=0, node_id, k;

//    fprintf(logs, "s1 %s\n", s1);
    for(i=p, j=0; i< s[1] && j<strlen(s1); i++, j++)
    {
        if( levs[i].base == 4)
            break;
	if(levs[i].node != levs[i-1].node && (nodes[levs[i].node].parent != levs[i-1].node || (nodes[levs[i-1].node].conf == 1 && i > p)))
	{
	    //dif++;

            if(i <= nodes[levs[i-1].node].level)
            {
                node_id = levs[i-1].node;
                
                for(k = i; k<= nodes[node_id].level && j < strlen(s1); k++, i++, j++)
                {
                    s2[i-p] = nodes[node_id].seq[nodes[node_id].len - (nodes[node_id].level -k)-1];
                    if(s2[i-p] != s1[j])
                        dif++;
                }
            }
	    break;
	}
        s2[i-p] = dnaChar[ levs[i].base ];
        if(s2[i-p] != s1[j])
            dif++;
    }
    s2[i-p]='\0';
    s1[i-p]='\0';
    //fprintf(logs, "%s, %s, %d, %d\n", s1, s2, p, sub_node);
    return dif;
}



int has_branches2(int level, Node* nodes, int nnum, int node)
{
    int i=0; 
    for(i=0; i< nnum; i++)
    {
        if(i == node)
            continue;
        if(nodes[i].level < level)
            continue;
        if(nodes[i].level - nodes[i].len >= level)
            break;
        if(nodes[i].level >= level && nodes[i].level - nodes[i].len < level && nodes[i].depth > 1)
            return nodes[i].level - nodes[i].len;
    }
    return 0;
}

int quick_decrease(Node* nodes, Lev* levs, int p, int seed_len)
{
    double base_depth = double(nodes[0].depth) / double(Read_length - seed_len + 1 ) * Read_length;
    double lev_dec = base_depth / double(Read_length);
    double node_dec = double(nodes[levs[p].node].len - (nodes[levs[p].node].level - p)) * lev_dec;
    int real_lev_dec = levs[p-1].depth - levs[p].depth;
    int real_node_dec = nodes[levs[p].node].depth - levs[p].depth;

    //fprintf(logs, " %.2f:[%.2f,%d;%.2f,%d]\n", base_depth, lev_dec, real_lev_dec, node_dec, real_node_dec);
    //decreasing too fast.
    if(real_lev_dec - lev_dec > Error_depth*2 )
         return 1;
    if(real_node_dec - node_dec > Error_depth*2 )// || (levs[p].depth <= 2* Error_depth && real_node_dec - node_dec > real_node_dec/2))
         return 1;
    
    //extending single node too long.
    if(p + seed_len >= depth_9_length && levs[p].depth < int(nodes[levs[p].node].depth/3))
        return 1;
    if(p + seed_len >= depth_6_length && levs[p].depth <= 2* Error_depth && levs[p].depth < int(nodes[levs[p].node].depth/2))
        return 1;
    
    return 0;
}

int is_short_cons(int len, int seed_len)
{
    return (len < seed_len/2 && len < MINOVERLAP) ? 1 : 0;
}

int get_consist_seq(int ii, int kk, char* s1, char* s2, char* seq, int seed_len)
{
    int j;
    int i = ii, k = kk;
    for(j=0; i>=1 && k>=1; i--, j++, k--)
    {
        if(s1[i] != s2[k])
        {
            if(is_short_cons(j, seed_len) == 1)
            {
                j=-1;
                continue;
            }
            break;
        }
        seq[j] = s1[i];
    }
    seq[j]='\0';
    return j;
}

int is_hete_bubble(char* s1, char* s2, int now_base_depth, int diff, int seed_len, int threadId)
{
    int nlen = strlen(s1);
    char seq[nlen + 1];

    int i=nlen-1, j, k=i;
    j = get_consist_seq(i, k, s1, s2, seq, seed_len);

    if(is_short_cons(j, seed_len) == 1)
    {
        if(s1[i] == s2[k-1])
        {   k--;
            j = get_consist_seq(i, k, s1, s2, seq, seed_len);
            k++;
        }

        if(is_short_cons(j, seed_len) == 1 && s1[i-1] == s2[k])
        {
           i--;
           j = get_consist_seq(i, k, s1, s2, seq, seed_len);
        }
    }else if(diff > 4)
        return 0;

    if(is_short_cons(j, seed_len) == 1)
        return 0;
    return dual_bubble_check(seq, now_base_depth, threadId);
}

int solve_by_seed(FILE* logs, char* cseq, char* st, char* s1, char* s2, int p, double d1, int d2, double diff, int seed_len, int max_level, int threadId)
{
    int nlen = p -1 ;
    int tlen = Read_length - 1 - nlen;//strlen(s1);

    int i,j, slen = strlen(s1);
    for(i=1; i< slen; i++)
    {
        if(s1[i] != s2[i])
            break;
        if(g_levs[threadId][p+i-1].c[ charMap[s1[i]] ] < Error_depth)
        {
            i = slen;
            break;
        }
    }

    if(i == slen || i+p-1 >= max_level)
        slen = 1;
    else
        slen = i + 1;

//    fprintf(logs, "i %d, p %d, max_level %d, slen %d\n", i, p, max_level, slen);
    tlen = Read_length - nlen - slen;
    if(strlen(cseq) < tlen)
        tlen = strlen(cseq);

    char seq1[tlen + nlen + slen+1];
    char seq2[tlen + nlen + slen+1];

    for(i=0, j=strlen(cseq)-tlen; i<tlen; i++, j++)
    {
        seq1[i] = complementMap[ cseq[j] ];
        seq2[i] = complementMap[ cseq[j] ];
    }

    //st starting from 1.
    for(j=1; j< p; j++, i++)
    {
        seq1[i] = complementMap[ st[j] ];
        seq2[i] = complementMap[ st[j] ];
    }
    
    for(j=0; j< slen; j++, i++)
    {
        seq1[i] = complementMap[ s1[j] ];
        seq2[i] = complementMap[ s2[j] ];
    }
    seq1[i]='\0'; seq2[i]='\0';
    //fprintf(logs, "seq1: %s seq2: %s\n", seq1, seq2); 
    int type = dual_seq_check(logs, seq1, seq2, nlen + seed_len + slen, d1, strlen(s1), diff, threadId);
    return type;
    
}

int is_long_repeat_seed(int depth, int len);

int high_confident_err(int seed_len, int depth0, int p, int prev_depth, int p_depth, int c_depth, int diff)
{
    if(seed_len < depth_half_length && seed_len + p < depth_9_length
        && depth0 < Expect_depth*(Read_length-seed_len+1)/Read_length  && c_depth == 1 
        && diff == 1 && prev_depth - p_depth < Error_depth)
        return 1;
    return 0;
}

int set_reads_cons(int *s, int p, int flag, int threadId)
{
    int i;
    for(i=0; i< s[2]; i++)
    {
        if(g_reads[threadId][i].level == p && g_reads[threadId][i].node == g_levs[threadId][p-1].node && g_reads[threadId][i].cons != flag)
        {
            //fprintf(stdout, "set cons %d\n", g_reads[threadId][i].id);
            g_reads[threadId][i].cons = flag;
        }
    }
}

int is_risk_tag(char* seed_seq, int flag)
{
    int i, j;
    if(flag == 1)
    {
        if(seed_seq[0] == 'C' && seed_seq[1] == 'G' && seed_seq[2] == 'G')
            return 1;
        for(i=0; i< 10; i++)
        {
            if(seed_seq[i] == 'C')
            {
                for(j=i+1; j<12; j++)
                {
                    if(seed_seq[j] != 'G')
                        break;
                }
                if(j-i > 2)
                    return 1;
            }
        }
        return 0;

    }else{
        if(seed_seq[0] == 'G' && seed_seq[1] == 'C' && seed_seq[2] == 'C')
            return 1;
        for(i=0; i< 10; i++)
        {
            if(seed_seq[i] == 'G')
            {
                for(j=i+1; j<12; j++)
                {
                    if(seed_seq[j] != 'C')
                        break;
                }
                if(j-i > 2)
                    return 1;
            }
        }
        return 0;
    }
    return 0;
}

int is_strand_bias(ULL plus, ULL minus)
{
    if((plus == 0 && minus > 3* Consensus_depth)|| (plus > 0 && plus < Error_depth/2 && plus * 3* Consensus_depth < minus))
        return 1;

    if((minus == 0 && plus > 3 * Consensus_depth)|| (minus > 0 && minus< Error_depth/2 && minus *3* Consensus_depth < plus))
        return 2;

    return 0;
}
int is_illumina_error(FILE* logs, char* seed_seq, Node &n, Node &nc, int now_level_depth, int branch_depth)
{
    int l = strlen(seed_seq);
    int i, j;

    int risk_plus = is_risk_tag(seed_seq, 1);
    int risk_minus = is_risk_tag(seed_seq, 0);

    if(risk_plus == 0 && risk_minus == 0)
        return 0;

    ULL plus_depth = get_strand_depth(n.b, 1);
    ULL minus_depth = get_strand_depth(n.b, 0);
  
    ULL cons_plus_depth = get_strand_depth(nc.b, 1);
    ULL cons_minus_depth = get_strand_depth(nc.b, 0);

    int strand_bias = is_strand_bias(plus_depth, minus_depth);
    int cons_strand_bias = is_strand_bias(cons_plus_depth, cons_minus_depth);

//    fprintf(logs, "risk plus %d risk minus %d; strand bias %d cons strand bias %d, level depth %d branch depth %d; nc depth %d nc len %d.\n", risk_plus, risk_minus, strand_bias, cons_strand_bias, now_level_depth, branch_depth, nc.depth, nc.len);

    if(cons_strand_bias == strand_bias)
        return 0;

    int risk_conf=0, risk_cons=0;
    if((risk_plus == 1 && strand_bias == 2) || (risk_minus == 1 && strand_bias == 1))
        risk_conf = 1;

    if((risk_plus == 1 && cons_strand_bias == 2) || (risk_minus == 1 && cons_strand_bias == 1))
        risk_cons = 1;

    if(risk_conf == 1 && risk_cons == 1)
    {
        if((plus_depth + minus_depth) == 1 && branch_depth > Error_depth && now_level_depth - (cons_plus_depth + cons_minus_depth) < Error_depth/2)
            return 2;
        return 1;
    }
    if(risk_conf == 1)
        return 1;

    if(risk_cons == 1 && now_level_depth - (cons_plus_depth + cons_minus_depth) < Error_depth/2 && (now_level_depth > Error_depth || nc.len == 1))
        return 2;

    return 0;
}

int find_last_lev(FILE* logs, int* s, char* ctg_seq, char* seed_seq, int seed_len, int ctg_len, int initial_depth, int* is_unique, int* has_repeat_break, int threadId)
{
    int i, j, k, is_depth = 0, has_other =0, has_other_depth=0;
    int solved[s[1]];
    
    int back_len = 0, sub_node=0, sub_node_level=0;
    char s1[Read_length], s2[Read_length], st[Read_length];
    int is_tip, is_solved;
    int is_low_depth = initial_depth < (double(Read_length-seed_len+1)*Expect_depth)/double(Read_length)/2.5 ? 1 : 0;
    int is_high_depth = initial_depth > (double(Read_length-seed_len+1)*Expect_depth)/double(Read_length) ? 1 : 0;
    double current_base_depth = double(initial_depth * Read_length)/double(Read_length-seed_len+1);
    if(current_base_depth > Expect_depth)
        current_base_depth = Expect_depth;

    for(i=1; i<s[1]; i++)
    {
        solved[i] = 1;
        double curr_depth = double(current_base_depth*(Read_length-(seed_len+i)+1))/double(Read_length);
        //checking no-branch issue.
        if(g_levs[threadId][i].depth < Error_depth //low depth.
           || (g_levs[threadId][i].nnum <= 3 && 3 * g_levs[threadId][i].depth <= initial_depth) //too many branches with high depth.
           || ( i>1
                && (
                     (g_levs[threadId][i].node != g_levs[threadId][i-1].node && g_nodes[threadId][ g_levs[threadId][i].node ].parent != g_levs[threadId][i-1].node) //path changed.
                     || (g_levs[threadId][i].nnum == 1 && curr_depth < 2 * Expect_depth && g_levs[threadId][i].depth < 2 * Error_depth) //low depth region, base by base extension.
                     || (g_levs[threadId][i].nnum <= 2 &&  //has confliction here.
                           ( quick_decrease(g_nodes[threadId], g_levs[threadId], i, seed_len) == 1 //quick decrease of depth.
                            || (g_levs[threadId][i].depth <= 3 * Error_depth && curr_depth <= 3*Error_depth )  //depth too low to solve this branch easily.
                            || (g_levs[threadId][i].depth <= 3 * Error_depth && 2 * g_levs[threadId][i].depth <= initial_depth ) //clear branch and low depth.
                            || (is_low_depth == 1 && 2 * g_levs[threadId][i].depth <= initial_depth)
                           )
                         )

                   )
              )
          )
        {
            is_depth = 1;
            break;
        }
        int changed_consensus = 0;
        //checking branching issue.
        for(j=0; j<4; j++)
        {
            if(g_levs[threadId][i].c[j] == 0 || j == g_levs[threadId][i].base)
                continue;
            //having one branch.
            solved[i] = 0;
if(DEBUG)   fprintf(logs, "seed_len %d, i %d, total len %d, depth %d, conf_depth %d, tobase %d, cur depth %.1f\n", seed_len, i, seed_len + i, g_levs[threadId][i].depth, g_levs[threadId][i].c[j], j, curr_depth);

            //long extension break, reduce mismatches.
            if( i> 1 && (g_levs[threadId][i].depth + g_levs[threadId][i].c[j]) * 2 < initial_depth && g_levs[threadId][i].depth <= 3* Error_depth)
            {
if(DEBUG)       fprintf(logs, "Decrease to low.\n");
                break;
            }
  
            //careful solving cases, increase length.
            if( is_low_depth == 1  && i > 1 && initial_depth - g_levs[threadId][i].depth > g_levs[threadId][i].depth/2 && initial_depth - g_levs[threadId][i].depth > 2*Error_depth && g_levs[threadId][i].c[j] >= Error_depth)
            {
if(DEBUG)       fprintf(logs, "Branches in low coverage regions.\n");
                is_depth = 1;
                break;
            }      

            //having been solved.
            if(g_nodes[threadId][ g_levs[threadId][i].node ].level - g_nodes[threadId][ g_levs[threadId][i].node ].len + 1 != i)
            {
if(DEBUG)       fprintf(logs, "Innode.\n");

                if(curr_depth <= 2* Error_depth || g_levs[threadId][i].depth <= 2*Error_depth)
                    break;

                solved[i] = 1;
                continue;
            }
            
            has_other = 0; has_other_depth=0;
            for(k=j+1; k< 4; k++)
            {
                if(k == g_levs[threadId][i].base)
                    continue;
                if(g_levs[threadId][i].c[k] >= Error_depth)
                {
                    has_other = 1; has_other_depth += g_levs[threadId][i].c[k];
                }
            }
            if(has_other == 1 && has_other_depth + g_levs[threadId][i].c[j] > g_levs[threadId][i].depth/2 && curr_depth < Expect_depth/5 && curr_depth < current_base_depth/3)
            {
if(DEBUG)       fprintf(logs, "Too much confliction and too long extension.\n");
                break;
            }
            //get current branch caused sub-node.
            s1[0]='\0';
            sub_node = get_branch_node(logs, g_nodes[threadId], g_levs[threadId], i, j, s[0]);
            
            
            //no new sub-node from the main branch.
            if(sub_node == 0)
            {
                s1[0]=dnaChar[j]; s1[1]='\0';
                s2[0]=dnaChar[g_levs[threadId][i].base]; s2[1]='\0';
                is_solved = solve_by_seed(logs, ctg_seq, st, s1, s2, i, current_base_depth, g_levs[threadId][i].depth, 0, seed_len, g_nodes[threadId][g_levs[threadId][i].node].level, threadId);
if(DEBUG)       fprintf(logs, "No sub node.\n");

                if(is_solved == -1 && has_other == 1)
                {
                    changed_consensus = 1;
if(DEBUG)           fprintf(logs, "ChangeCons3 at %d with num %d\n", i, g_levs[threadId][i].nnum);
                    set_reads_cons(s, i, 0, threadId);
                    g_levs[threadId][i].base = j;
                    g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                    g_levs[threadId][i].node = sub_node;
                    continue;
                }

                if(is_solved != 0 && is_solved < 3)
                {
                    if(is_solved != 2)
                    {
                        *is_unique = 0;
                        *has_repeat_break = 1;
                        if(is_solved == 1 && i > 1 && g_levs[threadId][i].c[j] == 1 && g_levs[threadId][i].depth <= 3 * Consensus_depth)
                            break;
                        if(is_solved == -1)
                        {
if(DEBUG)                   fprintf(logs, "ChangeCons1 at %d with num %d\n", i, g_levs[threadId][i].nnum);
                            set_reads_cons(s, i, 0, threadId);
                            g_levs[threadId][i].base = j;
                            g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                            g_levs[threadId][i].node = sub_node;
                            //set_reads_cons(s, i, 1, threadId);
                            return i;
                        }
                    }
                
                    solved[i] = 1;
                    continue;
                }
                break;
            }
            
            if(i == 1 && Expect_depth > Error_depth * 10 && g_levs[threadId][i].depth + g_levs[threadId][i].c[j] < curr_depth * Upper_bound)
            {
                //if(is_high_depth_sequencing_error(logs, g_nodes[threadId], g_levs[threadId], i, s, sub_node) == 1)
                //print_base(logs, g_nodes[threadId][sub_node].b);
                //print_base(logs, g_nodes[threadId][g_levs[threadId][i].node].b);
                int i_err = is_illumina_error(logs, seed_seq, g_nodes[threadId][sub_node], g_nodes[threadId][g_levs[threadId][i].node], g_levs[threadId][i].depth, g_levs[threadId][i].c[j]);
                if(i_err == 1)
                {
if(DEBUG)           fprintf(logs, "Illumina sequencing error: %s.\n", seed_seq);
                    solved[i] = 1;
                    continue;
                }
                if(i_err == 2 && is_low_depth == 0)
                {
if(DEBUG)           fprintf(logs, "ChangeCons3 at %d with num %d seed %s.\n", i, g_levs[threadId][i].nnum, seed_seq);
                    set_reads_cons(s, i, 0, threadId);
                    g_levs[threadId][i].base = j;
                    g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                    g_levs[threadId][i].node = sub_node;
                    return i;
                }
            }
 
            //get the sequences of two branches, one is main, the other is the new branch.
            int diff = get_pair_seq(logs, g_nodes[threadId], g_levs[threadId], i, s, s1, s2, sub_node);
            if(strlen(s1) == 0)
            {
if(DEBUG)       fprintf(logs, "treat as sequencing error.\n");
                solved[i] = 1;
                continue;
            }
            double dif_ratio = double(diff)/double(strlen(s1));
if(DEBUG)   fprintf(logs, "s1 %s, s2 %s diff %d len %d\n", s1, s2, diff, strlen(s1));


            if(dif_ratio < 0.2 && high_confident_err(seed_len, initial_depth, i, g_levs[threadId][i-1].depth, g_levs[threadId][i].depth, g_levs[threadId][i].c[j], diff) == 1)
            {
                solved[i] = 1;
                continue;
            }
if(DEBUG)   fprintf(logs, "lowconf\n");

            //double curr_depth = double(Expect_depth*(Read_length-(seed_len+i)+1))/double(Read_length);
            if(g_levs[threadId][i].c[j] <= Error_depth && g_levs[threadId][i].depth > 2*Error_depth && curr_depth > 2 * Error_depth
                && g_levs[threadId][i].c[j] < g_levs[threadId][i].depth * 0.1 
                && !((seed_len > depth_half_length && (curr_depth * 2 < g_levs[threadId][i].depth || (dif_ratio > 0.05))))// || (dif_ratio > 0.1 )))
                && !(diff >= 5)
              )
            {
                solved[i] = 1;
                continue;
            }
if(DEBUG)   fprintf(logs, "NotLargeDif&LowBranch\n");
	    
            //checking whether this is a tip or solved by increasing seed length.
            is_solved = solve_by_seed(logs, ctg_seq, st, s1, s2, i, current_base_depth , g_levs[threadId][i].depth, dif_ratio, seed_len, g_nodes[threadId][g_levs[threadId][i].node].level, threadId);

            if(is_solved == -1 && has_other == 1)
            {
                changed_consensus = 1;
                if(DEBUG)  fprintf(logs, "ChangeCons4 at %d with nnum %d \n", i, g_levs[threadId][i].nnum);
                set_reads_cons(s, i, 0, threadId);
                g_levs[threadId][i].base = j;
                g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                g_levs[threadId][i].node = sub_node ;
                continue;
            }
            if(is_solved != 0 && is_solved < 3 && !(is_solved == -1 && has_other == 1))
            {
                if(is_solved != 2)
                {
                    *is_unique = 0;
                    //break;
                    *has_repeat_break = 1;
                    if(is_solved == 1 && i > 1 && g_levs[threadId][i].c[j] == 1 && g_levs[threadId][i].depth <= 3 * Consensus_depth)
                        break;
                    if(is_solved == -1)
                    {
if(DEBUG)               fprintf(logs, "ChangeCons2 at %d with nnum %d \n", i, g_levs[threadId][i].nnum);
                        set_reads_cons(s, i, 0, threadId);
                        g_levs[threadId][i].base = j;
                        g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                        g_levs[threadId][i].node = sub_node ;
                        //set_reads_cons(s, i, 1, threadId);
                        return i;
                    }
                }
                solved[i] = 1;
                continue;
            }

            if(is_solved == 5 && has_other == 1) 
            {
                solved[i] = 1;
                continue;
            }
if(DEBUG)   fprintf(logs, "Unsolved.\n");

            if(Solve_Hete == 1 && curr_depth > 2 * Error_depth //&& is_low_depth == 0
               && (g_levs[threadId][i].depth + g_levs[threadId][i].c[j] <= Upper_bound * curr_depth || strlen(s1) > seed_len * 2)
               && is_short_cons(strlen(s1), seed_len) == 0
              )
            {
                if(is_hete_bubble(s1, s2, current_base_depth, diff, seed_len, threadId) == 1)
                {
                    solved[i] = 1;
                    continue;
                }else if(diff == 1 && i+ seed_len <= MINOVERLAP + 1 && strlen(s1) >= MINOVERLAP && g_levs[threadId][i].depth + g_levs[threadId][i].c[j] <= Upper_bound * curr_depth)
                {
                    solved[i] = 1;
                    continue;
                }
            }

if(DEBUG)   fprintf(logs, "NotHete\n");
	     //masked at Sep 18.2014 for YH assembly.
            if(Solve_Hete == 1 
               && (is_solved >= 3 && is_solved <5) 
               && ((diff < 3 && dif_ratio < 0.2) || (diff< 5 && dif_ratio <= 0.05))
               && g_levs[threadId][i].depth + g_levs[threadId][i].c[j] <= curr_depth
              )
            {
                if(is_solved == 4 && (g_levs[threadId][i].c[j] == g_levs[threadId][i].depth || g_levs[threadId][i].c[j] > Error_depth))
                {
if(DEBUG)           fprintf(logs, "ChangeCons3 at %d with nnum %d \n", i, g_levs[threadId][i].nnum);
                    set_reads_cons(s, i, 0, threadId);
                    g_levs[threadId][i].base = j;
                    g_levs[threadId][i].depth = g_levs[threadId][i].c[j];
                    g_levs[threadId][i].node = sub_node ;
                    return i;

                }
                solved[i] = 1;
                continue;
            }
            if(is_solved != 3 && is_solved != 4)
                *has_repeat_break = 1; 

if(DEBUG)   fprintf(logs, "NotHete2\n");
            //solving repeat by Pair end information.
            if(Solve_Conf == 1 && i == 1 )//&& !(g_levs[threadId][i].c[j] < Error_depth && g_levs[threadId][i].depth > 2 * Error_depth && g_levs[threadId][i].depth + g_levs[threadId][i].c[j] < 4 * Error_depth && diff >= 5))
            {
                is_solved = solve_conf(logs, g_nodes[threadId], g_levs[threadId], g_reads[threadId], s, seed_len, ctg_len, threadId, diff);
                if(is_solved == 1)
                {
if(DEBUG)   fprintf(logs, "conf_solved.\n");
                    solved[i] = 1;
                    return i;
                    //return find_last_lev(logs, s, ctg_seq, seed_seq, seed_len, ctg_len, initial_depth, is_unique, threadId);
                }
            }

if(DEBUG)   fprintf(logs, "CONF fail.\n");
            if(i == 1)
            {
                //fprintf(logs, " %d,%d,%d",seed_len, g_levs[threadId][i].depth, g_levs[threadId][i].c[j]);
                i=-1;
                return i;
            }
            break;
        }
        
        //checking have branch or not.
        if(solved[i] == 0)// || seed_len + i >= depth_6_length || levs[i].depth <= 2 * Error_depth)
        {
            break;
        }
        st[i] = dnaChar[g_levs[threadId][i].base];
    }
    i--;
    return i;
}

int iter_consensus(FILE* logs, int* s, int *len, char* ctg_seq, char* seed_seq, int seed_len, int ctg_len, int keep_used, int initial_depth, int iter, int* is_unique, int threadId)
{
    int i, j=0, k, keep_used_extension=0, has_used = 0, has_repeat_break = 0;

    i = find_last_lev(logs, s, ctg_seq, seed_seq, seed_len, ctg_len, initial_depth, is_unique, &has_repeat_break, threadId);
    if( i <= 0) //lowdepth2.
    {
       s[0]=0; s[1]=0; s[2]=0; *len = 0;
       if(i<0)
           return 5;
       return 3;
    }
    //s is updated here. Used reads.
    i = mask_last_reads(logs, i, s, ctg_seq, keep_used, &keep_used_extension, &has_used, ctg_len, iter, *is_unique, threadId, has_repeat_break, seed_len);

    if(i <=0 )
    {
       *len = 0;
       if(i == -100)
           return 10;

       if(i == 0 && keep_used_extension == 0)
           return 4;
       return 2;
    }

    i = g_levs[threadId][i].node; 
    j = 0;

    while(1)
    {
        for(k=g_nodes[threadId][i].len-1; k>=0; k--, j++)
            g_seq[threadId][j] = g_nodes[threadId][i].seq[k];
        
        if(i==0 && i == g_nodes[threadId][i].parent)
            break;
        i = g_nodes[threadId][i].parent;
    }
    *len = j;

    if(j==0)
    {
        fprintf(stderr, "[%d,%d]",i, g_nodes[threadId][i].len);
    }

    g_seq[threadId][j]='\0';

    if(keep_used_extension == 1)
        return 7;
    if(has_used == 1)
        return 4;

    return 1;
}

/*
    Update reads used in this contigs.
*/

ULL get_pair_id(ULL id)
{
    if(id%2 == 0)
        return id+1;
    return id-1;
}


int is_unique(int depth, int len)
{
	return depth > double(Read_length-len+1)/double(Read_length)*double(Expect_depth)*Upper_bound ? 0 : 1;
}

int is_repeat(int depth, int len)
{
    return depth > (3- Upper_bound)*double(Read_length-len+1)/double(Read_length)*double(Expect_depth) ? 1 : 0;
}

int is_long_repeat_seed(int depth, int len)
{
    return (len >= depth_9_length && depth >= 2 * double(Read_length-len+1)/double(Read_length)*double(Expect_depth));
}

void free_thread_tmpRead(int threadId, int flag)
{
    int i;
    if(flag == 1)
    {
        for(i=0; i< tmpCount[threadId]; i++)
        {
            if(tmpRead[threadId][i].keep == 0)
                continue;
            clean_use_bwt_sparse(tmpRead[threadId][i].id, tmpRead[threadId][i].size);
            endInf[threadId].ur--;
            tmpRead[threadId][i].keep = 0;
        }
    }else{
        for(i=0; i< tmpCount[threadId]; i++)
        {
            if(tmpRead[threadId][i].keep == 0)
                continue;
            int flag = clean_current_bwt_sparse(tmpRead[threadId][i].id, tmpRead[threadId][i].size);
        }
    }
    tmpCount[threadId] = 0;
}

/*
    Update contig and get the new iteration seed using the extended sequence.
*/

void str_copy(char* seq, char* from, int start, int len)
{
	int i, j, l=0;
	for(i=start,j=0,l=0; i< strlen(from) && l<len; i++,j++,l++)
	{
		seq[j]=from[i];
	}
	seq[j]='\0';
}

void reverse_com(char* from, char* to, int len)
{
    int i,j;
    for(i=0; i<len ; i++)
    {
        to[i] = complementMap[ from[len-i-1]];
    }
    to[i]='\0';
}

void clean_trimed_reads(int ctg_len, int read_num, int threadId)
{
    int i=0, j=0;
    for(i = tmpCount[threadId] -1; i>=0 && j < read_num; i--, j++)
    {
        if(tmpRead[threadId][i].pos < ctg_len)
            break;
        if(tmpRead[threadId][i].keep == 0)
            continue;

        clean_use_bwt_sparse(tmpRead[threadId][i].id, tmpRead[threadId][i].size);
        endInf[threadId].ur--;
        tmpRead[threadId][i].keep = 0;
    }
}


int update_ctg(int len, int cur_iter, int read_num, int threadId)
{
    //update the sequence of ctg.

    if(len + ctgs[threadId][ctg_num[threadId]].len >= MaxCtgLen)
    {
        fprintf(stderr, "Too long contigs with length %d, %d\n", ctgs[threadId][ctg_num[threadId]].len, len);
        return len;
    }
    int i=ctgs[threadId][ctg_num[threadId]].len, j=0;
    int oldctglen = i;
    for(; j < len ; i++, j++)
    {
        ctgs[threadId][ctg_num[threadId]].seq[i] = g_seq[threadId][len - 1 - j];
    }
    ctgs[threadId][ctg_num[threadId]].seq[i]='\0';
    ctgs[threadId][ctg_num[threadId]].len = i;

    //get the new seed.
    if(cur_iter < 1)
    {
        fprintf(stderr, "Please check: %d\n", cur_iter);
        exit(2);
    }
    char tmpseq[len + g_iters[threadId][cur_iter-1].slen];
    char rctmpseq[len + g_iters[threadId][cur_iter-1].slen];

    i=0; 
    for(j=0; j<len; i++, j++)
        tmpseq[i] = g_seq[threadId][j];

    for(j=0; j< g_iters[threadId][cur_iter-1].slen; i++, j++)
        tmpseq[i] = g_iters[threadId][cur_iter-1].seq[j];

    tmpseq[i]='\0';
    int tmplen = i;
    reverse_com(tmpseq, rctmpseq, tmplen);
    int end = tmplen - 1;
    for(i=0; i< bwtCount; i++)
    {
        tmpbase[threadId][i].saL = 0; tmpbase[threadId][i].saR = 0; 
        tmpbase[threadId][i].plus = 0; tmpbase[threadId][i].minus = 0; 
        tmpbase[threadId][i].saLL = 0; tmpbase[threadId][i].saRR = 0;
    }
    unsigned long long depth=0, tdepth;

    int prev_repeat = is_repeat(g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
    int trim_len = 0, rep_start = -1;

if(DEBUG)    fprintf(stdout, "[iter] %d\n", cur_iter-1);
    int start = extend_backwards_rep(rctmpseq, 0, end, tmpbase[threadId], &(depth), &rep_start);
if(DEBUG)    fprintf(stdout, "[iter] %d,d=%d,l=%d, r=%d, start=%d, depth=%d\n", cur_iter-1, g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen, prev_repeat, start, depth);

    if((start < 0 && depth > Error_depth) || prev_repeat ==1)
    {

        if(rep_start != -1)
        {
            end = tmplen - rep_start -1;
        }else
            end = strlen(tmpseq) - 1;
        
        if(end >= strlen(tmpseq))
	    fprintf(stdout, "T1: %s, %d; new seq %s; iter seq %s; newlen %d, seedlen %d, end %d; iter %d\n", tmpseq, strlen(tmpseq), g_seq[threadId], g_iters[threadId][cur_iter-1].seq, len, g_iters[threadId][cur_iter-1].slen, end, cur_iter-1);
        set_SAranges(tmpseq, 0, end, tmpbase[threadId], &(depth));
    }else{

        if(prev_repeat == 1)
        {
            int break_start;
            fprintf(stdout, "Iter %d PrevRepeat: %d, %d; ", cur_iter, g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
            while(start >= 0 && end > g_iters[threadId][cur_iter-1].slen)
            {
                end --;
                start = extend_backwards(rctmpseq, 0, end, tmpbase[threadId], &(depth), &break_start);
            }
            //fprintf(stdout, "Findstart %d ", start);

            if(start <0)
            {
                trim_len = tmplen - 1 - end;
                end = MINOVERLAP - 1 + trim_len;
                set_SAranges(tmpseq, trim_len, end, tmpbase[threadId], &(depth));
                while(is_repeat(depth, MINOVERLAP) == 0 && end < tmplen -1 && trim_len < len)
                {
                    trim_len ++;
                    end ++;
                    set_SAranges(tmpseq, trim_len, end, tmpbase[threadId], &(depth));
                }
                if(end == tmplen -1 || trim_len == len)
                {
                    fprintf(stdout, "all trimed0.\n");
                    ctgs[threadId][ctg_num[threadId]].len = oldctglen;
                    ctgs[threadId][ctg_num[threadId]].seq[ ctgs[threadId][ctg_num[threadId]].len ] = '\0';
                    return len;
                }
                fprintf(stdout, "FinalTrim%d and %d\n", trim_len, end);
                ctgs[threadId][ctg_num[threadId]].len -= trim_len;
                ctgs[threadId][ctg_num[threadId]].seq[ ctgs[threadId][ctg_num[threadId]].len ] = '\0';
            }else{
                fprintf(stdout, "all trimed1.\n");
                ctgs[threadId][ctg_num[threadId]].len = oldctglen;
                ctgs[threadId][ctg_num[threadId]].seq[ ctgs[threadId][ctg_num[threadId]].len ] = '\0';
                return len;
            }
        }else{
            //find a unique seed for the new iteration.
            if(depth <= Error_depth && g_iters[threadId][cur_iter-1].depth > Error_depth + 1)
            {
                trim_len = 1; int start_trim=0, now_rep_start = -1;
if(DEBUG)       fprintf(stdout, "Iter %d PrevRepeat: %d, %d; ", cur_iter, g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
                while(depth <= Error_depth && trim_len < len )
                {
                    now_rep_start = -1;
if(DEBUG)           fprintf(stdout, "Ending at %d\n", end-trim_len);
                    start_trim = extend_backwards_rep(rctmpseq, 0, end-trim_len, tmpbase[threadId], &(depth), &now_rep_start);
if(DEBUG)           fprintf(stdout, "End got depth %d, start %d\n", depth, start_trim);
                    trim_len++;
                }

                if(trim_len <= len && depth > Error_depth)
                {
                    trim_len--;
if(DEBUG)           fprintf(stdout, "LowTrim%d \n", trim_len);
                    ctgs[threadId][ctg_num[threadId]].len -= trim_len;
                    ctgs[threadId][ctg_num[threadId]].seq[ ctgs[threadId][ctg_num[threadId]].len ] = '\0';
                    if(start_trim >= 0)
                        end = tmplen -1 - start_trim;
                    else{
                        if(now_rep_start != -1)
                            end = tmplen -1 -  now_rep_start;
                        else
                            end = strlen(tmpseq) -1;
                    }
                    clean_trimed_reads(ctgs[threadId][ctg_num[threadId]].len, read_num, threadId);
                }else{
if(DEBUG)           fprintf(stdout, "Lowtrim fail with trim_len %d and start_trim %d, final depth %lld\n", trim_len, start_trim, depth);

                    trim_len = 0; end = tmplen - 1 - start;
                    if(start < 0)
                    {
                        if(rep_start != -1)
                            end = tmplen -1 - rep_start;
                        else
                            end = strlen(tmpseq) -1;
                    }
                }
            }else{
                trim_len = 0; end = tmplen - 1 - start;
                if(start < 0)
                {
                    if(rep_start != -1)
                    {
                        end = tmplen - rep_start -1;
//                        fprintf(stdout, "New seed length2: %d for iter %d\n", end+1, cur_iter);
                    }else
                        end = strlen(tmpseq) -1;
                }
            }
            if(end >= strlen(tmpseq))
                fprintf(stdout, "T2: %s, %d %d; iter %d\n", tmpseq, strlen(tmpseq), end, cur_iter-1);
            set_SAranges(tmpseq, trim_len, end, tmpbase[threadId], &(depth));
        }
    }
    
    set_base(g_iters[threadId][cur_iter].b, tmpbase[threadId]);

    if(depth > MAXDEPTH)
    	g_iters[threadId][cur_iter].depth = MAXDEPTH;
    else
        g_iters[threadId][cur_iter].depth = depth;

    g_iters[threadId][cur_iter].slen = end - trim_len + 1;
    str_copy(g_iters[threadId][cur_iter].seq, tmpseq, trim_len, g_iters[threadId][cur_iter].slen);

    g_iters[threadId][cur_iter].rnum = 0;
    g_iters[threadId][cur_iter].end = 0;
    g_iters[threadId][cur_iter].pe = 0;
    g_iters[threadId][cur_iter].len = 0;

    return trim_len;
}

/*
    Using PE information to improve iteration.
    iter_num >=0, read_num >=0.
    
    Output: 
    Fill in ctg, tmpRead, iters.
*/
void print_inf(FILE* fp, Node* nodes, Lev* levs, ReadId* reads, int* s, int cur_len)
{
    fprintf(fp, "current contig length: %d\n", cur_len);
    fprintf(fp, "Node:\n");
    int i;
    for(i=0; i< s[0]; i++)
        fprintf(fp, "%d\t%s,%d,%d,%d,%d,%d\n",i, nodes[i].seq, nodes[i].len, nodes[i].level, nodes[i].depth, nodes[i].parent, nodes[i].is_end);

    fprintf(fp, "Read:\n");
    for(i=0; i< s[2]; i++)
        fprintf(fp, "%d\t%d,%d,%d,%d,%d\n",i, reads[i].id, reads[i].node, reads[i].level, reads[i].cons, reads[i].strand);
    
    fprintf(fp, "Lev:\n");
    for(i=1; i< s[1]; i++)
        fprintf(fp, "%d\t%d,%d,%d,%d\n",i, levs[i].base, levs[i].depth, levs[i].node, levs[i].nnum);
   
}

void mask_branch_reads( Node* nodes, Lev* levs, ReadId* reads, int* s)
{
    int i;
    int node_type[s[0]];
    for(i=0; i< s[0]; i++)
        node_type[i]=0;

    for(i=1; i< s[1]; i++)
    {
        node_type[ levs[i].node ]=1;
    }

    for(i=0; i< s[2]; i++)
    {
        if(node_type[ reads[i].node ] == 0)
	{
            reads[i].node=0; reads[i].level=0;
        }
    }
}

int backward_extending_bwtpe(FILE* logs, int is_first, int threadId)
{
    int cur_iter = iterCount[threadId] + 1;
    //fprintf(stderr, "thread%d and cur %d\n", threadId, *iter_num);
    if(cur_iter > MaxIterNum-2)
    {
if(DEBUG)        fprintf(logs, " TooManyIteration");
        g_iters[threadId][cur_iter-1].end = 4; iterCount[threadId] = cur_iter;
        return 0;
    }

    if(g_iters[threadId][cur_iter-1].depth < Error_depth)
    {
if(DEBUG)        fprintf(logs, " LowDepth1");
        g_iters[threadId][cur_iter-1].end = 0; iterCount[threadId] = cur_iter;
        return 0;
    }

    if(g_iters[threadId][cur_iter-1].depth >= MAXDEPTH)
    {
if(DEBUG)        fprintf(logs, " HighDepth%d,%d", g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
        g_iters[threadId][cur_iter-1].end = 2; iterCount[threadId] = cur_iter;
	return 0;
    }

    if((cur_iter > 1 && strcmp(g_iters[threadId][cur_iter-1].seq, g_iters[threadId][cur_iter-2].seq) == 0) 
       || (cur_iter >2 && strcmp(g_iters[threadId][cur_iter-1].seq, g_iters[threadId][cur_iter-3].seq) == 0)
       || (cur_iter >3  && strcmp(g_iters[threadId][cur_iter-1].seq, g_iters[threadId][cur_iter-4].seq) == 0)
      )
    {
if(DEBUG)        fprintf(logs, " SameSeed: now seed %s, old1 %s\n", g_iters[threadId][cur_iter-1].seq, g_iters[threadId][cur_iter-2].seq);
        g_iters[threadId][cur_iter-1].end = 5; iterCount[threadId] = cur_iter;
        return 0;
    }

    if(cur_iter > 1)
    {
        double depth1 = double(Read_length*g_iters[threadId][cur_iter-1].depth)/double(Read_length-g_iters[threadId][cur_iter-1].slen+1);
        double depth2 = double(Read_length*g_iters[threadId][cur_iter-2].depth)/double(Read_length-g_iters[threadId][cur_iter-2].slen+1);
        if(depth1 * 6 < depth2 && depth1 < Expect_depth/3)
        {
            if(DEBUG)        fprintf(logs, " LargeDecrease%.1f,%.1f",depth1, depth2);
            g_iters[threadId][cur_iter-1].end = 2; iterCount[threadId] = cur_iter;
            return 0;
        }   
    }
    //define and initialize extending structures.
    int len = 0;
    int max_node_num = g_iters[threadId][ iterCount[threadId] ].depth * Read_length;
    
    int s[]={0,0,0};
    int i,j;
    clock_t start_time;

    for(i=0; i< max_node_num; i++)
    {
        g_nodes[threadId][i].seq[0] = '\0';
        g_nodes[threadId][i].len = 0;
        g_nodes[threadId][i].level = 0;
        g_nodes[threadId][i].parent = 0;
        g_nodes[threadId][i].is_end = 0;
        g_nodes[threadId][i].conf = 0;
        set_blank_base( g_nodes[threadId][i].b);
    }
    set_base(g_nodes[threadId][0].b, g_iters[threadId][ cur_iter - 1].b);

    g_nodes[threadId][0].depth = g_iters[threadId][ cur_iter - 1 ].depth;
    g_seq[threadId][0]='\0';
    s[0]++;
    
    for(i=0; i< Read_length; i++)
    {
        g_levs[threadId][i].nnum = 0;
        for(j=0; j<5; j++)
            g_levs[threadId][i].c[j]=0;
    }
    g_levs[threadId][0].node = 0;
    s[1]++;

    //fprintf(stderr, "Iter %d is fisrt %d with given depth %d from depth %d and calculated %d\n", cur_iter - 1, is_first, get_real_depth(g_nodes[threadId][0].b), get_real_depth(g_iters[threadId][cur_iter - 1].b), g_iters[threadId][cur_iter - 1].depth);
    //starting iterate.
    iterCount[threadId] = cur_iter;
    int is_continue = 1;
    start_time = clock();
    while(is_continue == 1 && s[1] + g_iters[threadId][cur_iter-1].slen <= Read_length + 1)
    {
        is_continue = backward_search(logs, s, threadId);
    }

if(DEBUG)
{
    fprintf(logs, "thread=%d;%d,%d,%d,%d,%d,%d\n",threadId, cur_iter-1, g_iters[threadId][cur_iter-1].slen, g_iters[threadId][cur_iter-1].depth, s[0], s[1], s[2]);
   
    if(cur_iter < 7000)
    {
    	fprintf(logs, "\nIter %d, depth %.1f\n", cur_iter-1, double(Read_length*g_iters[threadId][cur_iter-1].depth)/double(Read_length-g_iters[threadId][cur_iter-1].slen+1)); 
        if(s[1] > 0)
            print_inf(logs, g_nodes[threadId], g_levs[threadId], g_reads[threadId], s, ctgs[threadId][ctg_num[threadId]].len);
    }
}   
    //fflush(stdout);
   
    //local consensus.
    int extension_type =1;
    int unique = is_unique(g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
    int repeat = is_repeat(g_iters[threadId][cur_iter-1].depth, g_iters[threadId][cur_iter-1].slen);
    int initial_read_num=s[2], is_now_unique = unique;

    extension_type = iter_consensus(logs, s, &len, ctgs[threadId][ctg_num[threadId]].seq, g_iters[threadId][cur_iter-1].seq, g_iters[threadId][cur_iter-1].slen, ctgs[threadId][ctg_num[threadId]].len, 0, g_iters[threadId][cur_iter-1].depth, cur_iter-1, &is_now_unique, threadId );

    if(extension_type == 10)
    {
        //meet thread conflition.
        return 1;
    }

if(DEBUG)    fprintf(logs, "thread%d Seq: %s, len %d, %d, %d, %d\n", threadId, g_seq[threadId], len, s[0], s[1], s[2]); 

    int trim_len = 0;
    if(ctgs[threadId][ctg_num[threadId]].len > Read_length && extension_type == 4 && (is_now_unique == 1 || repeat == 0))// && ctgs[threadId][ctg_num[threadId]].len >= 2* Read_length)
    {
        ctgs[threadId][ctg_num[threadId]].len = ctgs[threadId][ctg_num[threadId]].len - Read_length;
        ctgs[threadId][ctg_num[threadId]].seq[ctgs[threadId][ctg_num[threadId]].len] = '\0';
        clean_trimed_reads(ctgs[threadId][ctg_num[threadId]].len, tmpCount[threadId], threadId);
        trim_len = len;
        extension_type = 2;
    }else if(len > 0)
    {
        trim_len = update_ctg(len, cur_iter, s[2], threadId);
        if(len == trim_len)
        {
            s[2]=0;
            len = 0;
            trim_len = 0;
            extension_type = 6;
        }
    }else if(extension_type == 1)
        extension_type = 0;

    //update ctg and iters.
    g_iters[threadId][cur_iter-1].len = len - trim_len;
    g_iters[threadId][cur_iter-1].rnum = s[2];
    g_iters[threadId][cur_iter-1].pe = s[0]>0;
    g_iters[threadId][cur_iter-1].end = 0;

    //further extending.
    switch (extension_type)
    {
        case 0:
            //fprintf(logs, " NOEXTEND");
            g_iters[threadId][cur_iter-1].end = 4;
            break;
        case 1:
            //fprintf(stderr, "SIN");
            return backward_extending_bwtpe(logs, 0, threadId);
            break;
        case 2:
if(DEBUG)   fprintf(logs, "UsedRead");
            g_iters[threadId][cur_iter-1].end = 3;
            //backward_extending_bwtpe(ctg, tmpRead, iters, iter_num, read_num);
            break;
        case 3:
if(DEBUG)   fprintf(logs, "LowDepth2");
            g_iters[threadId][cur_iter-1].end = 1;
            break;
        case 4:
if(DEBUG)   fprintf(logs, "UsedRead2");
            g_iters[threadId][cur_iter-1].end = 6;
            break;
        case 5:
if(DEBUG)   fprintf(logs, "CB");
            g_iters[threadId][cur_iter-1].end = 7;
            break;
        case 6:
if(DEBUG)   fprintf(logs, " R2U");
            g_iters[threadId][cur_iter-1].end = 8;
            break;
        case 7:
if(DEBUG)   fprintf(logs, "CurrentUse");
            g_iters[threadId][cur_iter-1].end = 9;
            break;
        case 8:
if(DEBUG)   fprintf(logs, "Cut");
            g_iters[threadId][cur_iter-1].end = 10;
            break;
        case 9:
if(DEBUG)   fprintf(logs, "CurrentUse0");
            g_iters[threadId][cur_iter-1].end = 11;
            break;
        default:
            break;
    }
    //fflush(logs);
    return 0;
}


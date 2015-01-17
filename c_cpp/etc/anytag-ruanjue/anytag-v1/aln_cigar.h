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
 
#ifndef __ALN_CIGAR_RJ_H
#define __ALN_CIGAR_RJ_H

#include <stdlib.h>
#include <stdio.h>

#define ALN_CIGAR_MAX_LEN	8191

#define ALN_CIGAR_TYPE_NULL	0
#define ALN_CIGAR_TYPE_MAT	3
#define ALN_CIGAR_TYPE_INS	1
#define ALN_CIGAR_TYPE_DEL	2
#define ALN_CIGAR_TYPE_SKIP	7
#define ALN_CIGAR_TYPE_CLIP1	5
#define ALN_CIGAR_TYPE_CLIP2	6

static const char aln_cigar_string[8] = "?IDM?SHN";

typedef struct {
	uint16_t len:13, type:3;
} AlnCigar;

static inline void cigars_lengths(AlnCigar *cigars, int n_cigar, int *aln_size, int *seq1_size, int *seq2_size){
	int i;
	if(aln_size) *aln_size = 0;
	if(seq1_size) *seq1_size = 0;
	if(seq2_size) *seq2_size = 0;
	for(i=0;i<n_cigar;i++){
		if(seq1_size &&  (cigars[i].type & 0x01)) *seq1_size += cigars[i].len;
		if(seq2_size &&  (cigars[i].type & 0x02)) *seq2_size += cigars[i].len;
		if(aln_size  && !(cigars[i].type & 0x04)) *aln_size += cigars[i].len;
	}
}

static inline int _aln_cigar_h_num_str_len(int n){
	int i;
	i = 0;
	while(n){
		i ++;
		n /= 10;
	}
	return i;
}

static inline int _aln_cigar_add_cigar(AlnCigar *cs, int n_cigar, int len, int type){
	while(len){
		if(len > ALN_CIGAR_MAX_LEN){
			cs[n_cigar].len = ALN_CIGAR_MAX_LEN;
			len -= ALN_CIGAR_MAX_LEN;
		} else {
			cs[n_cigar].len = len;
			len = 0;
		}
		cs[n_cigar++].type = type;
	}
	return n_cigar;
}

static inline char* cigars2string(AlnCigar *cigars, int n_cigar, char *str){
	int i, j, n, str_len, type;
	char *p;
	if(str == NULL){
		str_len = 0;
		for(i=0;i<n_cigar;i++) str_len += _aln_cigar_h_num_str_len(cigars[i].len) + 1;
		str = malloc(str_len + 1);
	}
	p = str;
	if(n_cigar){
		n = cigars[0].len;
		type = cigars[0].type;
		for(i=1;i<=n_cigar;i++){
			if(i == n_cigar || (type != cigars[i].type && n)){
				str_len = _aln_cigar_h_num_str_len(n) - 1;
				j = 0;
				while(n){
					p[str_len - j] = '0' + (n % 10);
					n /= 10;
					j ++;
				}
				p[str_len + 1] = aln_cigar_string[type];
				p = p + str_len + 1 + 1;
			}
			if(i == n_cigar) break;
			n += cigars[i].len;
			type = cigars[i].type;
		}
	}
	p[0] = 0;
	return str;
}

static inline int string2cigars(AlnCigar *cigars, char *str, int len){
	int i, n, x;
	n = 0;
	x = 0;
	for(i=0;i<len;i++){
		switch(str[i]){
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9': x = x * 10 + (str[i] - '0'); break;
			case 'M':
			case 'm': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_MAT); x = 0; break;
			case 'I':
			case 'i': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_INS); x = 0; break;
			case 'D':
			case 'd': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_DEL); x = 0; break;
			case 'S':
			case 's': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_CLIP1); x = 0; break;
			case 'H':
			case 'h': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_CLIP2); x = 0; break;
			case 'N':
			case 'n': n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_SKIP); x = 0; break;
			default : n = _aln_cigar_add_cigar(cigars, n, x, ALN_CIGAR_TYPE_NULL); x = 0; break;
		}
	}
	return n;
}

static inline int rank_cigars_seqlen(AlnCigar *cigars, int n_cigar, int len, int seq_idx){
	int i, ret;
	if(len < 0) len = 0x7FFFFFFF;
	for(i=0,ret=0;i<n_cigar&&len>=0;i++){
		if((cigars[i].type >> seq_idx) & 0x01){
			if(len > (int)cigars[i].len){
				ret += cigars[i].len;
			} else {
				ret += len;
			}
		}
		len -= cigars[i].len;
	}
	return ret;
}

static inline int rev_rank_cigars_seqlen(AlnCigar *cigars, int n_cigar, int len, int seq_idx){
	int i, ret;
	if(len < 0) len = 0x7FFFFFFF;
	for(i=n_cigar-1,ret=0;i>=0&&len>=0;i--){
		if((cigars[i].type >> seq_idx) & 0x01){
			if(len > (int)cigars[i].len){
				ret += cigars[i].len;
			} else {
				ret += len;
			}
		}
		len -= cigars[i].len;
	}
	return ret;
}

static inline int select_cigars_seqlen(AlnCigar *cigars, int n_cigar, int len, int seq_idx){
	int i, ret;
	for(i=0,ret=0;i<n_cigar;i++){
		if((cigars[i].type >> seq_idx) & 0x01){
			if(len > (int)cigars[i].len){
				ret += cigars[i].len;
			} else {
				ret += len;
				break;
			}
			len -= cigars[i].len;
		} else {
			ret += cigars[i].len;
		}
	}
	return ret;
}

static inline int rev_select_cigars_seqlen(AlnCigar *cigars, int n_cigar, int len, int seq_idx){
	int i, ret;
	for(i=n_cigar-1,ret=0;i>=0;i--){
		if((cigars[i].type >> seq_idx) & 0x01){
			if(len > (int)cigars[i].len){
				ret += cigars[i].len;
			} else {
				ret += len;
				break;
			}
			len -= cigars[i].len;
		} else {
			ret += cigars[i].len;
		}
	}
	return ret;
}

static inline void flip_cigars(AlnCigar *cigars, int n_cigar){
	int i;
	for(i=0;i<n_cigar;i++){
		cigars[i].type = (cigars[i].type & 0x04) | ((cigars[i].type & 0x01) << 1) | ((cigars[i].type & 0x02) >> 1);
	}
}

static inline int sub_cigars(AlnCigar *dst, AlnCigar *cigars, int n_cigar, int off, int len){
	int i, x, y, n_sub;
	n_sub = 0;
	if(len < 0) len = 0x3FFFFFFF;
	else if(len == 0) return 0;
	for(i=0,x=0;i<n_cigar;i++){
		y = x + cigars[i].len;
		if(x < off){
			if(y < off){
			} else if(y < off + len){
				n_sub = _aln_cigar_add_cigar(dst, n_sub, y - off, cigars[i].type);
			} else {
				n_sub = _aln_cigar_add_cigar(dst, n_sub, off + len - x, cigars[i].type);
				break;
			}
		} else if(x >= off + len){
			break;
		} else {
			if(y < off + len){
				n_sub = _aln_cigar_add_cigar(dst, n_sub, cigars[i].len, cigars[i].type);
			} else {
				n_sub = _aln_cigar_add_cigar(dst, n_sub, off + len - x, cigars[i].type);
				break;
			}
		}
		x = y;
	}
	return n_sub;
}

static inline int sub_seq_cigars(AlnCigar *dst, AlnCigar *c, int n, int seq_idx, int off, int len){
	int d, i, x, y;
	if(len < 0) len = 0x7FFFFFFF;
	d = 0;
	x = y = 0;
	for(i=0;i<n;i++){
		if(x > off + len) break;
		if((c[i].type >> seq_idx) & 0x01){
			y = x + c[i].len;
		} else {
			y = x;
		}
		if(y < off) continue;
		if(x < off){
			if(y < off + len){
				d = _aln_cigar_add_cigar(dst, d, y - off, c[i].type);
			} else {
				d = _aln_cigar_add_cigar(dst, d, off + len - x, c[i].type);
				break;
			}
		} else if(y < off + len){
			d = _aln_cigar_add_cigar(dst, d, c[i].len, c[i].type);
		} else {
			d = _aln_cigar_add_cigar(dst, d, off + len - x, c[i].type);
			break;
		}
		x = y;
	}
	return d;
}

static inline int cat_cigars(AlnCigar *cigars1, int n_cigar1, AlnCigar *cigars2, int n_cigar2){
	int i;
	if(n_cigar2 == 0) return n_cigar1;
	if(n_cigar1){
		if(cigars1[n_cigar1-1].type == cigars2[0].type){
			n_cigar1 = _aln_cigar_add_cigar(cigars1, n_cigar1 - 1, cigars1[n_cigar1 - 1].len + cigars2[0].len, cigars2[0].type);
			i = 1;
		} else i = 0;
	} else i = 0;
	while(i < n_cigar2){
		cigars1[n_cigar1].type = cigars2[i].type;
		cigars1[n_cigar1].len = cigars2[i].len;
		n_cigar1 ++;
		i ++;
	}
	return n_cigar1;
}

static inline int append_cigars(AlnCigar *cs, int n, int type, int len){
	if(n && cs[n-1].type == type){
		return _aln_cigar_add_cigar(cs, n-1, cs[n-1].len + len, type);
	} else {
		return _aln_cigar_add_cigar(cs, n, len, type);
	}
}

static inline void reverse_cigars(AlnCigar *cs, int n){
	int i, j;
	AlnCigar c;
	i = 0;
	j = n - 1;
	while(i < j){
		c = cs[i]; cs[i] = cs[j]; cs[j] = c;
		i ++; j --;
	}
}

// seq1_size of c1 must be no less than c2
// seq1_size of dst will equal length of c1
static inline int compile_cigars(AlnCigar *dst, AlnCigar *c1, int n1, AlnCigar *c2, int n2, int seq_idx){
	int i, j, n3, x1, x2, x3, f1, f2, f3;
	n3 = 0;
	x2 = 0;
	x3 = 0;
	f2 = ALN_CIGAR_TYPE_MAT;
	f3 = ALN_CIGAR_TYPE_INS;
	for(i=j=0;i<n1;i++){
		x1 = c1[i].len;
		f1 = (c1[i].type >> seq_idx) & 0x01;
		while(x1){
			if(x2 == 0){
				if(j < n2){
					x2 = c2[j].len;
					f2 = (c2[j].type >> seq_idx) & 0x01;
					j ++;
				} else {
					x2 = x1;
					f2 = f1;
				}
			}
			if(f1){
				if(f2){
					if(x3 && f3 != ALN_CIGAR_TYPE_INS){ n3 = _aln_cigar_add_cigar(dst, n3, x3, f3); x3 = 0; }
					if(x1 < x2){ x2 -= x1; x3 += x1; x1 = 0; }
					else { x1 -= x2; x3 += x2; x2 = 0; }
					f3 = ALN_CIGAR_TYPE_INS;
				} else {
					if(x3 && f3 != ALN_CIGAR_TYPE_DEL){ n3 = _aln_cigar_add_cigar(dst, n3, x3, f3); x3 = 0; }
					f3 = ALN_CIGAR_TYPE_DEL;
					x3 += x2; x2 = 0;
				}
			} else {
				if(f2){
					if(x3 && f3 != ALN_CIGAR_TYPE_INS){ n3 = _aln_cigar_add_cigar(dst, n3, x3, f3); x3 = 0; }
					f3 = ALN_CIGAR_TYPE_INS;
					x3 += x1; x1 = 0;
				} else {
					if(x3 && f3 != ALN_CIGAR_TYPE_INS){ n3 = _aln_cigar_add_cigar(dst, n3, x3, f3); x3 = 0; }
					f3 = ALN_CIGAR_TYPE_INS;
					if(x1 < x2){ x2 -= x1; x3 += x1; x1 = 0; }
					else { x1 -= x2; x3 += x2; x2 = 0; }
				}
			}
		}
	}
	if(x3){ n3 = _aln_cigar_add_cigar(dst, n3, x3, f3); }
	return n3;
}

// c2 is the dst in align_cigars
static inline int apply_cigars(AlnCigar *dst, AlnCigar *c1, int n1, AlnCigar *c2, int n2){
	int i, j, n3, x1, x2, x3, f1, f2;
	n3 = 0;
	x2 = 0;
	f2 = ALN_CIGAR_TYPE_MAT;
	for(i=j=0;i<n1;i++){
		x1 = c1[i].len;
		f1 = c1[i].type;
		while(x1){
			if(x2 == 0){
				if(j < n2){
					x2 = c2[j].len;
					f2 = c2[j].type;
					j ++;
				} else {
					x2 = x1;
					f2 = ALN_CIGAR_TYPE_MAT;
				}
			}
			if(f2 == ALN_CIGAR_TYPE_DEL){
				n3 = _aln_cigar_add_cigar(dst, n3, x2, ALN_CIGAR_TYPE_DEL);
				x2 = 0;
			} else {
				x3 = (x1 < x2)? x1 : x2;
				n3 = _aln_cigar_add_cigar(dst, n3, x3, f1);
				x1 -= x3;
				x2 -= x3;
			}
		}
	}
	return n3;
}

static inline int refine_cigars(AlnCigar *c, int n){
	int i, j, x, f;
	x = 0;
	f = ALN_CIGAR_TYPE_NULL;
	for(i=j=0;i<n;i++){
		if(c[i].type == f){
			x += c[i].len;
		} else {
			if(x){
				j = _aln_cigar_add_cigar(c, j, x, f);
			}
			x = c[i].len;
			f = c[i].type;
		}
	}
	if(x){ j = _aln_cigar_add_cigar(c, j, x, f); }
	return j;
}

static inline int cigars_seq2aln(char *dst, AlnCigar *c, int n, int seq_idx, char *seq){
	int i, j, k, m;
	for(i=m=k=0;i<n;i++){
		if((c[i].type >> seq_idx) & 0x01){
			if(c[i].type & 0x04){
				k += c[i].len;
			} else {
				for(j=0;j<c[i].len;j++) dst[m++] = seq[k++];
			}
		} else if(!(c[i].type & 0x04)){
			for(j=0;j<c[i].len;j++) dst[m++] = '-';
		}
	}
	dst[m] = '\0';
	return m;
}

#endif

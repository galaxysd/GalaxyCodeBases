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
 
#ifndef __SORT_RJ_H
#define __SORT_RJ_H

#include <stdio.h>
#include <stdlib.h>

#define num_cmp_script(e1, e2, obj, val_macro) ((val_macro(e1, obj) == val_macro(e2, obj))? 0 : ((val_macro(e1, obj) < val_macro(e2, obj))? -1 : 1))

#define define_bubble_sort(name, e_type, cmp_func)	\
static inline void name(e_type* list, size_t size, void *ref){	\
	size_t i, j, n;	\
	e_type t;	\
	i = 0;	\
	while(i < size){	\
		n = 0;	\
		for(j=size-1;j>i;j--){	\
			if(cmp_func(list[j-1], list[j], ref) == 1){	\
				t = list[j-1]; list[j-1] = list[j]; list[j] = t;	\
				n = 1;	\
			}	\
		}	\
		if(n == 0) break;	\
		i ++;	\
	}	\
	if(ref == ref) return;	\
}

#define define_quick_sort(name, e_type, cmp_func)	\
static inline void name##_core(e_type* list, size_t s, size_t e, void *ref){	\
	e_type pivot, tmp;	\
	size_t i, j, m, ok;	\
	int cmp;	\
	if(s + 1 >= e) return;	\
	i = s;	\
	if(cmp_func(list[i], list[e - 1], ref) == -1) i = e - 1;	\
	if(cmp_func(list[i], list[(s + e)/2], ref) == -1) i = (s + e) / 2;	\
	if(i != e - 1){ pivot = list[i]; list[i] = list[e - 1]; list[e - 1] = pivot; }	\
	else pivot = list[i];	\
	j = s;	\
	ok = 0;	\
	for(m=s;m<e-1;m++){	\
		cmp = cmp_func(list[m], pivot, ref);	\
		if(cmp) ok = 1;	\
		if(cmp != 1){	\
			if(m != j){ tmp = list[m]; list[m] = list[j]; list[j] = tmp; }	\
			j ++;	\
		}	\
	}	\
	if(j != e - 1){ tmp = list[j]; list[j] = list[e - 1]; list[e - 1] = tmp; }	\
	if(ok == 0) return;	\
	name##_core(list, s, j, ref);	\
	name##_core(list, j + 1, e, ref);	\
}	\
static inline void name(e_type* list, size_t size, void *ref){ name##_core(list, 0, size, ref); }

#define define_merge(name, e_type, cmp_func, output_func)	\
static inline void name(e_type *list1, size_t size1, e_type *list2, size_t size2, void *ref){	\
	size_t i, j;	\
	i = j = 0;	\
	while(i < size1 && j < size2){	\
		if(cmp_func(list1[i], list2[j], ref) != 1){	\
			output_func(list1[i], ref);	\
			i ++;	\
		} else {	\
			output_func(list2[j], ref);	\
			j ++;	\
		}	\
	}	\
	while(i < size1){ output_func(list1[i++], ref); }	\
	while(j < size2){ output_func(list2[j++], ref); }	\
}	\
	\
static inline size_t name##_files(FILE **files, int n, void *ref){	\
	e_type *es;	\
	int *flags, i, min;	\
	size_t ret;	\
	ret = 0;	\
	es = malloc(sizeof(e_type) * n);	\
	flags = malloc(sizeof(int) * n);	\
	for(i=0;i<n;i++) flags[i] = 0;	\
	while(1){	\
		min = -1;	\
		for(i=0;i<n;i++){	\
			if(flags[i] == 0){	\
				flags[i] = (fread(es + i, sizeof(e_type), 1, files[i]) == 1)? 1 : 2;	\
			}	\
			if(flags[i] == 1){	\
				if(min == -1){	\
					min = i;	\
				} else if(cmp_func(es[i], es[min], ref) != 1){	\
					min = i;	\
				}	\
			}	\
		}	\
		if(min == -1) break;	\
		output_func(es[min], ref);	\
		flags[min] = 0;	\
		ret ++;	\
	}	\
	free(es);	\
	free(flags);	\
	return ret;	\
}

#define define_unique_merge(name, e_type, cmp_func, output_func)	\
static inline void name(e_type *list1, size_t size1, e_type *list2, size_t size2, void *ref){	\
	size_t i, j;	\
	i = j = 0;	\
	while(i < size1 && j < size2){	\
		switch(cmp_func(list1[i], list2[j])){	\
			case 0:  output_func(list1[i++], ref); j ++; break;	\
			case -1: output_func(list1[i++], ref); break;	\
			default: output_func(list2[j++], ref);	\
		}	\
	}	\
	while(i < size1){ output_func(list1[i++], ref); }	\
	while(j < size2){ output_func(list2[j++], ref); }	\
}

#define define_reverse_array(name, e_type)	\
static inline void name(e_type *list, size_t size){	\
	size_t i, j;	\
	e_type t;	\
	if(size == 0) return;	\
	i = 0;	\
	j = size - 1;	\
	while(i < j){	\
		t = list[i]; list[i] = list[j]; list[j] = t;	\
		i ++; j --;	\
	}	\
}

#define define_apply_array(name, e_type, apply_func)	\
static inline size_t name(e_type *list, size_t size, void *ref){	\
	size_t i, ret;	\
	ret = 0;	\
	for(i=0;i<size;i++){	\
		ret += apply_func(list[i], ref);	\
	}	\
	return ret;	\
	ref = NULL;	\
}

#define define_search_array(name, e_type, cmp_func)	\
static inline long long name(e_type *array, long long size, e_type key, void *ref){	\
	long long i, j, m;	\
	i = 0;	\
	j = size;	\
	while(i < j){	\
		m = i + (j - i) / 2;	\
		if(cmp_func(array[m], key, ref) < 0){	\
			i = m + 1;	\
		} else {	\
			j = m;	\
		}	\
	}	\
	if(i < size && cmp_func(array[i], key, ref) == 0) return i;	\
	else return - (i + 1);	\
	if(0 && ref) return 0;	\
}

#endif

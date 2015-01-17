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
 
#ifndef __LIST_RJ_H
#define __LIST_RJ_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

/**
 * * Common sttaic functions
 * */

#define num_min(n1, n2) (((n1) < (n2))? (n1) : (n2))
#define num_max(n1, n2) (((n1) > (n2))? (n1) : (n2))
/**
 * List
 */

#define define_list_core(list_type, e_type, size_type, inc_size)	\
	\
typedef struct { e_type* buffer; size_type size; size_type cap; } list_type;	\
	\
static inline list_type* init_##list_type(size_type init_size){	\
	if(init_size == 0) init_size = 2;	\
	list_type *list = (list_type*)malloc(sizeof(list_type));	\
	list->size = 0;	\
	list->cap  = init_size;	\
	list->buffer = (e_type*)malloc(sizeof(e_type) * list->cap);	\
	return list;	\
}	\
	\
static inline void list_type##_init(list_type *list, size_type init_size){	\
	if(init_size == 0) init_size = 2;	\
	list->size = 0;	\
	list->cap  = init_size;	\
	list->buffer = (e_type*)malloc(sizeof(e_type) * list->cap);	\
}	\
	\
static inline size_type count_##list_type(list_type *list){ return list->size; }	\
	\
static inline void clear_##list_type(list_type *list){ list->size = 0; }	\
	\
static inline void encap_##list_type(list_type *list, size_type n){	\
	size_type next_size;	\
	if(list->size + n <= list->cap) return;	\
	if(list->size + n < list->size){	\
		fprintf(stderr, " -- elements size exceed %s's data type %s in %s -- %s:%d --\n", #list_type, #size_type, __FUNCTION__, __FILE__, __LINE__);	\
		fflush(stderr);	\
		abort();	\
	}	\
	next_size = list->cap;	\
	while(list->size + n > list->cap){	\
		if(list->cap < inc_size){	\
			list->cap <<= 1;	\
		} else {	\
			list->cap += inc_size;	\
		}	\
	}	\
	list->buffer = realloc(list->buffer, list->cap * sizeof(e_type));	\
}	\
	\
static inline void trunc_##list_type(list_type *list, size_type size){	\
	if(size > count_##list_type(list)) size = count_##list_type(list);	\
	list->size -= size;	\
}	\
	\
static inline void set_##list_type##_size(list_type *list, size_type size){ list->size = size; }	\
	\
static inline void incre_##list_type(list_type *list, size_type size){	\
	if(size + list->size > list->cap) list->size = list->cap;	\
	else list->size += size;	\
}	\
	\
static inline void push_##list_type(list_type *list, e_type e){	\
	encap_##list_type(list, 1);	\
	list->buffer[list->size++] = e;	\
}	\
	\
static inline int pop_##list_type(list_type *list, e_type*e){	\
	if(count_##list_type(list)){	\
		list->size --;	\
		*e = list->buffer[list->size];	\
		return 1;	\
	} else return 0;	\
}	\
	\
static inline void insert_##list_type(list_type *list, size_type idx, e_type e){	\
	if(idx > list->size) return;	\
	encap_##list_type(list, 1);	\
	if(idx == list->size){	\
		list->buffer[list->size] = e;	\
	} else {	\
		memmove(list->buffer + idx + 1, list->buffer + idx, (list->size - idx) * sizeof(e_type));	\
		list->buffer[idx] = e;	\
	}	\
	list->size ++;	\
}	\
	\
static inline void remove_##list_type(list_type *list, size_type idx){	\
	if(idx >= list->size) return;	\
	if(idx + 1 < list->size){	\
		memmove(list->buffer + idx, list->buffer + idx + 1, (list->size - idx - 1) * sizeof(e_type));	\
	}	\
	list->size --;	\
}	\
	\
static inline void set_##list_type(list_type *list, size_type idx, e_type e){ list->buffer[idx] = e; }	\
	\
static inline e_type get_##list_type(list_type *list, size_type idx){ return list->buffer[idx]; }	\
	\
static inline e_type* ref_##list_type(list_type *list, size_type idx){ return list->buffer + idx; }	\
	\
static inline e_type* next_ref_##list_type(list_type *list){ encap_##list_type(list, 1); list->size ++; return list->buffer + list->size - 1; }	\
	\
static inline e_type* ref_next_##list_type(list_type *list){ list->size ++; return list->buffer + list->size - 1; }	\
	\
static inline e_type* as_array_##list_type(list_type *list){ return list->buffer; }	\
	\
static inline void reverse_##list_type(list_type *list){	\
	size_type i, j;	\
	e_type t;	\
	if(count_##list_type(list) == 0) return;	\
	i = 0;	\
	j = count_##list_type(list) - 1;	\
	while(i < j){	\
		t = get_##list_type(list, i);	\
		set_##list_type(list, i, get_##list_type(list, j));	\
		set_##list_type(list, j, t);	\
		i ++;	\
		j --;	\
	}	\
}	\
	\
static inline void append_##list_type(list_type *list1, list_type *list2){	\
	encap_##list_type(list1, count_##list_type(list2));	\
	memcpy(list1->buffer + list1->size, list2->buffer, sizeof(e_type) * list2->size);	\
	list1->size += list2->size;	\
}	\
	\
static inline size_type dump_##list_type(list_type *list, FILE *out){	\
	return fwrite(list->buffer, sizeof(e_type), count_##list_type(list), out);	\
}	\
	\
static inline void free_##list_type(list_type *list){ free(list->buffer); free(list); }	\
	\
static inline void list_type##_free(list_type *list){ free(list->buffer); list->buffer = NULL; }	\

#define define_list_ext(list_type, e_type, size_type, equals)	\
static inline size_type delete_##list_type(list_type *list, e_type e){	\
	size_type i, ret;	\
	ret = 0;	\
	for(i=list->size;i>0;i--){	\
		if(equals(list->buffer[i-1], e)){	\
			if(i < list->size){	\
				memmove(list->buffer + i - 1, list->buffer + i, (list->size - i) * sizeof(e_type));	\
			}	\
			list->size --;	\
			ret ++;	\
		}	\
	}	\
	return ret;	\
}	\
	\
static inline size_type occ_##list_type(list_type *list, e_type e){	\
	size_type i, n;	\
	for(i=0,n=0;i<list->size;i++){	\
		if(equals(list->buffer[i], e)) n++;	\
	}	\
	return n;	\
}	\
	\
static inline size_type replace_##list_type(list_type *list, e_type from, e_type to){	\
	size_type i, ret;	\
	ret = 0;	\
	for(i=0;i<list->size;i++){	\
		if(equals(list->buffer[i], from)){	\
			list->buffer[i] = to;	\
			ret ++;	\
		}	\
	}	\
	return ret;	\
}	\
	\
static inline size_type locate_##list_type(list_type *list, e_type e, size_type start){	\
	size_type i;	\
	for(i=start;i<list->size;i++){	\
		if(equals(list->buffer[i], e)) return i;	\
	}	\
	return i;	\
}

#define define_list(name, e_type) define_list_core(name, e_type, size_t, 0xFFFFFU)

#define native_number_equals(e1, e2) ((e1) == (e2))

#define define_native_list(name, e_type)	\
define_list_core(name, e_type, size_t, 0xFFFFFU);	\
define_list_ext(name, e_type, size_t, native_number_equals);

define_native_list(u8list,  uint8_t);
define_native_list(u16list, uint16_t);
define_native_list(u32list, uint32_t);
define_native_list(u64list, uint64_t);

define_native_list(b8list,  int8_t);
define_native_list(b16list, int16_t);
define_native_list(b32list, int32_t);
define_native_list(b64list, int64_t);

define_list(vplist, void*);

#endif

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
 
#ifndef __THEAD_RJ_H
#define __THEAD_RJ_H

#include <pthread.h>
#include <assert.h>
#include <stdlib.h>
#include "timer.h"

#define THREAD_LOCATION_FUNC	0
#define THREAD_LOCATION_INIT	1
#define THREAD_LOCATION_MAIN	2
#define THREAD_LOCATION_CLOSE	3

#define thread_begin_def(tname)	\
struct tname##_struct {		\
	struct tname##_struct *tname##_params;	\
	int n_cpu;	\
	int t_idx;	\
	int running;	\
	int state;	\
	pthread_mutex_t *mutex_lock;	\
	int sig_pipe
#define thread_end_def(tname) }

#define thread_beg_def(tname) thread_begin_def(tname)

#define thread_begin_func(tname) void* thread_##tname##_func(void *obj){\
	struct tname##_struct * tname = (struct tname##_struct *)obj;\
	int tname##_where = THREAD_LOCATION_FUNC;\
	int tname##_var_i;\
	struct tname##_struct * tname##_params
#define thread_beg_func(tname) thread_begin_func(tname)
#define thread_begin_loop(tname) tname##_params = tname->tname##_params;\
	while(tname->running){\
		if(tname->state == 0){ microsleep(1); continue; }\
		for(tname##_var_i=0;tname##_var_i<1;tname##_var_i++){
#define thread_beg_loop(tname) thread_begin_loop(tname)
#define thread_begin_syn(tname) pthread_mutex_lock(tname->mutex_lock)
#define thread_beg_syn(tname) thread_begin_syn(tname)
#define thread_end_syn(tname) pthread_mutex_unlock(tname->mutex_lock)
#define thread_end_loop(tname) } tname->state = 0; }
#define thread_end_func(tname) return NULL; tname##_where = tname##_where; }

#define thread_preprocess(tname) struct tname##_struct * tname##_params = NULL;\
	struct tname##_struct * tname = NULL;\
	pthread_mutex_t tname##_mlock = PTHREAD_MUTEX_INITIALIZER;\
	pthread_t * tname##_pids = NULL;\
	int tname##_i, tname##_j;\
	int tname##_where = 0;	\
	Timer * tname##_timer = NULL
#define thread_begin_init(tname, n_thread) tname##_where = THREAD_LOCATION_INIT;\
	assert(n_thread > 0);\
	tname##_params = (struct tname##_struct *)malloc(sizeof(struct tname##_struct) * n_thread);\
	tname##_pids = (pthread_t *)malloc(sizeof(pthread_t) * n_thread); \
	tname##_timer = init_timer();\
	for(tname##_i=0,tname##_j=0;tname##_i<n_thread;tname##_i++){ \
		tname = tname##_params + tname##_i;\
		tname->mutex_lock = &tname##_mlock;\
		tname->n_cpu      = n_thread;\
		tname->t_idx      = tname##_i;\
		tname->running    = 1;\
		tname->state      = 0;\
		tname->sig_pipe   = tname##_timer->pipes[1];\
		tname->tname##_params = tname##_params
#define thread_beg_init(tname, n_thread) thread_begin_init(tname, n_thread)

#define thread_end_init(tname) if(pthread_create(tname##_pids + tname##_i, NULL, thread_##tname##_func, tname) != 0){\
			fprintf(stderr, " -- Failed to create thread [%s, %04d] in %s -- %s:%d --\n", #tname, tname##_i, __FUNCTION__, __FILE__, __LINE__);\
			abort();\
			}\
		}\
		tname##_where = THREAD_LOCATION_MAIN

#define thread_begin_operate(tname, idx) tname = tname##_params + idx
#define thread_beg_operate(tname, idx) thread_begin_operate(tname, idx)
#define thread_sleep(tname) tname->state = 0
#define thread_wake(tname)  tname->state = 1
#define thread_waitfor_idle(tname) while(tname->state > 0){ microsleep(1); }
#define thread_end_operate(tname, idx)   tname = NULL
#define thread_begin_iter(tname) for(tname##_i=0;tname##_i<tname##_params[0].n_cpu;tname##_i++){ tname = tname##_params + tname##_i
#define thread_beg_iter(tname) thread_begin_iter(tname)
#define thread_is_idle(tname) (tname->state == 0)
#define thread_n_cpus(tname) (tname->n_cpu)
#define thread_index(tname) (tname->t_idx)
#define thread_end_iter(tname) tname = NULL; }
#define thread_access(tname, idx) (tname##_params + idx)
#define thread_waitfor_one_idle(tname)	\
while(1){	\
	for(tname##_j=0;tname##_j<tname##_params[0].n_cpu;tname##_j++){	\
		if(tname##_params[tname##_j].state == 0){	\
			tname = tname##_params + tname##_j;	\
			break;	\
		}	\
	}	\
	if(tname##_j == tname##_params[0].n_cpu){	\
		microsleep(1);	\
	} else break;	\
}

#define thread_waitfor_all_idle(tname) thread_begin_iter(tname); thread_waitfor_idle(tname); thread_end_iter(tname)

#define thread_begin_close(tname) for(tname##_i=0;tname##_i<tname##_params[0].n_cpu;tname##_i++){ \
		tname= tname##_params + tname##_i;\
		tname->running = 0; \
		pthread_join(tname##_pids[tname##_i], NULL);\
		tname##_where = THREAD_LOCATION_CLOSE
#define thread_beg_close(tname) thread_begin_close(tname)
#define thread_end_close(tname) } free(tname##_params); free(tname##_pids); free_timer(tname##_timer)


#endif

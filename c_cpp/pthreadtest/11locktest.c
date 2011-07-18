/* mutexex.c   */
/* Simple pthread example using pthread_mutex to ensure mutual exclusion */
/* This corrects the bug from raceexample.c                              */
/* To compile me for Linux, type:  gcc -o filename filename.c -lpthread  */
/* To execute, type:  filename                                           */
#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>   //getrusage, gettimeofday
#include <sys/resource.h>

void * simplemux1(void *);
void * simplemux2(void *);
void * simplemux3(void *);
void * simplemux4(void *);
void * simplemux5(void *);
void * simplespin(void *);
void * simplespins(void *);
void * simplecas1(void *);
void * simplecas2(void *);
void * simplecas3(void *);
void * simplecas4(void *);
void * simplecas5(void *);
void * simplena(void *);

int bignum,NumThreads,SumCount;
pthread_mutex_t mut1=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mut2=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
pthread_mutex_t mut3=PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP;
pthread_mutex_t mut4=PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP;
pthread_mutex_t mut5={ { 0, 0, 0, 0, PTHREAD_MUTEX_FAST_NP, 0, { 0, 0 } } };
pthread_spinlock_t spin,spins;

struct timeval __g_timeofday_start, __g_timeofday_end, __g_timeofday_diff,\
       __g_resource_usage_tmp,__g_resource_usage_other_time;
struct rusage __g_resource_usage;

#define G_TIMER_START \
    gettimeofday(&__g_timeofday_start,NULL)

#define G_TIMER_END \
    do {\
        getrusage(RUSAGE_SELF, &__g_resource_usage);\
        gettimeofday(&__g_timeofday_end,NULL);\
        timersub(&__g_timeofday_end,&__g_timeofday_start,&__g_timeofday_diff);\
    } while (0)

#define G_TIMER_PRINT \
    fprintf(stderr,\
        "U: %ld.%06ld, S: %ld.%06ld, R: %ld.%06ld (s)\n",\
        __g_resource_usage.ru_utime.tv_sec, __g_resource_usage.ru_utime.tv_usec,\
        __g_resource_usage.ru_stime.tv_sec, __g_resource_usage.ru_stime.tv_usec,\
        __g_timeofday_diff.tv_sec, __g_timeofday_diff.tv_usec);

int main( int argc, char *argv[] ) {
  int i;

  if (argc != 3) {
    printf("Usage: %s threads to_count\n",argv[0]);
    exit(1);
  } else {
    SumCount=atoi(argv[2]);
    NumThreads=atoi(argv[1]);
	fprintf(stderr, "Threads:%d, To_Count:%d\n\n",NumThreads,SumCount);
  }
  pthread_t *tid=malloc(NumThreads*sizeof(pthread_t));      /* array of thread IDs */
  //pthread_mutex_init(&mut1, NULL);
  pthread_spin_init(&spin, PTHREAD_PROCESS_PRIVATE);
  pthread_spin_init(&spins, PTHREAD_PROCESS_SHARED);

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplena, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "Non-lock:\n[%d]%d   ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
fputs("\n", stderr);

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplespin, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "Spin-lock:\n[%d]%d PROCESS_PRIVATE ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplespins, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d PROCESS_SHARED  ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
fputs("\n", stderr);
G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux1, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "Mux-lock:\n[%d]%d NULL       ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux2, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d RECURSIVE  ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux3, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d ERRORCHECK ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux4, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d ADAPTIVE   ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux5, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d FAST       ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
fputs("\n", stderr);
G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas1, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "Atomic:\n[%d]%d fetch_and_add ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas2, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d add_and_fetch ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas3, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d val cmpxchg   ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas4, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  fprintf(stderr, "[%d]%d bool cmpxchg  ", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
exit(0);
}  /* main */

  
void * simplemux1(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut1);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut1);
  }
  return NULL;
}
  
void * simplemux2(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut2);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut2);
  }
  return NULL;
}
  
void * simplemux3(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut3);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut3);
  }
  return NULL;
}
  
void * simplemux4(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut4);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut4);
  }
  return NULL;
}
  
void * simplemux5(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut5);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut5);
  }
  return NULL;
}

void * simplespin(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_spin_lock(&spin);   
       bignum++;  /* critical section */
     pthread_spin_unlock(&spin);
  }
  return NULL;
}

void * simplespins(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_spin_lock(&spins);   
       bignum++;  /* critical section */
     pthread_spin_unlock(&spins);
  }
  return NULL;
}

void * simplecas1(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
       __sync_fetch_and_add(&bignum,1);  /* critical section */
  }
  return NULL;
}
void * simplecas2(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
       __sync_add_and_fetch(&bignum,1);  /* critical section */
  }
  return NULL;
}
void * simplecas3(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
	  int old,ret;
	  do {
		  old=bignum;
		  ret=__sync_val_compare_and_swap(&bignum,old,old+1);
		  //if (ret!=old) printf("O:%d R:%d\n",old,ret);
	  } while (ret!=old);
  }
  return NULL;
}
void * simplecas4(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
      int old;
      do {
        old=bignum;
      } while (!__sync_bool_compare_and_swap(&bignum,old,old+1));
  }
  return NULL;
}

#pragma GCC optimize 0
void * simplena(void * parm)
{ 
  register int i;
  for(i=0;i<SumCount;++i) {
       ++bignum;  /* critical section */
  }
  return NULL;
}
#pragma GCC optimize 3


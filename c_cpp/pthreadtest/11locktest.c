/* mutexex.c   */
/* Simple pthread example using pthread_mutex to ensure mutual exclusion */
/* This corrects the bug from raceexample.c                              */
/* To compile me for Linux, type:  gcc -o filename filename.c -lpthread  */
/* To execute, type:  filename                                           */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>   //getrusage, gettimeofday
#include <sys/resource.h>

void * simplemux(void *);
void * simplespin(void *);
void * simplespins(void *);
void * simplecas1(void *);
void * simplecas2(void *);
void * simplecas3(void *);
void * simplena(void *);

int bignum,NumThreads,SumCount;
pthread_mutex_t mut;
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
        "User: %ld.%06ld(s), System: %ld.%06ld(s). Real: %ld.%06ld(s).\n\n",\
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
	printf("Threads:%d, To_Count:%d\n",NumThreads,SumCount);
  }
  pthread_t *tid=malloc(NumThreads*sizeof(pthread_t));      /* array of thread IDs */
  pthread_mutex_init(&mut, NULL);
  pthread_spin_init(&spin, PTHREAD_PROCESS_PRIVATE);
  pthread_spin_init(&spins, PTHREAD_PROCESS_SHARED);

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplena, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]Non-lock: bignum=%d\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplespin, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]Spin-lock: bignum=%d   PROCESS_PRIVATE\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplespins, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]Spin-lock: bignum=%d   PROCESS_SHARED\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplemux, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]Mux-lock: bignum=%d\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas1, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]CAS1: bignum=%d\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas2, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]CAS2: bignum=%d\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;

G_TIMER_START;
  bignum = 0;
  for (i=0; i<NumThreads; i++) {
    pthread_create(&tid[i], NULL, simplecas3, NULL);
  }
  for ( i = 0; i < NumThreads; i++)
    pthread_join(tid[i], NULL);
  printf("[%d]CAS3: bignum=%d\n", i, bignum);
G_TIMER_END;
G_TIMER_PRINT;
exit(0);
}  /* main */

  

void * simplemux(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
     pthread_mutex_lock(&mut);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut);
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
	  int value,old,ret;
	  do {
		  old=bignum;
		  //value=old+1;
		  ret=__sync_val_compare_and_swap(&bignum,old,old+1);
		  if (ret!=old) printf("O:%d N:%d\n");
	  } while (ret!=old);
  }
  return NULL;
}
void * simplena(void * parm)
{ 
  int i;
  for(i=0;i<SumCount;i++) {
       bignum++;  /* critical section */
  }
  return NULL;
}

/* mutexex.c   */
/* Simple pthread example using pthread_mutex to ensure mutual exclusion */
/* This corrects the bug from raceexample.c                              */
/* To compile me for Linux, type:  gcc -o filename filename.c -lpthread  */
/* To execute, type:  filename                                           */

#include <pthread.h>
#include <stdio.h>

void * simple(void *);

#define NUM_THREADS 2
pthread_t tid[NUM_THREADS];      /* array of thread IDs */

int bignum = 0;
pthread_mutex_t mut;

main( int argc, char *argv[] ) 
{
  int i, ret;

  pthread_mutex_init(&mut, NULL);

  for (i=0; i<NUM_THREADS; i++) {
    pthread_create(&tid[i], NULL, simple, NULL);
  }
  for ( i = 0; i < NUM_THREADS; i++)
    pthread_join(tid[i], NULL);

  printf("main() reporting that all %d threads have terminated\n", i);
  printf("I am main! bignum=%d\n", bignum);

}  /* main */

  

void * simple(void * parm)
{ 
  int i;
  for(i=0;i<10000;i++) {
     pthread_mutex_lock(&mut);   
       bignum++;  /* critical section */
     pthread_mutex_unlock(&mut);
  }
}    


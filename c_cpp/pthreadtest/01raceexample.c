/* raceexample.c               */
/* Simple pthread race example */
/* THIS CODE HAS A BUG IN IT!  */
/* To compile me for Linux, type:  gcc -o filename filename.c -lpthread */
/* To execute, type:  filename */

#include <pthread.h>
#include <stdio.h>

void * adder(void *);

#define NUM_THREADS 2
pthread_t tid[NUM_THREADS];      /* array of thread IDs */

int bignum = 0;

main( int argc, char *argv[] ) 
{
  int i, ret;

  //thr_setconcurrency(3);
  for (i=0; i<NUM_THREADS; i++) {
    pthread_create(&tid[i], NULL, adder, NULL);
  }
  for ( i = 0; i < NUM_THREADS; i++)
    pthread_join(tid[i], NULL);

  printf("main() reporting that all %d threads have terminated\n", i);
  printf("I am main! bignum=%d\n", bignum);

}  /* main */

  

void * adder(void * parm)
{
   int i;
   printf("I am a new thread!\n");
   for(i=0;i<100000;i++) {
       bignum++;   /* BUG HERE: THIS IS NOT IN A MUTEX AND IS INCORRECT!! */
   }
   pthread_exit(0);
}    


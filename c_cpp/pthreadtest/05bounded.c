/* bounded.c                                                              */
/* Code for Producer/Consumer problem using mutex and condition variables */
/* To compile me for Unix, type:  gcc -o filename filename.c -lpthread */

#include <pthread.h>
#include <stdio.h>

#define BSIZE 4
#define NUMITEMS 30

typedef struct {
  char buf[BSIZE];
  int occupied;
  int nextin, nextout;
  pthread_mutex_t mutex;
  pthread_cond_t more;
  pthread_cond_t less;
} buffer_t;

buffer_t buffer;

void * producer(void *);
void * consumer(void *);


#define NUM_THREADS 2
pthread_t tid[NUM_THREADS];      /* array of thread IDs */

main( int argc, char *argv[] ) 
{
  int i;

  pthread_cond_init(&(buffer.more), NULL);
  pthread_cond_init(&(buffer.less), NULL);

  pthread_create(&tid[1], NULL, consumer, NULL);
  pthread_create(&tid[0], NULL, producer, NULL);
  for ( i = 0; i < NUM_THREADS; i++)
    pthread_join(tid[i], NULL);

  printf("\nmain() reporting that all %d threads have terminated\n", i);

}  /* main */

  

void * producer(void * parm)
{
  char item[NUMITEMS]="IT'S A SMALL WORLD, AFTER ALL.";
  int i;

  printf("producer started.\n");

  for(i=0;i<NUMITEMS;i++)
    { /* produce an item, one character from item[] */

      if (item[i] == '\0') break;  /* Quit if at end of string. */

      pthread_mutex_lock(&(buffer.mutex));

      if (buffer.occupied >= BSIZE) printf("producer waiting.\n");
      while (buffer.occupied >= BSIZE)
	pthread_cond_wait(&(buffer.less), &(buffer.mutex) );
      printf("producer executing.\n");

      buffer.buf[buffer.nextin++] = item[i];
      buffer.nextin %= BSIZE;
      buffer.occupied++;

      /* now: either buffer.occupied < BSIZE and buffer.nextin is the index
	 of the next empty slot in the buffer, or
	 buffer.occupied == BSIZE and buffer.nextin is the index of the
	 next (occupied) slot that will be emptied by a consumer
	 (such as buffer.nextin == buffer.nextout) */

      pthread_cond_signal(&(buffer.more));
      pthread_mutex_unlock(&(buffer.mutex));
    }
  printf("producer exiting.\n");
  pthread_exit(0);
}    

void * consumer(void * parm)
{
  char item;
  int i;

  printf("consumer started.\n");

  for(i=0;i<NUMITEMS;i++){

  pthread_mutex_lock(&(buffer.mutex) );

  if (buffer.occupied <= 0) printf("consumer waiting.\n");
  while(buffer.occupied <= 0)
    pthread_cond_wait(&(buffer.more), &(buffer.mutex) );
  printf("consumer executing.\n");

  item = buffer.buf[buffer.nextout++];
  printf("%c\n",item);
  buffer.nextout %= BSIZE;
  buffer.occupied--;

  /* now: either buffer.occupied > 0 and buffer.nextout is the index
     of the next occupied slot in the buffer, or
     buffer.occupied == 0 and buffer.nextout is the index of the next
     (empty) slot that will be filled by a producer (such as
     buffer.nextout == buffer.nextin) */

  pthread_cond_signal(&(buffer.less));
  pthread_mutex_unlock(&(buffer.mutex));
  }
  printf("consumer exiting.\n");
  pthread_exit(0);
}



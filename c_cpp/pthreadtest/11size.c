/* mutexex.c   */
/* Simple pthread example using pthread_mutex to ensure mutual exclusion */
/* This corrects the bug from raceexample.c                              */
/* To compile me for Linux, type:  gcc -o filename filename.c -lpthread  */
/* To execute, type:  filename                                           */
#define _GNU_SOURCE
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>


pthread_mutex_t mut1=PTHREAD_MUTEX_INITIALIZER;
pthread_spinlock_t spin1;
pthread_cond_t cond1;

int main(void) {
    puts("Sizes:");
    printf("pthread_cond_t\t%zd\n",sizeof(cond1));
    printf("pthread_spinlock_t\t%zd\n",sizeof(spin1));
    printf("pthread_mutex_t\t%zd\n",sizeof(mut1));
    return(0);
}

/*
 * cond.c
 *
 * Demonstrate a simple condition variable wait.
 * 
 * to compile me for Linux, type:  gcc -o conditionexample conditionexample.c -lpthread
 */
#include <pthread.h>
#include <stdio.h>

typedef struct {
    pthread_mutex_t     mutex;  /* Protects access to value */
    pthread_cond_t      cond;   /* Signals change to value */
    int                 value;  /* Access protected by mutex */
} cond_struct_t;

cond_struct_t data = {
    PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, 0};


/*
 * Thread start routine. It will set the main thread's predicate
 * and signal the condition variable.
 */
void * wait_thread (int *hibernation)
{
    int status;

    sleep (*hibernation);

    /****** ENTER THE MONITOR HERE *******/
    status = pthread_mutex_lock (&data.mutex);
    if (status != 0) {fprintf(stderr,"Lock error!\n"); exit(-1);}

    data.value = 1;             /* Set predicate */

    status = pthread_cond_signal (&data.cond);
    if (status != 0) {fprintf(stderr,"Signal error!\n"); exit(-1);}

    status = pthread_mutex_unlock (&data.mutex);
    if (status != 0) {fprintf(stderr,"Unlock error!\n"); exit(-1);}
    /******* EXIT the MONITOR ************/

    return NULL;
}

int main (int argc, char *argv[])
{
    int status;
    pthread_t wait_thread_id;
    struct timespec timeout;
    int hibernation = 1;

    /*
     * If an argument is specified, interpret it as the number
     * of seconds for wait_thread to sleep before signaling the
     * condition variable.  You can play with this to see the
     * condition wait below time out or wake normally.
     */
    //thr_setconcurrency(2);
    if (argc > 1)
        hibernation = atoi (argv[1]);

    /*
     * Create wait_thread.
     */
    status = 
       pthread_create (&wait_thread_id, NULL, (void *) wait_thread, &hibernation);
    if (status != 0) {fprintf(stderr,"Create error!\n"); exit(-1);}

    /*
     * Wait on the condition variable until signaled by 
     * the wait_thread.
     */

    /***** ENTER the MONITOR HERE **********/
    status = pthread_mutex_lock (&data.mutex);
    if (status != 0) {fprintf(stderr,"Lock error!\n"); exit(-1);}

    while (data.value == 0) {
        status = pthread_cond_wait(&data.cond, &data.mutex);
        }
    if (status != 0) {fprintf(stderr,"Wait error!\n"); exit(-1);}

    status = pthread_mutex_unlock (&data.mutex);
    if (status != 0) {fprintf(stderr,"Unlock error!\n"); exit(-1);}
    /******* EXIT the MONITOR **************/

    printf("Condition was signalled, data.value=%d\n",data.value);
    return 0;
}

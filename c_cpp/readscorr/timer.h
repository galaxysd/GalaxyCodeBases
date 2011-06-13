// by Hu Xuesong
#ifndef _G_TIMER_H
#define _G_TIMER_H

#include <sys/time.h>   //getrusage, gettimeofday
#include <sys/resource.h>

struct timeval __g_timeofday_start, __g_timeofday_end, __g_timeofday_diff;
struct rusage __g_resource_usage;

#define G_TIMER_START \
do { gettimeofday(&__g_timeofday_start,NULL); } while (0)

#define G_TIMER_END \
do {\
    getrusage(RUSAGE_SELF, &__g_resource_usage);\
    gettimeofday(&__g_timeofday_end,NULL);\
    timersub(&__g_timeofday_end,&__g_timeofday_start,&__g_timeofday_diff);\
    fprintf(stderr,"\n--------------------------------------------------------------------------------\n"\
        "Resource Usage Measures:\n"\
        "   User:%ld.%06ld s, System:%ld.%06ld s. Real:%ld.%06ld s\n"\
        "   Block I/O times: %ld/%ld. MaxRSS: %ld KiB\n"\
        "   Wait times: %ld(nvcsw) + %ld(nivcsw). "\
        "Page Faults times: %ld(minflt) + %ld(majflt)\n",\
        __g_resource_usage.ru_utime.tv_sec, __g_resource_usage.ru_utime.tv_usec,\
        __g_resource_usage.ru_stime.tv_sec, __g_resource_usage.ru_stime.tv_usec,\
        __g_timeofday_diff.tv_sec, __g_timeofday_diff.tv_usec,\
        __g_resource_usage.ru_inblock, __g_resource_usage.ru_oublock,\
        __g_resource_usage.ru_maxrss,\
        __g_resource_usage.ru_nvcsw, __g_resource_usage.ru_nivcsw,\
        __g_resource_usage.ru_minflt, __g_resource_usage.ru_majflt);\
} while (0)

#endif /* timer.h */

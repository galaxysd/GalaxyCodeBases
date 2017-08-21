// by Hu Xuesong
#ifndef _G_TIMER_H
#define _G_TIMER_H

#include <sys/time.h>   //getrusage, gettimeofday
#include <sys/resource.h>

#ifdef __cplusplus
extern "C" {
#endif

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
        timersub(&__g_timeofday_diff,&__g_resource_usage.ru_utime,&__g_resource_usage_tmp);\
        timersub(&__g_resource_usage_tmp,&__g_resource_usage.ru_stime,&__g_resource_usage_other_time);\
    } while (0)

#define __G_TIMER_PRINT(__sustype__) \
    fprintf(stderr,"\n--------------------------------------------------------------------------------\n"\
        "Resource Usage Measures:\n"\
        "   User: %ld."__sustype__"(s), System: %ld."__sustype__"(s). Real: %ld."__sustype__"(s).\n"\
        "   Sleep: %ld."__sustype__"(s). Block I/O times: %ld/%ld. MaxRSS: %ld kiB.\n"\
        "   Wait(s): %ld(nvcsw) + %ld(nivcsw). "\
        "Page Fault(s): %ld(minflt) + %ld(majflt).\n",\
        __g_resource_usage.ru_utime.tv_sec, __g_resource_usage.ru_utime.tv_usec,\
        __g_resource_usage.ru_stime.tv_sec, __g_resource_usage.ru_stime.tv_usec,\
        __g_timeofday_diff.tv_sec, __g_timeofday_diff.tv_usec,\
        __g_resource_usage_other_time.tv_sec, __g_resource_usage_other_time.tv_usec,\
        __g_resource_usage.ru_inblock, __g_resource_usage.ru_oublock,\
        __g_resource_usage.ru_maxrss,\
        __g_resource_usage.ru_nvcsw, __g_resource_usage.ru_nivcsw,\
        __g_resource_usage.ru_minflt, __g_resource_usage.ru_majflt)

#ifdef __linux__
	#define G_TIMER_PRINT __G_TIMER_PRINT("%06ld")
#elif __APPLE__
	#define G_TIMER_PRINT __G_TIMER_PRINT("%06d")
#else
#   error "Unknown compiler !"
#endif

#ifdef __cplusplus
}
#endif

#endif /* timer.h */

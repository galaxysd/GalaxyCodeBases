#ifndef _G_SDLEFT_TF_H
#define _G_SDLEFT_TF_H
/*
Well, let's try "Template" in C with #define and pointer to function.
Since the function body is the same, I prefer to call it "template" instead of "overload".

Here, we select typeof(Count) on SUM(Count)==dleftobj->ItemInsideAll
*/
#include "gtypendef.h"
#include <stddef.h> //size_t
#include <stdint.h> //uint64_t
#include <stdio.h>  //FILE
#include "sdleft.h"

typedef struct __SDLeftStat_t {
    uint64_t HistMaxCntVal;
    uint64_t HistMaxHistVal;
    double HistMean;
    double HistSStd;
} SDLeftStat_t;

typedef SDLeftStat_t *(G_SDLeftArray_IN)(SDLeftArray_t * const, FILE *);    //*G_SDLeftArray_IN() is OK,too .

SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream);

#endif  // sdleftTF.h

#ifndef _G_SDLEFT_TF_H
#define _G_SDLEFT_TF_H
/*
Well, let's try "Template" in C with #define and pointer to function.
Since the function body is the same, I prefer to call it "template" instead of "overload".
*/

#ifndef _G_UINT128_T
#define _G_UINT128_T
typedef unsigned int uint128_t __attribute__((mode(TI)));
#endif

typedef struct __SDLeftStat_t {
    size_t Count;
    uint64_t Sum;
} SDLeftStat_t;

typedef SDLeftStat_t * (*G_SDLeftArray_IN)(SDLeftArray_t * const);

#endif  // sdleftTF.h

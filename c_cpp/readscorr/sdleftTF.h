#ifndef _G_SDLEFT_TF_H
#define _G_SDLEFT_TF_H
/*
Well, let's try "Template" in C with #define and pointer to function.
Since the function body is the same, I prefer to call it "template" instead of "overload".
*/

typedef void (*G_SDLeftArray_IN)(const SDLeftArray_t * const);
uint16_t *dleft_statu16(SDLeftArray_t *dleftobj);
uint32_t *dleft_statu32(SDLeftArray_t *dleftobj);
uint64_t *dleft_statu64(SDLeftArray_t *dleftobj);

#endif  // sdleftTF.h

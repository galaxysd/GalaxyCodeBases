#ifndef _G_TYPEnDEF_H
#define _G_TYPEnDEF_H

#ifndef _G_UINT128_T
#define _G_UINT128_T
//typedef unsigned __int128 uint128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));
#endif

#ifndef _G_FLOAT128_T
#define _G_FLOAT128_T
typedef __float128 float128 __attribute__((mode(TF)));
#endif

#if defined(_MSC_VER)
#define FORCE_INLINE	__forceinline
#else	// defined(_MSC_VER)
#define	FORCE_INLINE static inline __attribute__((always_inline))
#endif // !defined(_MSC_VER)




/*
http://stackoverflow.com/questions/1188939/representing-128-bit-numbers-in-c
typedef unsigned int uint128_t __attribute__((mode(TI)));

http://gcc.gnu.org/gcc-4.4/changes.html
IA-32/IA64:
Support for __float128 (TFmode) IEEE quad type and corresponding TCmode IEEE complex quad type is available via the soft-fp library on IA-32/IA64 targets. This includes basic arithmetic operations (addition, subtraction, negation, multiplication and division) on __float128 real and TCmode complex values, the full set of IEEE comparisons between __float128 values, conversions to and from float, double and long double floating point types, as well as conversions to and from signed or unsigned integer, signed or unsigned long integer and signed or unsigned quad (TImode, IA64 only) integer types. Additionally, all operations generate the full set of IEEE exceptions and support the full set of IEEE rounding modes.
MIPS:
    asm ("dmultu\t%1,%2" : "=h" (result) : "r" (x), "r" (y));
You can now achieve the same effect using 128-bit types:
    typedef unsigned int uint128_t __attribute__((mode(TI)));
    result = ((uint128_t) x * y) >> 64;

http://gcc.gnu.org/gcc-4.3/changes.html
IA-32/x86-64:
Support for __float128 (TFmode) IEEE quad type and corresponding TCmode IEEE complex quad type is available via the soft-fp library on x86_64 targets. This includes basic arithmetic operations (addition, subtraction, negation, multiplication and division) on __float128 real and TCmode complex values, the full set of IEEE comparisons between __float128 values, conversions to and from float, double and long double floating point types, as well as conversions to and from signed or unsigned integer, signed or unsigned long integer and signed or unsigned quad (TImode) integer types. Additionally, all operations generate the full set of IEEE exceptions and support the full set of IEEE rounding modes.

http://gcc.gnu.org/ml/gcc-patches/2010-04/msg00354.html
[patch]: Add support of new __int128 type for targets having 128-bit integer scalar support
http://gcc.gnu.org/ml/gcc-patches/2010-04/msg00354/int128doc.diff
+  static const char *const suffixes[] = { "", "U", "L", "UL", "LL", "ULL", "I128", "UI128" };

http://en.wikipedia.org/wiki/Floating_point
Type    Sign 	Exponent 	Significand 	Total bits 		Exponent bias 	Bits precision
Single 	1   	8 	        23          	32 	        	127 	        24
Double 	1   	11      	52          	64 	        	1023        	53
Quad 	1   	15      	112          	128     		16383       	113
*/
#endif  // gtypendef.h


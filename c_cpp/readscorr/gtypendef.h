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












#endif  // gtypendef.h


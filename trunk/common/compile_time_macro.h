#ifndef __compile_time_macro_h__
#define __compile_time_macro_h__

#include <climits> // needed for CHAR_BIT definition

#define UNIQ_NAME_HLP2(x,y) x##y
#define UNIQ_NAME_HLP1(x,y) UNIQ_NAME_HLP2(x,y)
#define UNIQ_NAME(p) UNIQ_NAME_HLP1(p,__LINE__)

#define COMPILE_TIME_ASSERT(pred)             switch(0){case 0:case pred:;}

#define ASSERT_MIN_BITSIZE(type, size)        COMPILE_TIME_ASSERT(sizeof(type) * CHAR_BIT >= (size))

#define ASSERT_EXACT_BITSIZE(type, size)      COMPILE_TIME_ASSERT(sizeof(type) * CHAR_BIT == (size))

#define COMPILE_TIME_ASSERT_GLOB(pred, l)     namespace UNIQ_NAME (ns_) { static void UNIQ_NAME (ff_) () { COMPILE_TIME_ASSERT(pred) }  }

#define ASSERT_MIN_BITSIZE_GLOB(type, size)   COMPILE_TIME_ASSERT_GLOB(sizeof(type) * CHAR_BIT >= (size), __LINE__)

#define ASSERT_EXACT_BITSIZE_GLOB(type, size) COMPILE_TIME_ASSERT_GLOB(sizeof(type) * CHAR_BIT == (size), __LINE__)

#endif

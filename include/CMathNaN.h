#ifndef CMathNaN_H
#define CMathNaN_H

#ifndef __linux__

/*****************************************/
/* double Not-a-Number macros            */
/*****************************************/

typedef union {
  struct {
#ifdef _MIPSEL
    unsigned fraction_low:32;
    unsigned bits        :20;
    unsigned exponent    :11;
    unsigned sign        : 1;
#else
    unsigned sign        : 1;
    unsigned exponent    :11;
    unsigned bits        :20;
    unsigned fraction_low:32;
#endif
  } inf_parts;
  struct {
#ifdef _MIPSEL
    unsigned fraction_low:32;
    unsigned bits        :19;
    unsigned qnan_bit    : 1;
    unsigned exponent    :11;
    unsigned sign        : 1;
#else
    unsigned sign        : 1;
    unsigned exponent    :11;
    unsigned qnan_bit    : 1;
    unsigned bits        :19;
    unsigned fraction_low:32;
#endif
  } nan_parts;
  double d;
} dnan;

#define IsInf(X)     (IsNaN(X) && \
                      ((dnan *)&(X))->inf_parts.bits == 0 && \
                     ((dnan *)&(X))->inf_parts.fraction_low == 0)
#define IsPosInf(X)  (((dnan *)&(X))->nan_parts.sign == 0 && IsInf(X))
#define IsNegInf(X)  (((dnan *)&(X))->nan_parts.sign == 1 && IsInf(X))
#define IsNaN(X)     (((dnan *)&(X))->nan_parts.exponent == 0x7ff)

/* IsPosNaN and IsNegNaN can be used to check the sign of infinities too */

#define IsPosNaN(X)  (((dnan *)&(X))->nan_parts.sign == 0 && IsNaN(X))
#define IsNegNaN(X)  (((dnan *)&(X))->nan_parts.sign == 1 && IsNaN(X))

#define SetNaN(X)    (((dnan *)&(X))->nan_parts.exponent = 0x7ff,\
                      ((dnan *)&(X))->inf_parts.bits = 1)
#define SetPosInf(X) (((dnan *)&(X))->nan_parts.exponent = 0x7ff,\
                      ((dnan *)&(X))->inf_parts.bits =\
                      ((dnan *)&(X))->inf_parts.fraction_low =\
                      ((dnan *)&(X))->nan_parts.sign = 0)
#define SetNegInf(X) (((dnan *)&(X))->nan_parts.exponent =  0x7ff,\
                      ((dnan *)&(X))->inf_parts.bits =\
                      ((dnan *)&(X))->inf_parts.fraction_low = 0,\
                      ((dnan *)&(X))->nan_parts.sign = 1)

/* GETNaNPC gets the leftmost 32 bits of the fraction part */

#define GETNaNPC(dval) \
  (((dnan *)&(dval))->inf_parts.bits        << 12 | \
   ((dnan *)&(dval))->nan_parts.fraction_low>> 20)

#define KILLFPE()    (void) kill(getpid(), 8)
#define KILLNaN(X)   if (NaN(X)) KILLFPE()
#define NaN(X)       IsNaN(X)

/*****************************************/
/* float Not-a-Number macros             */
/*****************************************/

typedef union {
  struct {
#ifdef _MIPSEL
    unsigned bits        :23;
    unsigned exponent    : 8;
    unsigned sign        : 1;
#else
    unsigned sign        : 1;
    unsigned exponent    : 8;
    unsigned bits        :23;
#endif
  } nan_parts;
  float f;
} fnan;

#define IsInfF(X)     (IsNaNF(X) && \
                       ((fnan *)&(X))->nan_parts.bits == 0)
#define IsPosInfF(X)  (((fnan *)&(X))->nan_parts.sign == 0 && IsInfF(X))
#define IsNegInfF(X)  (((fnan *)&(X))->nan_parts.sign == 1 && IsInfF(X))
#define IsNaNF(X)     (((fnan *)&(X))->nan_parts.exponent == 0xff)

/* IsPosNaNF and IsNegNaNF can be used to check the sign of infinities too */

#define IsPosNaNF(X)  (((fnan *)&(X))->nan_parts.sign == 0 && IsNaNF(X))
#define IsNegNaNF(X)  (((fnan *)&(X))->nan_parts.sign == 1 && IsNaNF(X))

#define SetNaNF(X)    (((fnan *)&(X))->nan_parts.exponent = 0xff,\
                       ((fnan *)&(X))->nan_parts.bits = 1)
#define SetPosInfF(X) (((fnan *)&(X))->nan_parts.exponent = 0xff,\
                       ((fnan *)&(X))->nan_parts.bits =\
                       ((fnan *)&(X))->nan_parts.sign = 0)
#define SetNegInfF(X) (((fnan *)&(X))->nan_parts.exponent = 0xff,\
                       ((fnan *)&(X))->nan_parts.bits =\
                       ((fnan *)&(X))->nan_parts.sign = 1)

#define NaNF(X)       IsNaNF(X)
#define FloatNaN(X)   IsNaNF(X)

#else
#include <bits/nan.h>

#define IsNaN(X) isnan(X)
#define IsNaNF(X) isnan((double) (X))
#define IsInf(X) (isinf(X)!=0)
#define IsPosInf(X) (isinf(X)==1)
#define IsNegInf(X) (isinf(X)==-1)
#define NaN(X) isnan(X)
#define SetNaN(X) (X = NAN)
#define SetPosInf(X) (X = HUGE_VAL)
#define SetNegInf(X) (X = -HUGE_VAL)

#endif

#endif

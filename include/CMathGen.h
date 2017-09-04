#ifndef CMATH_GEN_H
#define CMATH_GEN_H

#include <cmath>
#include <sys/types.h>
#include <NaN.h>
#include <CIPoint2D.h>

// Fast conversion from a IEEE 32-bit floating point number F in [0,1] to a
// a 32-bit integer I in [0,2^L-1].

#define SCALED_FLOAT_TO_INT(F,L,I) \
{ float fFloat = F; \
  int iShift = 150 - L - ((*(int*)(&fFloat) >> 23) & 0xFF); \
  if (iShift < 24) { \
    I = ((*(int*)(&fFloat) & 0x007FFFFF) | 0x00800000) >> iShift; \
    if (I == (1 << L)) { \
      I--; \
    } \
  } \
  else { \
    I = 0; \
  } \
}

namespace CMathGen {
  enum Handedness {
    LEFT_HANDEDNESS,
    RIGHT_HANDEDNESS
  };

  enum AxisType2D {
    X_AXIS_2D = (1<<0),
    Y_AXIS_2D = (1<<1),

    XY_AXIS_2D = (X_AXIS_2D | Y_AXIS_2D)
  };

  enum AxisType3D {
    X_AXIS_3D = (1<<0),
    Y_AXIS_3D = (1<<1),
    Z_AXIS_3D = (1<<2),

    XY_AXIS_3D = (X_AXIS_3D | Y_AXIS_3D),
    XZ_AXIS_3D = (X_AXIS_3D | Z_AXIS_3D),
    YZ_AXIS_3D = (Y_AXIS_3D | Z_AXIS_3D),

    XYZ_AXIS_3D = (X_AXIS_3D | Y_AXIS_3D | Z_AXIS_3D)
  };

  enum AxisOrder3D {
    XYZ_AXIS_3D_ORDER,
    XZY_AXIS_3D_ORDER,
    YXZ_AXIS_3D_ORDER,
    YZX_AXIS_3D_ORDER,
    ZXY_AXIS_3D_ORDER,
    ZYX_AXIS_3D_ORDER
  };

  enum IntersectType {
    INTERSECT_NONE    = 0,
    INTERSECT_ALL     = (1<<0),
    INTERSECT_INSIDE  = (1<<1),
    INTERSECT_OUTSIDE = (1<<2),
    INTERSECT_VALID   = (INTERSECT_INSIDE | INTERSECT_OUTSIDE)
  };

  extern double EPSILON_E4;
  extern double EPSILON_E5;
  extern double EPSILON_E6;

  extern float EPSILON_E4_F;
  extern float EPSILON_E5_F;
  extern float EPSILON_E6_F;

  double DSin(double angle);
  double DCos(double angle);
  double DTan(double angle);

  int    fastDistance(int x, int y);
  double fastDistance(double x, double y, double z);

  double pow(double real1, double real2);
  double pow(long real1, long real2);
  double ipow(double value, uint power);

  double log(double real);
  double log10(double real);
  double logN(int base, double real);

  bool isPowerOf(uint base, uint value);
  bool isPowerOf2(uint value);

  int numSetBits(uint x);

  uint clearLastBit(uint x);

  double sqrt(double real);

  double acos(double x);
  double asin(double x);

  double atan2(double x, double y);
  float  atan2(float x, float y);

  long   modulus(long   x, long   y);
  double modulus(double x, double y);

  int sign(long x);
  int sign(double x);

  uint  GetHCF(uint a, uint b);
  ulong GetHCF(ulong a, ulong b);
  uint  GetLCM(uint a, uint b);
  ulong GetLCM(ulong a, ulong b);

  double DegToRad(double deg);
  double RadToDeg(double rad);

  bool isInteger(double d);

  bool cmp(double a, double b, double prec=1E-5);

  uint   twosCompliment(int value);
  ushort twosCompliment(short value);

  bool range(double *x, double *y, uint num_xy,
             double *xmin, double *ymin, double *xmax, double *ymax);

  double pythag(double a, double b);

  float  FastInvSqrtF(float fValue);
  double FastInvSqrt(double dValue);
  double FastSin0(double fAngle);
  double FastSin1(double fAngle);
  double FastCos0(double fAngle);
  double FastCos1(double fAngle);
  double FastTan0(double fAngle);
  double FastTan1(double fAngle);
  double FastInvSin0(double fValue);
  double FastInvSin1(double fValue);
  double FastInvCos0(double fValue);
  double FastInvCos1(double fValue);
  double FastInvTan0(double fValue);
  double FastInvTan1(double fValue);

  //-----

  template<typename T>
  T avg(T v1, T v2) {
    return (v1 + v2)/2;
  }

  //-----

  bool realEq(double v1, double v2, double tol=1E-6);

  //-----

  template<typename T>
  T clamp(T val, T low, T high) {
    if (val < low ) return low;
    if (val > high) return high;
    return val;
  }

  //-----

  template<typename T>
  bool solveQuadratic(T a, T b, T c, T *r1, T *r2) {
    if (a == 0.0) return false;

    T b2_4ac = b*b - 4.0*a*c;

    if (b2_4ac < 0) return false;

    T sqrt_b2_4ac = ::sqrt(b2_4ac);

    T q;

    if (b < 0.0)
      q = -0.5*(b - sqrt_b2_4ac);
    else
      q = -0.5*(b + sqrt_b2_4ac);

    *r1 = q/a;
    *r2 = c/q;

    if (*r1 > *r2) std::swap(*r1, *r2);

    return true;
  }

  //-----

  double mapToReal(char i);
  double mapToReal(unsigned char i);
  double mapToReal(short i);
  double mapToReal(ushort i);
  double mapToReal(int i);
  double mapToReal(uint i);
  double mapToReal(long i);
  double mapToReal(ulong i);

  //-----

  ulong  factorial(ushort f);
  double factorial(double f);

  ulong binomialCoeff(uint n, uint k);

  //-----

  void getXYVals(const CIPoint2D *points, uint npoints,
                 double **xvals, int *num_xvals,
                 double **yvals, int *num_yvals);

  void getXYVals(double *x, double *y, int num_xy,
                 double **xvals, int *num_xvals,
                 double **yvals, int *num_yvals);

  //-----

  double getNaN();

  double noise(double x, double y, double z);
}

#endif

#ifndef CMathGen_H
#define CMathGen_H

#include <NaN.h>

#include <cerrno>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <sys/types.h>

#define PI2       6.283185407
#define PI_DIV_2  1.570796327
#define PI_DIV_4  0.785398163
#define PI3_DIV_2 4.712388980
#define PI_INV    0.318309886

#define PI_F        3.141592654f
#define PI2_F       6.283185407f
#define PI_DIV_2_F  1.570796327f
#define PI_DIV_4_F  0.785398163f
#define PI3_DIV_2_F 4.712388980f
#define PI_INV_F    0.318309886f

#define DEG_TO_RAD(d) (M_PI*(d)/180.0)
#define RAD_TO_DEG(r) ((180.0*(r))/M_PI)

#define REAL_EQ(r1,r2) (fabs((r1)-(r2))<1E-6)

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

  enum IntersectType {
    INTERSECT_NONE    = 0,
    INTERSECT_ALL     = (1<<0),
    INTERSECT_INSIDE  = (1<<1),
    INTERSECT_OUTSIDE = (1<<2),
    INTERSECT_VALID   = (INTERSECT_INSIDE | INTERSECT_OUTSIDE)
  };

  enum Rounding {
    ROUND_DOWN,
    ROUND_UP,
    ROUND_NEAREST
  };

  inline int RoundNearest(double x) {
    double x1;

    if (x <= 0.0)
      x1 = (x - 0.499999);
    else
      x1 = (x + 0.500001);

    if (x1 < INT_MIN || x1 > INT_MAX)
      errno = ERANGE;

    return int(x1);
  }

  inline int RoundUp(double x) {
    double x1;

    if (x <= 0.0)
      x1 = (x       - 1E-6);
    else
      x1 = (x + 1.0 - 1E-6);

    if (x1 < INT_MIN || x1 > INT_MAX)
      errno = ERANGE;

    return int(x1);
  }

  inline int RoundDown(double x) {
    double x1;

    if (x >= 0.0)
      x1 = (x       + 1E-6);
    else
      x1 = (x - 1.0 + 1E-6);

    if (x1 < INT_MIN || x1 > INT_MAX)
      errno = ERANGE;

    return int(x1);
  }

  inline int Round(double x, Rounding rounding=ROUND_NEAREST) {
    switch (rounding) {
      case ROUND_UP  : return RoundUp(x);
      case ROUND_DOWN: return RoundDown(x);
      default        : return RoundNearest(x);
    }
  }

  inline bool cmp(double a, double b, double prec=1E-5) {
    return (std::fabs(a - b) < prec);
  }

  //-----

  inline double DSin(double angle) {
    return sin(DEG_TO_RAD(angle));
  }

  inline double DCos(double angle) {
    return cos(DEG_TO_RAD(angle));
  }

  inline double DTan(double angle) {
    return tan(DEG_TO_RAD(angle));
  }

  //-----

  inline void math_throw(const char *msg) {
    std::cerr << msg << std::endl;
  }

  //-----

  inline long sign(long x) {
    if      (x > 0) return 1;
    else if (x < 0) return -1;
    else            return 0;
  }

  inline long sign(double x) {
    if      (x > 1E-6) return 1;
    else if (x < 1E-6) return -1;
    else               return 0;
  }

  inline double DegToRad(double d) {
    return DEG_TO_RAD(d);
  }

  inline double RadToDeg(double r) {
    return RAD_TO_DEG(r);
  }

  inline double sqrt(double real) {
     double result;

    if (IsNaN(real)) {
      SetNaN(result);
      return result;
    }

    if (real < 0.0)
      math_throw("Negative Value for Sqrt");

    return ::sqrt(real);
  }

  inline bool is_integer(double r) {
    return (long(r) == r);
  }

  inline double pow(double real1, double real2) {
    double real;

    if (IsNaN(real1) || IsNaN(real2)) {
      SetNaN(real);
      return real;
    }

    if (real1 < 0.0 && ! is_integer(real2))
      math_throw("Non Integer Power of Negative");

    if (real1 == 0.0 && real2 < 0.0)
      math_throw("Zero to Negative Power is Undefined");

    errno = 0;

    if (real2 < 0.0)
      real = 1.0/::pow(real1, -real2);
    else
      real = ::pow(real1, real2);

    if (errno != 0)
      math_throw("Power Failed");

    return real;
  }

  inline double pow(long integer1, long integer2) {
    if (integer1 == 0 && integer2 < 0)
      math_throw("Zero to Negative Power is Undefined");

    errno = 0;

    double real;

    if (integer2 < 0.0)
      real = 1.0/::pow((double) integer1, (double) -integer2);
    else
      real = ::pow((double) integer1, (double) integer2);

    if (errno != 0)
      math_throw("Power Failed");

    return real;
  }

  inline double acos(double x) {
    if      (x >=  1.0) return 0.0;
    else if (x <= -1.0) return M_PI;
    else                return ::acos(x);
  }

  inline double asin(double x) {
    if      (x >=  1.0) return  M_PI/2.0;
    else if (x <= -1.0) return -M_PI/2.0;
    else                return ::asin(x);
  }

  inline double log(double real) {
    double result;

    if (IsNaN(real)) {
      SetNaN(result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log(real);

    return result;
  }

  inline double logN(int n, double real) {
    double result;

    if (IsNaN(real)) {
      SetNaN(result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log(real)/log(n);

    return result;
  }

  inline double log10(double real) {
    double result;

    if (IsNaN(real)) {
      SetNaN(result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log10(real);

    return result;
  }

  inline double modulus(double real1, double real2) {
    double result = 0;

    if (IsNaN(real1) || IsNaN(real2)) {
      SetNaN(result);
      return result;
    }

    if (real2 == 0.0)
      math_throw("Divide by Zero");

    int factor = (int) (real1/real2);

    result = real1 - (real2*factor);

    return result;
  }

  inline bool isPowerOf2(uint x) {
    // Complement and Compare
    return ((x != 0) && ((x & (~x + 1)) == x));
  }

  inline bool isPowerOf(uint base, uint value) {
    assert(base > 1);

    if (base == 2) return isPowerOf2(value);

    if (value == 0) return false;

    double n = logN(base, value);

    double r = n - int(n);

    return (r < 1E-10);
  }

  inline double Hypot(double a, double b) {
    double r;

    if      (abs(a) > abs(b)) {
      r = b/a;
      r = abs(a)*sqrt(1 + r*r);
    }
    else if (b != 0) {
      r = a/b;
      r = abs(b)*sqrt(1 + r*r);
    }
    else {
      r = 0.0;
    }

    return r;
  }

  // note: usual argument order of atan2 flipped - ugh
  inline float atan2(float x, float y) {
    if      (fabsf(x) < 1E-6) {
      if      (fabsf(y) < 1E-6)
        return 0.0f;
      else if (y > 0.0f)
        return PI_DIV_2_F;
      else
        return PI3_DIV_2_F;
    }
    else if (fabsf(y) < 1E-6) {
      if      (fabsf(x) < 1E-6)
        return 0.0f;
      else if (x > 0.0f)
        return 0.0f;
      else
        return PI_F;
    }
    else {
      float angle = ::atan2(y, x);

      if (angle >= 0.0f)
        return angle;
      else
        return (PI2_F + angle);
    }
  }

  // note: usual argument order of atan2 flipped - ugh
  inline double atan2(double x, double y) {
    if      (fabs(x) < 1E-6) {
      if      (fabs(y) < 1E-6)
        return 0.0;
      else if (y > 0.0)
        return PI_DIV_2;
      else
        return PI3_DIV_2;
    }
    else if (fabs(y) < 1E-6) {
      if      (fabs(x) < 1E-6)
        return 0.0;
      else if (x > 0.0)
        return 0.0;
      else
        return M_PI;
    }
    else {
      double angle = ::atan2(y, x);

      if (angle >= 0.0)
        return angle;
      else
        return (PI2 + angle);
    }
  }

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

  inline double getNaN() {
    double r;

    SetNaN(r);

    return r;
  }
}

#endif

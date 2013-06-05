#ifndef CMathGen_H
#define CMathGen_H

#include <COSNaN.h>

#include <cerrno>
#include <cmath>
#include <cassert>
#include <iostream>

#define DEG_TO_RAD(d) (M_PI*(d)/180.0)
#define RAD_TO_DEG(r) ((180.0*(r))/M_PI)

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

  inline int Round(double x) {
    if (x >= 0.0) return int(x + 0.5);
    else          return int(x - 0.5);
  }

  inline int RoundDown(double x) {
    if (x >= 0.0) return int(x);
                  return int(x);
  }

  inline int RoundUp(double x) {
    if (x >= 0.0) return int(x + 0.999999);
                  return int(x - 0.999999);
  }

  //-----

  inline void math_throw(const char *msg) {
    std::cerr << msg << std::endl;
  }

  inline long sign(long x) {
    if      (x > 0) return 1;
    else if (x < 0) return -1;
    else            return 0;
  }

  inline double sqrt(double real) {
     double result;

    if (COSNaN::is_nan(real)) {
      COSNaN::set_nan(&result);
      return real;
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

    if (COSNaN::is_nan(real1) || COSNaN::is_nan(real2)) {
      COSNaN::set_nan(&real);
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

    if (COSNaN::is_nan(real)) {
      COSNaN::set_nan(&result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log(real);

    return result;
  }

  inline double logN(int n, double real) {
    double result;

    if (COSNaN::is_nan(real)) {
      COSNaN::set_nan(&result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log(real)/log(n);

    return result;
  }

  inline double log10(double real) {
    double result;

    if (COSNaN::is_nan(real)) {
      COSNaN::set_nan(&result);
      return real;
    }

    if (real <= 0.0)
      math_throw("Negative or Zero Value for Log");

    result = ::log10(real);

    return result;
  }

  inline double modulus(double real1, double real2) {
    double result = 0;

    if (COSNaN::is_nan(real1) || COSNaN::is_nan(real2)) {
      COSNaN::set_nan(&result);
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
}

#endif

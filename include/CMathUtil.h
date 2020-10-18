#ifndef CMathUtil_H
#define CMathUtil_H

#include <COSNaN.h>
#include <cmath>

namespace CMathUtil {

inline bool isNaN(double r) { return COSNaN::is_nan(r) && ! COSNaN::is_inf(r); }
inline bool isInf(double r) { return COSNaN::is_inf(r); }

inline double getNaN() { double r; COSNaN::set_nan(r); return r; }

inline double getNegInf() { double r; COSNaN::set_neg_inf(r); return r; }
inline double getPosInf() { double r; COSNaN::set_pos_inf(r); return r; }

inline bool isInteger(double r) {
  if (isNaN(r)) return false;

  return std::abs(r - int(r)) < 1E-3;
}

inline bool realEq(double r1, double r2, double tol=1E-6) {
  if (isNaN(r1) || isNaN(r2)) return false;

  return (std::abs(r1 - r2) < tol);
}

#if 0
inline bool realEq(double r1, double r2, double tol=1E-6) {
  if (r1 == r2 || std::abs(r1 - r2) < tol)
    return true;

  return (std::abs((r1 - r2)/(std::abs(r2) > std::abs(r1) ? r2 : r1)) <= tol);
}
#endif

inline bool isZero(double r) {
  if (isNaN(r)) return false;

  return realEq(r, 0.0);
}

//------

// sign of value
template<typename T>
int sign(T v) {
  return (T(0) < v) - (v < T(0));
}

// average of two reals
inline double avg(double x1, double x2) {
  return (x1 + x2)/2.0;
}

// map x in low->high to 0->1
inline double norm(double x, double low, double high) {
  if (high != low)
    return (x - low)/(high - low);
  else
    return low;
}

// map x in 0->1 to low->high
inline double lerp(double x, double low, double high) {
  return low + (high - low)*x;
}

// map value in range low1->high2 to low2->high2
inline double map(double value, double low1, double high1, double low2, double high2) {
  return lerp(norm(value, low1, high1), low2, high2);
}

// clamp value to range
template<typename T>
inline T clamp(const T &val, const T &low, const T &high) {
  if (val < low ) return low;
  if (val > high) return high;
  return val;
}

//------

inline double Deg2Rad(double d) { return M_PI*d/180.0; }
inline double Rad2Deg(double r) { return 180.0*r/M_PI; }

inline double normalizeAngle(double a, bool isEnd=false) {
  while (a < 0.0) a += 360.0;

  if (! isEnd) {
    while (a >= 360.0) a -= 360.0;
  }
  else {
    while (a > 360.0) a -= 360.0;
  }

  return a;
}

}

#endif

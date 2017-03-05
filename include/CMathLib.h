#ifndef CMATH_LIB_H
#define CMATH_LIB_H

namespace CMath {

// test if type is real
template<typename T> struct is_real { enum { Value=-1 }; };

template<> struct is_real<double> { enum { Value=1 }; };
template<> struct is_real<float>  { enum { Value=1 }; };

#define ASSERT_IS_REAL(T) { int a[is_real<T>::Value]; }

//---

template<typename T>
int sign(T v) {
  return (T(0) < v) - (v < T(0));
}

template<typename T>
T abs(T v) {
  return (v >= T(0) ? v : -v);
}

template<typename T>
T avg(T v1, T v2) {
  return (v1 + v2)/T(2);
}

template<typename T>
int round(T v) {
  ASSERT_IS_REAL(T)

  return (v >= T(0) ? int(v+0.5) : int(v-0.5));
}

template<typename T>
bool realEq(T r1, T r2, T tol=1E-6) {
  ASSERT_IS_REAL(T)

  if (r1 == r2 || abs(r1 - r2) < tol) return true;

  return (abs((r1 - r2)/(abs(r2) > abs(r1) ? r2 : r1)) <= tol);
}

template<typename T>
double radToDeg(T v) {
  return (double(v)*180.0/M_PI);
}

template<typename T>
double degToRad(T v) {
  return (double(v)/180.0*M_PI);
}

#undef ASSERT_IS_REAL

}

#endif

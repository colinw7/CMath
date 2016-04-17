#ifndef CSPHERICAL_3D_H
#define CSPHERICAL_3D_H

#include <CPoint3D.h>

template<typename T>
class CSpherical3DT {
 private:
  typedef CPoint3DT<T> Point;

 private:
  T rho_;    // radius
  T theta_;  // azimuth
  T phi_;    // zenith

 public:
  CSpherical3DT(T rho, T theta, T phi) :
   rho_(rho), theta_(theta), phi_(phi) {
  }

  CSpherical3DT(const Point &point) {
    rho_   = ::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    theta_ = ::atan(point.y/point.x);

    T r = ::sqrt(point.x*point.x + point.y*point.y);

    phi_ = ::asin(r/rho_);
  }

  T getRho  () const { return rho_  ; }
  T getTheta() const { return theta_; }
  T getPhi  () const { return phi_  ; }

  void setRho  (T rho  ) { rho_   = rho  ; }
  void setTheta(T theta) { theta_ = theta; }
  void setPhi  (T phi  ) { phi_   = phi  ; }

  Point toPoint() const {
    T r = rho_*::sin(phi_);

    Point(r   *::cos(theta_),
          r   *::sin(theta_),
          rho_*::cos(phi_));
  }

  void toXYZ(T *x, T *y, T *z) {
    T r = rho_*::sin(phi_);

    *x = r   *::cos(theta_);
    *y = r   *::sin(theta_);
    *z = rho_*::cos(phi_);
  }

  static void xyzToRhoThetaPhi(T x, T y, T z,
                               T *rho, T *theta, T *phi) {
    *rho   = ::sqrt(x*x + y*y + z*z);
    *theta = ::atan(y/x);

    T r = ::sqrt(x*x + y*y);

    *phi = ::asin(r/(*rho_));
  }
};

typedef CSpherical3DT<double> CSpherical3D;

#endif

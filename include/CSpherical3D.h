#ifndef CSPHERICAL_3D_H
#define CSPHERICAL_3D_H

#include <CPoint3D.h>

class CSpherical3D {
 private:
  typedef CPoint3D Point;

 public:
  CSpherical3D(double rho, double theta, double phi) :
   rho_(rho), theta_(theta), phi_(phi) {
  }

  CSpherical3D(const Point &point) {
    rho_   = ::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    theta_ = ::atan(point.y/point.x);

    double r = ::sqrt(point.x*point.x + point.y*point.y);

    phi_ = ::asin(r/rho_);
  }

  double getRho  () const { return rho_  ; }
  double getTheta() const { return theta_; }
  double getPhi  () const { return phi_  ; }

  void setRho  (double rho  ) { rho_   = rho  ; }
  void setTheta(double theta) { theta_ = theta; }
  void setPhi  (double phi  ) { phi_   = phi  ; }

  Point toPoint() const {
    double r = rho_*::sin(phi_);

    Point(r   *::cos(theta_),
          r   *::sin(theta_),
          rho_*::cos(phi_));
  }

  void toXYZ(double *x, double *y, double *z) {
    double r = rho_*::sin(phi_);

    *x = r   *::cos(theta_);
    *y = r   *::sin(theta_);
    *z = rho_*::cos(phi_);
  }

  static void xyzToRhoThetaPhi(double x, double y, double z,
                               double *rho, double *theta, double *phi) {
    *rho   = ::sqrt(x*x + y*y + z*z);
    *theta = ::atan(y/x);

    double r = ::sqrt(x*x + y*y);

    *phi = ::asin(r/(*rho_));
  }

 private:
  double rho_;    // radius
  double theta_;  // azimuth
  double phi_;    // zenith
};

#endif

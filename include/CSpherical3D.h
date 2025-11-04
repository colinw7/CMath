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
    rho_   = std::sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    theta_ = std::atan(point.y/point.x);

    double r = std::sqrt(point.x*point.x + point.y*point.y);

    phi_ = std::asin(r/rho_);
  }

  double getRho  () const { return rho_  ; }
  double getTheta() const { return theta_; }
  double getPhi  () const { return phi_  ; }

  void setRho  (double rho  ) { rho_   = rho  ; }
  void setTheta(double theta) { theta_ = theta; }
  void setPhi  (double phi  ) { phi_   = phi  ; }

  Point toPoint() const {
    double r = rho_*::sin(phi_);

    Point(r   *std::cos(theta_),
          r   *std::sin(theta_),
          rho_*std::cos(phi_));
  }

  void toXYZ(double *x, double *y, double *z) {
    double r = rho_*::sin(phi_);

    *x = r   *std::cos(theta_);
    *y = r   *std::sin(theta_);
    *z = rho_*std::cos(phi_);
  }

  static void xyzToRhoThetaPhi(double x, double y, double z,
                               double *rho, double *theta, double *phi) {
    *rho   = std::sqrt(x*x + y*y + z*z);
    *theta = std::atan(y/x);

    double r = std::sqrt(x*x + y*y);

    *phi = std::asin(r/(*rho_));
  }

 private:
  double rho_;   // radius
  double theta_; // azimuth
  double phi_;   // zenith
};

#endif

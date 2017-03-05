#ifndef CCYLINDRICAL_3D_H
#define CCYLINDRICAL_3D_H

// Cylindrical Coordinate System

#include <CPoint3D.h>

class CCylindrical3D {
 public:
  CCylindrical3D(double r, double theta, double z) :
   r_(r), theta_(theta), z_(z) {
  }

  CCylindrical3D(const CPoint3D &point) :
   r_    (::sqrt((point.x*point.x) + (point.y*point.y))),
   theta_(::atan(point.y/point.x)),
   z_    (rect->z) {
  }

  double getR    () const { return r_    ; }
  double getTheta() const { return theta_; }
  double getZ    () const { return z_    ; }

  void setR    (double r    ) { r_     = r    ; }
  void setTheta(double theta) { theta_ = theta; }
  void setZ    (double z    ) { z_     = z    ; }

  CPoint3D toPoint() const {
    return CPoint3D(r_*::cos(theta_), r_*::sin(theta_), z_);
  }

  void toXYZ(double *x, double *y, double *z) const {
    *x = r_*::cos(theta_);
    *y = r_*::sin(theta_);
    *z = z_;
  }

  static void xyzToRThetaZ(double x, double y, double z, double *r1, double *theta1, double *z1) {
    *r1     = ::sqrt((x*x) + (y*y));
    *theta1 = ::atan(y/x);
    *z1     = z;
  }

 private:
  double r_ { 1 };
  double theta_ { 0 };
  double z_ { 0 };
};

#endif

#ifndef CCYLINDRICAL_3D_H
#define CCYLINDRICAL_3D_H

// Cylindrical Coordinate System

#include <CPoint3D.h>

template<typename T>
class CCylindrical3DT {
 private:
  T r_;
  T theta_;
  T z_;

 public:
  CCylindrical3DT(T r, T theta, T z) :
   r_(r), theta_(theta), z_(z) {
  }

  CCylindrical3DT(const CPoint3DT<T> &point) :
   r_    (::sqrt((point.x*point.x) + (point.y*point.y))),
   theta_(::atan(point.y/point.x)),
   z_    (rect->z) {
  }

  T getR    () const { return r_    ; }
  T getTheta() const { return theta_; }
  T getZ    () const { return z_    ; }

  void setR    (T r    ) { r_     = r    ; }
  void setTheta(T theta) { theta_ = theta; }
  void setZ    (T z    ) { z_     = z    ; }

  CPoint3DT<T> toPoint() const {
    return CPoint3DT<T>(r_*::cos(theta_),
                        r_*::sin(theta_),
                        z_);
  }

  void toXYZ(T *x, T *y, T *z) const {
    *x = r_*::cos(theta_);
    *y = r_*::sin(theta_);
    *z = z_;
  }

  static void xyzToRThetaZ(T x, T y, T z, T *r1, T *theta1, T *z1) {
    *r1     = ::sqrt((x*x) + (y*y));
    *theta1 = ::atan(y/x);
    *z1     = z;
  }
};

typedef CCylindrical3DT<double> CCylindrical3D;

#endif

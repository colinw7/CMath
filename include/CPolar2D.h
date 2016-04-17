#ifndef CPOLAR_2D_H
#define CPOLAR_2D_H

#include <CPoint2D.h>

template<typename T>
class CPolar2DT {
 private:
  T r_;
  T theta_;

 public:
  CPolar2DT(T r, T theta) :
   r_(r), theta_(theta) {
  }

  CPolar2DT(const CPoint2DT<T> &point) :
   r_    (::sqrt(point.x*point.x + point.y*point.y)),
   theta_(::atan2(point.y, point.x)) {
  }

  T getR    () const { return r_    ; }
  T getTheta() const { return theta_; }

  void setR    (T r    ) { r_     = r    ; }
  void setTheta(T theta) { theta_ = theta; }

  CPoint2D toPoint() const {
    return CPoint2D(r_*::cos(theta_),
                    r_*::sin(theta_));
  }

  void toXY(T *x, T *y) const {
    *x = r_*::cos(theta_);
    *y = r_*::sin(theta_);
  }

  static void xyToRTheta(T x, T y, T *r, T *theta) {
    *r     = ::sqrt(x*x + y*y);
    *theta = ::atan2(y, x);
  }
};

#endif

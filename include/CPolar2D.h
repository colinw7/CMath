#ifndef CPOLAR_2D_H
#define CPOLAR_2D_H

#include <CPoint2D.h>

class CPolar2D {
 public:
  CPolar2D(double r, double theta) :
   r_(r), theta_(theta) {
  }

  CPolar2D(const CPoint2D &point) :
   r_    (::sqrt(point.x*point.x + point.y*point.y)),
   theta_(::atan2(point.y, point.x)) {
  }

  double getR    () const { return r_    ; }
  double getTheta() const { return theta_; }

  void setR    (double r    ) { r_     = r    ; }
  void setTheta(double theta) { theta_ = theta; }

  CPoint2D toPoint() const {
    return CPoint2D(r_*::cos(theta_),
                    r_*::sin(theta_));
  }

  void toXY(double *x, double *y) const {
    *x = r_*::cos(theta_);
    *y = r_*::sin(theta_);
  }

  static void xyToRTheta(double x, double y, double *r, double *theta) {
    *r     = ::sqrt(x*x + y*y);
    *theta = ::atan2(y, x);
  }

 private:
  double r_;
  double theta_;
};

#endif

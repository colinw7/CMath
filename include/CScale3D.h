#ifndef CScale3D_H
#define CScale3D_H

#include <CMatrix3D.h>

class CScale3D {
 public:
  CScale3D() { }

  CScale3D(double x, double y, double z) :
   x_(x), y_(y), z_(z) {
  }

  const CMatrix3D &matrix() const {
    if (! matrixValid_) {
      auto *th = const_cast<CScale3D *>(this);

      th->matrix_ = CMatrix3D::scale(x_, y_, z_);

      th->matrixValid_ = true;
    }

    return matrix_;
  }

  double x() const { return x_; }
  void setX(double r) { x_ = r; matrixValid_ = false; }

  double y() const { return y_; }
  void setY(double r) { y_ = r; matrixValid_ = false; }

  double z() const { return z_; }
  void setZ(double r) { z_ = r; matrixValid_ = false; }

  CPoint3D point() const { return CPoint3D(x_, y_, z_); }

 private:
  double x_ { 1.0 };
  double y_ { 1.0 };
  double z_ { 1.0 };

  CMatrix3D matrix_;
  bool      matrixValid_ { false };
};

#endif

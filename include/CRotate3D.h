#ifndef CRotate3D_H
#define CRotate3D_H

#include <CMatrix3D.h>
#include <CQuaternion.h>
#include <CPoint4D.h>

class CRotate3D {
 public:
  CRotate3D() { }

  CRotate3D(double w, double x, double y, double z) :
   w_(w), x_(x), y_(y), z_(z) {
  }

  explicit CRotate3D(const CQuaternion &q) :
    w_(q.getW()), x_(q.getX()), y_(q.getY()), z_(q.getZ()) {
  }

  const CMatrix3D &matrix() const {
    if (! matrixValid_) {
      auto *th = const_cast<CRotate3D *>(this);

      auto q = CQuaternion(w_, x_, y_, z_);

      q.toRotationMatrix(th->matrix_);

      th->matrixValid_ = true;
    }

    return matrix_;
  }

  double w() const { return z_; }
  void setW(double r) { w_ = r; matrixValid_ = false; }

  double x() const { return x_; }
  void setX(double r) { x_ = r; matrixValid_ = false; }

  double y() const { return y_; }
  void setY(double r) { y_ = r; matrixValid_ = false; }

  double z() const { return z_; }
  void setZ(double r) { z_ = r; matrixValid_ = false; }

  CPoint4D point() const { return CPoint4D(w_, x_, y_, z_); }

 private:
  double w_ { 0.0 };
  double x_ { 0.0 };
  double y_ { 0.0 };
  double z_ { 0.0 };

  CMatrix3D matrix_;
  bool      matrixValid_ { false };
};

#endif
